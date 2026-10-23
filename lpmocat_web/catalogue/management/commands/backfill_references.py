"""Backfill structure-paper citations for each protein from RCSB.

The Protein model already has `crystal_reference` (a DOI) and
`other_crystal_references` (JSON), and the compare page already renders them in a
References column — they were just never populated (the original populate_db fetch
used the wrong JSON key, `pdbx_database_id_doi`, which RCSB spells
`pdbx_database_id_DOI`, so it always came back empty).

This command fetches, for every PDB code on a protein, the primary citation
(author, year, title, DOI, PubMed) from RCSB and stores it as a richer JSON dict in
`other_crystal_references` ({pdb: {doi, pubmed, author, n_authors, year, title,
journal, selected}}), plus the selected structure's DOI in `crystal_reference` for
backward compatibility. AlphaFold-only proteins (AF-… codes) have no structure
paper and are left empty.

Usage:
    python manage.py backfill_references                 # all not-yet-populated
    python manage.py backfill_references --force          # refetch everything
    python manage.py backfill_references --proteins 2VTC AsAA10B
    python manage.py backfill_references --dry-run --limit 5
"""
import json
import re
import time

from django.core.management.base import BaseCommand

from catalogue.models import Protein

RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/{}"


def _fetch_citation(pdb_id, timeout=15):
    """Return the primary-citation dict for a PDB entry, or None on miss/error."""
    try:
        import requests
        r = requests.get(RCSB_ENTRY.format(pdb_id), timeout=timeout)
        if r.status_code != 200:
            return None
        cit = (r.json() or {}).get("rcsb_primary_citation") or {}
        # NB: RCSB spells these with uppercase DOI / mixed-case PubMed.
        doi = cit.get("pdbx_database_id_DOI")
        pmid = cit.get("pdbx_database_id_PubMed")
        authors = cit.get("rcsb_authors") or []
        if not (doi or pmid or authors):
            return None
        # RCSB authors are "Surname, F.M." — take the surname (the last token before
        # the comma, which also salvages malformed entries where the depositor put the
        # full given names first, e.g. "Rafael Fanchini Terrasan, C." -> "Terrasan").
        first = None
        if authors:
            before_comma = authors[0].split(",")[0].strip()
            first = before_comma.split()[-1] if before_comma else None
        return {
            "doi": doi,
            "pubmed": pmid,
            "author": first,
            "n_authors": len(authors),
            "year": cit.get("year"),
            "title": cit.get("title"),
            "journal": cit.get("journal_abbrev"),
        }
    except Exception:
        return None


class Command(BaseCommand):
    help = "Backfill structure-paper citations (from RCSB) into the Protein reference fields."

    def add_arguments(self, parser):
        parser.add_argument("--proteins", nargs="+", metavar="CODE",
                            help="Only these protein codes (default: all).")
        parser.add_argument("--force", action="store_true",
                            help="Refetch even if references are already stored.")
        parser.add_argument("--dry-run", action="store_true",
                            help="Report what would change without writing to the DB.")
        parser.add_argument("--limit", type=int, default=None,
                            help="Process at most this many proteins.")
        parser.add_argument("--delay", type=float, default=0.2,
                            help="Seconds to pause between RCSB requests (default 0.2).")

    def handle(self, *args, **o):
        qs = Protein.objects.all().order_by("protein_code")
        if o["proteins"]:
            qs = qs.filter(protein_code__in=o["proteins"])
        if o["limit"]:
            qs = qs[:o["limit"]]

        cache = {}            # pdb_id -> citation dict (dedupe network calls in this run)
        updated = skipped = empty = 0

        for p in qs:
            if p.other_crystal_references and not o["force"]:
                skipped += 1
                continue

            pdbs = [x.strip() for x in re.split(r"[;,\s]+", p.pdb_code or "") if x.strip()]
            sel = (p.selected_structure or "").strip()
            # Selected structure first, then the rest (order is preserved in the JSON).
            ordered = ([sel] if sel in pdbs else []) + [x for x in pdbs if x != sel]

            refs = {}
            for pdb in ordered:
                if pdb.upper().startswith(("AF-", "AF_", "C3D-")):
                    continue
                if pdb not in cache:
                    cache[pdb] = _fetch_citation(pdb)
                    time.sleep(o["delay"])
                cit = cache[pdb]
                if cit:
                    entry = dict(cit)
                    entry["selected"] = (pdb == sel)
                    refs[pdb] = entry

            if not refs:
                empty += 1
                self.stdout.write(f"  {p.protein_code:14} — no citations  ({p.pdb_code})")
                continue

            sel_doi = (refs.get(sel) or {}).get("doi")
            if not sel_doi:
                sel_doi = next((r.get("doi") for r in refs.values() if r.get("doi")), None)

            uniq = len({(r.get("doi") or r.get("pubmed") or r.get("title")) for r in refs.values()})
            self.stdout.write(
                f"  {p.protein_code:14} — {len(refs)} structure(s), {uniq} paper(s)")

            if not o["dry_run"]:
                p.crystal_reference = sel_doi
                p.other_crystal_references = json.dumps(refs)
                p.save(update_fields=["crystal_reference", "other_crystal_references"])
            updated += 1

        verb = "would update" if o["dry_run"] else "updated"
        self.stdout.write(self.style.SUCCESS(
            f"Done. {verb}={updated}  skipped(already had refs)={skipped}  no-citation={empty}"))
