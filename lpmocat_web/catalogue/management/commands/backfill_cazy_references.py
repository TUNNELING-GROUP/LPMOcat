"""Backfill characterization-paper references from CAZy into the proteins.

The scraped CAZy table (`data/01_proteins/cazy_lpmos.csv`) has a `references`
column of PubMed links to the papers that characterized each LPMO. This command
matches those rows to catalogue proteins by UniProt accession, resolves each PubMed
ID to a citation (author, year, title, DOI) via NCBI E-utilities, and stores them as
JSON in `Protein.characterization_references`.

The entries of the other lists of pipeline step 1 (`01_proteins/other_lpmos.csv`: Swiss-Prot,
BRENDA, curation/extra_proteins.csv) have the same column, with the papers that show their
activity; they are read too. A strain variant standing in for a CAZy entry without a structure
(pdb_selection.csv `variant_of`: Swiss-Prot's Q4X066 for AfAA9A) also gets that entry's papers.

This is the natural reference for proteins that have no crystal structure (the 64
AlphaFold-only entries), and complements the RCSB structure papers for the rest —
`_format_references` merges the two and dedupes by DOI/PubMed, so a paper that is
both a structure and a characterization paper is only shown once.

Usage:
    python manage.py backfill_cazy_references
    python manage.py backfill_cazy_references --force
    python manage.py backfill_cazy_references --proteins AsAA10A TrAA9
    python manage.py backfill_cazy_references --dry-run
    python manage.py backfill_cazy_references --csv /path/to/cazy_lpmos.csv
"""
import csv
import json
import os
import re
import time

from django.conf import settings
from django.core.management.base import BaseCommand, CommandError

from catalogue.models import Protein

ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
_PMID_RE = re.compile(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)")


def _resolve_pmids(pmids, delay=0.34, timeout=20):
    """Map PMID -> {pmid, doi, author, n_authors, year, title} via NCBI esummary."""
    import requests
    out = {}
    pmids = list(dict.fromkeys(str(p) for p in pmids))   # unique, order-preserving
    for i in range(0, len(pmids), 100):                  # esummary takes many ids at once
        chunk = pmids[i:i + 100]
        try:
            r = requests.get(ESUMMARY, timeout=timeout, params={
                "db": "pubmed", "id": ",".join(chunk), "retmode": "json"})
            res = (r.json() or {}).get("result", {})
        except Exception:
            res = {}
        for uid in chunk:
            d = res.get(uid)
            if not d:
                continue
            authors = d.get("authors") or []
            first = authors[0]["name"].split()[0] if authors else None
            doi = next((a["value"] for a in d.get("articleids", [])
                        if a.get("idtype") == "doi"), None)
            year = (d.get("pubdate") or "")[:4] or None
            out[uid] = {
                "pmid": uid,
                "doi": doi,
                "author": first,
                "n_authors": len(authors),
                "year": year,
                "title": d.get("title"),
            }
        time.sleep(delay)
    return out


class Command(BaseCommand):
    help = "Backfill characterization-paper references (from CAZy PMIDs) into the proteins."

    def add_arguments(self, parser):
        parser.add_argument("--proteins", nargs="+", metavar="CODE",
                            help="Only these protein codes (default: all).")
        parser.add_argument("--force", action="store_true",
                            help="Refetch even if characterization refs are already stored.")
        parser.add_argument("--dry-run", action="store_true",
                            help="Report what would change without writing to the DB.")
        parser.add_argument("--csv", default=None,
                            help="Path to cazy_lpmos.csv (default: auto-detect under the repo).")
        parser.add_argument("--delay", type=float, default=0.34,
                            help="Seconds between NCBI requests (default 0.34, ~3/s limit).")

    def _find_csv(self, override):
        if override:
            if not os.path.exists(override):
                raise CommandError(f"--csv path not found: {override}")
            return override
        base = str(settings.BASE_DIR)
        candidates = [
            os.path.join(base, "data", "01_proteins", "cazy_lpmos.csv"),
            os.path.join(base, "..", "data", "01_proteins", "cazy_lpmos.csv"),
        ]
        for c in candidates:
            if os.path.exists(c):
                return c
        raise CommandError("cazy_lpmos.csv not found; pass --csv. Tried:\n  "
                           + "\n  ".join(candidates))

    def handle(self, *args, **o):
        csv_path = self._find_csv(o["csv"])
        self.stdout.write(f"Reading CAZy table: {csv_path}")
        with open(csv_path, encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        other_path = os.path.join(os.path.dirname(csv_path), "other_lpmos.csv")
        if os.path.exists(other_path):
            with open(other_path, encoding="utf-8") as f:
                other = list(csv.DictReader(f))
            self.stdout.write(f"Reading the other lists: {other_path} ({len(other)} entries)")
            rows += other

        # UniProt accession -> Protein, also for deleted accessions that CAZy still lists
        uni2prot = {}
        for p in Protein.objects.all():
            for u in re.split(r"[;,\s]+", p.uniprot_accessions):
                u = u.strip()
                if u:
                    uni2prot.setdefault(u.upper(), p)

        # GenBank accession -> Protein, for structures found through UniParc or CAZyme3D, whose
        # UniProt ID (if any) is not the one CAZy lists; pdb_selection.csv records the GenBank ID.
        # A duplicate entry (identical sequence, no structure of its own) points to the protein
        # holding its UniProt ID, so its characterization papers still reach that protein.
        gb2prot = {}
        sel_path = os.path.join(os.path.dirname(os.path.dirname(csv_path)),
                                "02_data", "pdb_selection.csv")
        if os.path.exists(sel_path):
            by_structure = {p.selected_structure: p for p in Protein.objects.all()}
            with open(sel_path, encoding="utf-8") as f:
                for s in csv.DictReader(f):
                    gb = (s.get("genbank_id") or "").strip().upper()
                    us = [u.strip().upper() for u in re.split(r"[;,\s]+", s.get("uniprot_ids") or "") if u.strip()]
                    prot = (by_structure.get(s.get("selected_pdb"))
                            or next((uni2prot[u] for u in us if u in uni2prot), None))
                    if gb and prot:
                        gb2prot[gb] = prot
                    variant_of = (s.get("variant_of") or "").strip().upper()
                    if variant_of and prot and s.get("selected_pdb"):
                        gb2prot.setdefault(variant_of, prot)
                        uni2prot.setdefault(variant_of, prot)

        # struct_id -> Protein, for the local models of step 1e (no accession of their own)
        by_structure = {p.selected_structure: p for p in Protein.objects.all()}

        # protein_code -> ordered unique PMIDs from all matching CAZy rows
        prot_pmids = {}
        for x in rows:
            us = [u.strip().upper() for u in re.split(r"[;,\s]+", x.get("uniprot_ids") or "") if u.strip()]
            gbs = [g.strip().upper() for g in (x.get("genbank_ids") or "").split(";") if g.strip()]
            prot = (next((uni2prot[u] for u in us if u in uni2prot), None)
                    or next((gb2prot[g] for g in gbs if g in gb2prot), None)
                    or by_structure.get((x.get("struct_id") or "").strip()))
            if not prot:
                continue
            pmids = _PMID_RE.findall(x.get("references") or "")
            if not pmids:
                continue
            bucket = prot_pmids.setdefault(prot.protein_code, [])
            for m in pmids:
                if m not in bucket:
                    bucket.append(m)

        # Which proteins to (re)fill
        qs = Protein.objects.all().order_by("protein_code")
        if o["proteins"]:
            qs = qs.filter(protein_code__in=o["proteins"])
        targets = {}
        for p in qs:
            if p.characterization_references and not o["force"]:
                continue
            if p.protein_code in prot_pmids:
                targets[p.protein_code] = prot_pmids[p.protein_code]

        all_pmids = sorted({m for lst in targets.values() for m in lst})
        self.stdout.write(f"{len(targets)} proteins to fill, {len(all_pmids)} unique PMIDs to resolve...")
        resolved = _resolve_pmids(all_pmids, delay=o["delay"]) if all_pmids else {}

        updated = empty = 0
        for code, pmids in targets.items():
            cites = [resolved[m] for m in pmids if m in resolved]
            if not cites:
                empty += 1
                self.stdout.write(f"  {code:14} — {len(pmids)} PMID(s) but none resolved")
                continue
            self.stdout.write(f"  {code:14} — {len(cites)} characterization ref(s)")
            if not o["dry_run"]:
                p = Protein.objects.filter(protein_code=code).first()
                if p:
                    p.characterization_references = json.dumps(cites)
                    p.save(update_fields=["characterization_references"])
            updated += 1

        verb = "would update" if o["dry_run"] else "updated"
        self.stdout.write(self.style.SUCCESS(
            f"Done. {verb}={updated}  unresolved={empty}  "
            f"(proteins with no CAZy PMID are left as-is)"))
