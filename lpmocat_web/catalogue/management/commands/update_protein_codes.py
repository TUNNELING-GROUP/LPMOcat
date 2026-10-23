"""
update_protein_codes.py — Assign human-readable LPMO abbreviations as protein_code.

Reads pdb_selection.csv (which has exactly one row per DB protein), extracts
the standard LPMO name from the protein_name field (e.g. 'LsAA9A', 'CjLPMO10B'),
or constructs one from organism + family when no standard name is present.

When a species has at least one protein with an extracted prefix, that prefix
is inherited by proteins of the same species that have no extractable name
(e.g. Thermothelomyces thermophilus uses 'Mt' from MtAA13/MtLPMO9A, not 'Tt').
A name without an organism prefix, as Swiss-Prot gives many ('LPMO9C'), takes
the species' prefix (MtAA9C) unless that code is taken in either form.

A code in lpmocat_pipeline/curation/protein_codes.csv (struct_id, code, reason) wins
over both: it is the name the literature uses where CAZy lists another one first
(AF-A0A2I6QAZ5: PcAA14A, not PcoAA14A).

Regioselectivity is not set here: populate_db derives it from the EC numbers.

After updating protein_code, all related prot_code fields in linked models are
updated to match.

Usage:
  python manage.py update_protein_codes [--dry-run] [--proteins STRUCT_ID ...]
"""

import csv
import re
from collections import defaultdict
from pathlib import Path

from django.core.management.base import BaseCommand
from django.db import transaction

from catalogue.models import (
    Electrostatics, Flexibility, HoleHopping,
    Protein, SurfaceAccesability,
)
from catalogue.overview import species_name
from catalogue.management.commands.rename_lpmo_to_aa import CONFLICT_OVERRIDES

PROJECT_ROOT = Path(__file__).parents[4]
PDB_SEL_CSV  = PROJECT_ROOT / "data" / "02_data" / "pdb_selection.csv"
# curated codes (struct_id, code, reason): the name the literature uses, when it is not the
# first name CAZy lists
CODES_CSV    = PROJECT_ROOT / "lpmocat_pipeline" / "curation" / "protein_codes.csv"


def curated_codes(path=CODES_CSV):
    """{struct_id: code} of the curated codes; {} without the file."""
    if not Path(path).exists():
        return {}
    with open(path, newline='') as fh:
        return {row['struct_id'].strip(): row['code'].strip()
                for row in csv.DictReader(fh) if row.get('struct_id') and row.get('code')}

# ── Regex patterns for abbreviation extraction ────────────────────────────────

STRICT_RE     = re.compile(r'^([A-Z][a-z]{1,4})(AA|LPMO)(\d{1,2})([A-Za-z]?)$')
HYPH_DIGIT_RE = re.compile(r'^([A-Z][a-z]{1,4})(AA|LPMO)(\d{1,2})-(\d)$')
BROAD_RE      = re.compile(r'^([A-Z][a-z]{1,4})(AA|LPMO|Lpmo)(\d{1,2})([_-]?[A-Za-z]?)', re.I)
PREFIX_RE     = re.compile(r'^([A-Z][a-z]{1,4})(AA|LPMO)', re.I)
BARE_RE       = re.compile(r'^(AA|LPMO)(\d{1,2})([A-Z])$')      # 'LPMO9C', 'AA9A'

SUFFIX_LETTERS = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
SKIP_SPECIES   = {'sp.', 'cf.', 'aff.', 'nr.', 'subsp.', 'var.'}


def _tokens(name):
    """Return candidate name tokens from parenthetical groups and bare words."""
    tokens = []
    for part in re.findall(r'\(([^)]+)\)', name):
        for t in re.split(r'[;,\s]+', part):
            t = t.strip()
            if t:
                tokens.append(t)
    bare = re.sub(r'\([^)]*\)', '', name)
    for t in re.split(r'[\s;,]+', bare):
        t = t.strip()
        if t:
            tokens.append(t)
    return tokens


def extract_all_abbreviations(protein_name_full):
    """Return ordered list of all valid LPMO abbreviations found in the name.

    Tries strict match first, then hyphenated-digit, then broad case-insensitive.
    Duplicate normalised forms are suppressed.
    """
    toks  = _tokens(protein_name_full)
    found = []
    seen  = set()

    def _add(s):
        if s and s not in seen:
            found.append(s)
            seen.add(s)

    for t in toks:
        if STRICT_RE.match(t):
            _add(t)

    for t in toks:
        m = HYPH_DIGIT_RE.match(t)
        if m:
            org, lt, num, dig = m.groups()
            _add(f'{org}{lt}{num}{chr(ord("A") + int(dig) - 1)}')

    for t in toks:
        m = BROAD_RE.match(t)
        if m:
            org, lt, num, grp4 = m.group(1), m.group(2), m.group(3), m.group(4)
            lt_norm = 'AA' if lt.upper() == 'AA' else 'LPMO'
            # group 4 is ([_-]?[A-Za-z]?); strip the optional prefix character
            suf = re.sub(r'[^A-Za-z]', '', grp4).upper()[:1]
            _add(f'{org}{lt_norm}{num}{suf}')

    return found


def organism_prefix(organism, org_prefix_map=None):
    """The species' prefix from extracted names, else genus + species initials.

    org_prefix_map: species (overview.species_name) → known prefix from extracted names.
    """
    if org_prefix_map and species_name(organism) in org_prefix_map:
        return org_prefix_map[species_name(organism)]
    parts = organism.split()
    if len(parts) >= 2 and parts[1].lower() not in SKIP_SPECIES:
        return parts[0][0].upper() + parts[1][0].lower()
    # Species is "sp." etc. — use first two letters of genus
    return parts[0][0].upper() + parts[0][1].lower() if len(parts[0]) > 1 else parts[0].upper()


def construct_base(organism, family, org_prefix_map=None):
    """Construct a two-letter prefix + family (no suffix)."""
    fam_m = re.match(r'(AA|LPMO)(\d{1,2})', family)
    fam_str = fam_m.group(0) if fam_m else family
    return organism_prefix(organism, org_prefix_map) + fam_str


def bare_candidates(protein_name_full, family, prefix):
    """Codes from the names without an organism prefix ('LPMO9C' of an AA9 → prefix + 'AA9C')."""
    number = re.sub(r'\D', '', family)
    found = []
    for t in _tokens(protein_name_full):
        m = BARE_RE.match(t)
        if m and m.group(2) == number and f'{prefix}AA{number}{m.group(3)}' not in found:
            found.append(f'{prefix}AA{number}{m.group(3)}')
    return found


def _aa_form(code):
    return re.sub(r'LPMO(\d)', r'AA\1', code, flags=re.I)


def _taken(used_codes):
    """Codes a new one must avoid: those used, in their AA form too (rename_lpmo_to_aa makes
    MtLPMO9C MtAA9C), and the targets of its conflict overrides."""
    return set(used_codes) | {_aa_form(code) for code in used_codes} | set(CONFLICT_OVERRIDES.values())


def assign_codes(entries, curated):
    """{struct_id: (code, source)} for entries {struct_id: {candidates, name, org, family}}:
    curated codes first, then greedily in struct_id order an extracted name, a name without
    organism prefix with the species' prefix, else a constructed code with suffix letters."""
    # Any species that has at least one protein with an extracted LPMO name
    # tells us the preferred prefix for all proteins of that species.
    org_prefix_map = {}
    for info in entries.values():
        if info['candidates']:
            m = PREFIX_RE.match(info['candidates'][0])
            if m:
                # Take first seen prefix per species (stable iteration order)
                org_prefix_map.setdefault(species_name(info['org']), m.group(1))

    used_codes = {}        # code → struct_id that claimed it
    plan = {}
    for sid, code in curated.items():
        if sid in entries:
            used_codes[code] = sid
            plan[sid] = (code, 'curated')

    for sid in sorted(entries.keys()):
        if sid in plan:
            continue
        info = entries[sid]
        assigned = None
        source = 'extracted'

        # Try each extracted candidate first
        for cand in info['candidates']:
            if cand not in used_codes:
                assigned = cand
                break

        # Then a name without organism prefix, with the species' prefix
        if assigned is None:
            taken = _taken(used_codes)
            prefix = organism_prefix(info['org'], org_prefix_map)
            assigned = next((cand for cand in bare_candidates(info['name'], info['family'], prefix)
                             if cand not in taken), None)

        # Fall back to constructed code with suffix letters (not a target of the
        # rename_lpmo_to_aa conflict overrides)
        if assigned is None:
            source = 'constructed'
            base = construct_base(info['org'], info['family'], org_prefix_map)
            assigned = next((base + letter for letter in [''] + SUFFIX_LETTERS
                             if base + letter not in used_codes
                             and base + letter not in CONFLICT_OVERRIDES.values()), base + '_ERR')
        used_codes[assigned] = sid
        plan[sid] = (assigned, source)
    return plan


class Command(BaseCommand):
    help = "Assign human-readable LPMO abbreviations as protein_code"

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true',
                            help="Print planned changes without updating the DB")
        parser.add_argument('--proteins', nargs='+', metavar='STRUCT_ID',
                            help="Limit to these current protein_codes (struct_ids)")

    def handle(self, *args, **options):
        dry_run = options['dry_run']

        # ── 1. Load pdb_selection.csv ─────────────────────────────────────────
        if not PDB_SEL_CSV.exists():
            self.stderr.write(f"Not found: {PDB_SEL_CSV}")
            return
        with open(PDB_SEL_CSV, newline='') as fh:
            sel_rows = {row['selected_pdb']: row for row in csv.DictReader(fh)}
        self.stdout.write(f"Loaded {len(sel_rows)} rows from pdb_selection.csv")

        # ── 2. Collect proteins to process ───────────────────────────────────
        qs = Protein.objects.all()
        if options['proteins']:
            qs = qs.filter(protein_code__in=options['proteins'])
        proteins = list(qs.order_by('protein_code'))
        self.stdout.write(f"Processing {len(proteins)} proteins\n")

        # ── 3. Build per-protein candidate list ───────────────────────────────
        # struct_id → {candidates, org, family, protein}
        entries = {}
        for p in proteins:
            sid = p.protein_code
            row = sel_rows.get(sid)
            if row is None:
                self.stdout.write(self.style.WARNING(
                    f"  {sid}: no row in pdb_selection.csv — skipping"))
                continue
            entries[sid] = {
                'candidates': extract_all_abbreviations(row['protein_name']),
                'name':       row['protein_name'],
                'org':        row['organism'],
                'family':     row['family'],
                'protein':    p,
            }

        # ── 4–5. Assign codes ─────────────────────────────────────────────────
        plan = {}              # struct_id → {code, source, protein}
        for sid, (code, source) in assign_codes(entries, curated_codes()).items():
            if code.endswith('_ERR'):
                self.stdout.write(self.style.ERROR(f"  {sid}: could not find a unique code ({code})"))
            plan[sid] = {'code': code, 'source': source, 'protein': entries[sid]['protein']}

        # ── 6. Print plan ─────────────────────────────────────────────────────
        self.stdout.write(f"\n{'STRUCT_ID':<22} {'NEW CODE':<20} SOURCE")
        self.stdout.write('-' * 60)
        n_extracted = n_constructed = 0
        for sid in sorted(plan.keys()):
            info = plan[sid]
            self.stdout.write(f"  {sid:<20} {info['code']:<20} {info['source']}")
            if info['source'] == 'constructed':
                n_constructed += 1
            else:
                n_extracted += 1

        self.stdout.write(f"\nExtracted: {n_extracted}  Constructed: {n_constructed}")

        if dry_run:
            self.stdout.write(self.style.WARNING('\n[DRY RUN] DB not updated'))
            return

        # ── 7. Apply updates in a transaction ────────────────────────────────
        n_updated = 0
        with transaction.atomic():
            for sid, info in plan.items():
                old_code  = sid
                new_code  = info['code']
                p         = info['protein']

                if old_code == new_code:
                    self.stdout.write(f"  {old_code}: unchanged")
                else:
                    SurfaceAccesability.objects.filter(prot_code=old_code).update(prot_code=new_code)
                    Flexibility.objects.filter(prot_code=old_code).update(prot_code=new_code)
                    Electrostatics.objects.filter(prot_code=old_code).update(prot_code=new_code)
                    HoleHopping.objects.filter(prot_code=old_code).update(prot_code=new_code)

                    p.protein_code = new_code
                    p.save(update_fields=['protein_code'])
                    self.stdout.write(f"  {old_code}  →  {new_code}  [{info['source']}]")
                    n_updated += 1

        self.stdout.write(self.style.SUCCESS(
            f"\nDone. Updated {n_updated} protein_codes."))
