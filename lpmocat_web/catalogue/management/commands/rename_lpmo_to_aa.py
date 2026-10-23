"""
Rename proteins whose protein_code contains 'LPMO' to use the AA naming
convention (e.g. NcLPMO9A → NcAA9A) and add the old name to other_names.

Conflicts (where the target AA code already exists for a different protein)
are resolved by assigning the next available letter suffix.

Run:
    python manage.py rename_lpmo_to_aa          # dry-run (preview only)
    python manage.py rename_lpmo_to_aa --apply  # actually write to DB
"""

import re
from django.core.management.base import BaseCommand
from django.db import transaction
from catalogue.models import Protein


# Explicit conflict overrides: old_code -> resolved_new_code
CONFLICT_OVERRIDES = {
    'LmLPMO10': 'LmAA10B',   # LmAA10 already exists (5L2V)
    'MtLPMO9A': 'MtAA9D',    # MtAA9A/B/C already exist
    'NcLPMO9A': 'NcAA9B',    # NcAA9A already exists (AF-Q7SCJ5)
}


def _lpmo_to_aa(code):
    """Convert an LPMO-style code to AA-style. Returns None if no match."""
    m = re.match(r'^(.+?)LPMO(\d+)(.*)$', code, re.IGNORECASE)
    if not m:
        return None
    prefix, num, suffix = m.group(1), m.group(2), m.group(3)
    return f'{prefix}AA{num}{suffix}'


def _add_to_other_names(existing, name_to_add):
    """Append name_to_add to the other_names string if not already present."""
    if not name_to_add:
        return existing
    parts = [p.strip() for p in re.split(r'[;,]', existing) if p.strip()]
    if name_to_add in parts:
        return existing
    if existing and existing.strip():
        return existing.rstrip('; ') + '; ' + name_to_add
    return name_to_add


class Command(BaseCommand):
    help = "Rename LPMO-style protein_codes to AA-style and move old names to other_names."

    def add_arguments(self, parser):
        parser.add_argument(
            '--apply',
            action='store_true',
            help='Write changes to the database (default is dry-run).',
        )

    def handle(self, *args, **options):
        apply = options['apply']
        if not apply:
            self.stdout.write(self.style.WARNING(
                'DRY RUN — pass --apply to write changes to the database.'
            ))

        all_codes = set(Protein.objects.values_list('protein_code', flat=True))
        lpmo_proteins = list(
            Protein.objects.filter(protein_code__iregex=r'LPMO').order_by('protein_code')
        )
        self.stdout.write(f'Found {len(lpmo_proteins)} proteins with LPMO in code.\n')

        renames = []
        skipped = []

        for protein in lpmo_proteins:
            old_code = protein.protein_code

            if old_code in CONFLICT_OVERRIDES:
                new_code = CONFLICT_OVERRIDES[old_code]
                if new_code in all_codes:
                    skipped.append((old_code, new_code, 'override target already exists'))
                    continue
            else:
                new_code = _lpmo_to_aa(old_code)
                if new_code is None:
                    skipped.append((old_code, None, 'no LPMO pattern matched'))
                    continue
                if new_code in all_codes:
                    skipped.append((old_code, new_code, 'target code already exists — add to CONFLICT_OVERRIDES'))
                    continue

            new_other_names = _add_to_other_names(protein.other_names or '', old_code)
            renames.append((protein, old_code, new_code, new_other_names))

        self.stdout.write('=== RENAMES ===')
        for protein, old, new, names in renames:
            self.stdout.write(f'  {old:30s} → {new}')
            if names != (protein.other_names or ''):
                self.stdout.write(f'    other_names: {names}')

        if skipped:
            self.stdout.write('\n=== SKIPPED ===')
            for old, new, reason in skipped:
                self.stdout.write(self.style.WARNING(f'  {old} → {new}: {reason}'))

        self.stdout.write(f'\n{len(renames)} to rename, {len(skipped)} skipped.')

        if not apply:
            self.stdout.write(self.style.WARNING('\nNo changes written. Pass --apply to execute.'))
            return

        with transaction.atomic():
            for protein, old_code, new_code, new_other_names in renames:
                protein.protein_code = new_code
                protein.other_names = new_other_names
                protein.save(update_fields=['protein_code', 'other_names'])
                self.stdout.write(f'  renamed {old_code} → {new_code}')

        self.stdout.write(self.style.SUCCESS(f'\nDone. Renamed {len(renames)} proteins.'))
