"""
Remove functional descriptions (regioselectivity, activity notes) from
other_names. These belong in dedicated DB fields, not in other_names.

Removed patterns:
  - Regioselectivity:  C1-hydroxylating, C4-dehydrogenating, C1/C4-hydroxylating,
                       C4-oxidizing, C1-specific, and all slash/space variants
  - Activity notes:    active on ..., high activity on ..., weak activity,
                       lignocellulose-active
  - Misc notes:        formerly ..., probable fragment, partial
  - The protein's own protein_code (CAZy lists it among the aliases)

Run:
    python manage.py clean_other_names          # dry-run
    python manage.py clean_other_names --apply  # write to DB
"""

import re
from django.core.management.base import BaseCommand
from catalogue.models import Protein


# Matches any part that is a functional/activity description rather than a name/ID.
# Each alternative is anchored to the full stripped token (case-insensitive).
_FUNCTIONAL_RE = re.compile(
    r'^('
    # Regioselectivity — any combination of C1/C4 + action word(s)
    r'C[1-4](?:[/-]\s*C?-?[1-4])?\s*[-/]?\s*(?:hydroxylating|dehydrogenating|oxidizing|specific)'
    r'(?:\s*/\s*C?-?[1-4]\s*-?\s*(?:hydroxylating|dehydrogenating))?'
    r'|'
    # "C1-hydroxylating / C4-dehydrogenating" and similar spaced variants
    r'C[1-4]-(?:hydroxylating|dehydrogenating|oxidizing|specific)'
    r'\s*/\s*C?-?[1-4]-?\s*(?:hydroxylating|dehydrogenating|oxidizing)?'
    r'|'
    # Activity notes
    r'active\s+on\b.*'
    r'|'
    r'high\s+activity\s+on\b.*'
    r'|'
    r'weak\s+activity'
    r'|'
    r'lignocellulose-active'
    r'|'
    # Misc notes
    r'probable\s+fragment'
    r'|'
    r'partial'
    r'|'
    r'formerly\b.*'
    r')$',
    re.IGNORECASE,
)


def _clean(other_names: str, protein_code: str = None) -> str:
    """Return other_names with functional descriptions and the protein's own code stripped out."""
    if not other_names or not other_names.strip():
        return other_names

    parts = [p.strip() for p in other_names.split(';') if p.strip()]
    kept = [p for p in parts if not _FUNCTIONAL_RE.match(p) and p != protein_code]
    return '; '.join(kept)


class Command(BaseCommand):
    help = "Strip regioselectivity/activity notes and the protein's own code from Protein.other_names."

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
                'DRY RUN — pass --apply to write changes.\n'
            ))

        changed = 0
        for protein in Protein.objects.exclude(other_names__isnull=True).exclude(other_names__exact=''):
            cleaned = _clean(protein.other_names, protein.protein_code)
            if cleaned != protein.other_names:
                self.stdout.write(
                    f'  {protein.protein_code}\n'
                    f'    before: {protein.other_names!r}\n'
                    f'    after:  {cleaned!r}\n'
                )
                if apply:
                    protein.other_names = cleaned
                    protein.save(update_fields=['other_names'])
                changed += 1

        verb = 'Updated' if apply else 'Would update'
        self.stdout.write(self.style.SUCCESS(f'\n{verb} {changed} proteins.'))
        if not apply:
            self.stdout.write(self.style.WARNING('Pass --apply to execute.'))
