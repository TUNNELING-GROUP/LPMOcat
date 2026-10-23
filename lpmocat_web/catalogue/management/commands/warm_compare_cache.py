"""Build the compare page of whole families ahead of time.

A whole family is the one selection above the comparison limit that the browse page accepts
(catalogue.compare_limits). This command fills its cache -- the 3D viewer colouring and every
plot section, catalytic-domain and full-length view -- so the first visitor does not wait for
the server. Run it after a deployment that emptied the plot cache and after a data reload:

    python manage.py warm_compare_cache              # the families above the limit
    python manage.py warm_compare_cache --all        # every family
    python manage.py warm_compare_cache --family AA9
"""
import time

from django.core.management.base import BaseCommand

from catalogue.compare_limits import MAX_COMPARE_PROTEINS, family_codes
from catalogue.views import _COMPARE_PLOT_TYPES, _get_molstar_coloring, _get_section_plots


class Command(BaseCommand):
    help = 'Build and cache the compare page of whole families'

    def add_arguments(self, parser):
        parser.add_argument('--family', action='append', help='a family to build (repeatable)')
        parser.add_argument('--all', action='store_true', help='every family, not only those above the limit')

    def handle(self, *args, **options):
        families = family_codes()
        if options['family']:
            unknown = [f for f in options['family'] if f not in families]
            if unknown:
                self.stderr.write(f'Unknown family: {", ".join(unknown)} (have {", ".join(sorted(families))})')
                return
            chosen = options['family']
        elif options['all']:
            chosen = sorted(families)
        else:
            chosen = sorted(f for f, codes in families.items() if len(codes) > MAX_COMPARE_PROTEINS)
        if not chosen:
            self.stdout.write(f'No family has more than {MAX_COMPARE_PROTEINS} proteins; nothing to build.')
            return

        for family in chosen:
            codes = families[family]
            start = time.time()
            self.stdout.write(f'{family}: {len(codes)} proteins')
            _get_molstar_coloring(codes)
            self.stdout.write(f'  3D viewer colouring  {time.time() - start:6.1f} s')
            for ptype in _COMPARE_PLOT_TYPES:
                for full in (False, True):
                    step = time.time()
                    _get_section_plots(codes, ptype, full=full)
                    self.stdout.write(f'  {ptype:<32} {"full length" if full else "catalytic":<11} '
                                      f'{time.time() - step:6.1f} s')
            self.stdout.write(f'  done in {time.time() - start:.0f} s')
