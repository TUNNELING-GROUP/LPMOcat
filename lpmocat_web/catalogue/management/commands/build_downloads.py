"""Build the files of the Downloads page into data/django_outputs/downloads/ (so they travel in
the data bundle), and release.json, the data release every page's footer shows.

Run after a data reload, as the last command of the Django import (after pregenerate_plots):

    python manage.py build_downloads                     # release = today
    python manage.py build_downloads --release 2026-10-21
"""
import datetime
import io
import json
import re
import shutil
import subprocess
import zipfile
from pathlib import Path

from django.core.management.base import BaseCommand
from plotly.subplots import make_subplots

from catalogue import downloads, views
from catalogue.makePlots import makePlots
from catalogue.models import Protein

DATA_DIR = Path(__file__).resolve().parents[4] / 'data'


def _family(protein):
    return int(re.sub(r'\D', '', protein.clade.name if protein.clade else '') or 0)


class Command(BaseCommand):
    help = 'Build the catalogue downloads and release.json (data/django_outputs/downloads/)'

    def add_arguments(self, parser):
        parser.add_argument('--release', default=datetime.date.today().isoformat(),
                            help='the data release name shown on the site (default: today)')

    def handle(self, *args, **options):
        out_dir = downloads.DOWNLOADS_DIR
        out_dir.mkdir(parents=True, exist_ok=True)
        proteins = sorted(Protein.objects.select_related('clade', 'organism', 'electrostatics')
                          .prefetch_related('substrate', 'regioselectivity'),
                          key=lambda p: (_family(p), p.protein_code))
        proteins = [p for p in proteins if p.aligned_sequence]
        self.stdout.write(f'{len(proteins)} proteins')

        # the release first, so the READMEs written below name it
        release = self._release_info(options['release'], proteins)
        (out_dir / 'release.json').write_text(json.dumps(release, indent=1) + '\n')

        files = {}   # name in the catalogue zip -> text or bytes
        files['proteins.csv'] = downloads.csv_text(downloads.protein_rows(proteins))
        files['substrate_evidence.csv'] = downloads.csv_text(downloads.substrate_rows(proteins))
        self.stdout.write('  tables')

        # per-residue data: the compare page's plot export over the whole catalogue
        msa = {'name': [p.protein_code for p in proteins],
               'seq': [p.aligned_sequence for p in proteins]}
        pm = makePlots(proteins, msa)
        pm.makeSequencePlots(make_subplots(rows=len(proteins), cols=1), row=1, col=1)
        export = pm.makeExportData()
        scale = pm.hydrophobicityScaleNames[0]
        for view, name in downloads.VIEWS:
            files[f'residues_{name}.csv'] = downloads.residue_csv(export[view], scale)
            files[f'alignment_{name}.fasta'] = views.alignment_fasta(full=view == 'full')
        self.stdout.write('  per-residue data and alignments')

        matrix = DATA_DIR / '09_distance_matrix' / 'distance_matrix.csv'
        if matrix.exists():
            files['distance_matrix.csv'] = matrix.read_text()
        paths = {}
        for p in proteins:
            hh = getattr(p, 'hole_hopping', None)
            raw = getattr(hh, 'paths_json', None) if hh else None
            if raw:
                try:
                    paths[p.protein_code] = json.loads(raw)
                except ValueError:
                    pass
        files['hole_hopping_paths.json'] = json.dumps(paths, indent=1) + '\n'

        structures = {}
        for view, name in downloads.VIEWS:
            for filename, text in views.structure_files(proteins, catalytic=view == 'catalytic'):
                structures[f'structures/{name}/{filename}'] = text
        protonated = {f'structures_protonated/full_length/{filename}': text
                      for filename, text in views.structure_files(proteins, protonated=True)}
        self.stdout.write('  structures')

        # the single files of the page
        for name, content in files.items():
            self._write(out_dir / f'lpmocat_{name}', content)
        self._write(out_dir / 'lpmocat_structures.zip', self._zip(structures))
        self._write(out_dir / 'lpmocat_structures_protonated_full_length.zip', self._zip(protonated))

        # everything in one archive, with the README and the licence
        readme = downloads.readme_text(
            'LPMOcat catalogue', proteins, built=True,
            extra_lines=['Also: distance_matrix.csv (BLOSUM62 distances between the catalytic '
                         'domains), hole_hopping_paths.json (EHPath paths to the copper), '
                         'structures_protonated/full_length/ (whole structures with hydrogens added by PDB2PQR; '
                         'PDB and PQR files).',
                         'Data sources and their retrieval: https://lpmo-cat.polsl.pl/docs/reference/versions/'])
        everything = {'README.txt': readme, **files, **structures, **protonated}
        self._write(out_dir / 'lpmocat_catalogue.zip', self._zip(everything))

        for path in sorted(out_dir.iterdir()):
            self.stdout.write(f'  {path.name:42s} {path.stat().st_size / 1e6:7.1f} MB')
        self.stdout.write(self.style.SUCCESS(f"Release {release['release']} built in {out_dir}"))

    @staticmethod
    def _zip(members):
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, 'w', zipfile.ZIP_DEFLATED) as zf:
            for name, content in members.items():
                zf.writestr(f'lpmocat/{name}', content)
        return buffer.getvalue()

    @staticmethod
    def _write(path, content):
        tmp = path.with_suffix(path.suffix + '.tmp')
        if isinstance(content, str):
            tmp.write_text(content)
        else:
            tmp.write_bytes(content)
        shutil.move(tmp, path)

    @staticmethod
    def _release_info(release, proteins):
        """What the release holds: its name, when it was built, the proteins per family, the
        BRENDA release and the code version. The retrieval dates of the other sources are in the
        documentation (docs/reference/versions.md): file dates do not tell them reliably."""
        families = {}
        for p in proteins:
            families[p.clade.name] = families.get(p.clade.name, 0) + 1
        brenda = sorted(DATA_DIR.glob('brenda_*.json.tar.gz'))
        try:
            commit = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=DATA_DIR.parent,
                                    capture_output=True, text=True, check=True).stdout.strip()
        except (OSError, subprocess.CalledProcessError):
            commit = None
        return {
            'release': release,
            'built': datetime.datetime.now().isoformat(timespec='seconds'),
            'proteins': len(proteins),
            'families': families,
            # brenda_2026_1.json.tar.gz -> 2026.1
            'brenda_release': (brenda[-1].name.split('.json')[0].replace('brenda_', '').replace('_', '.')
                               if brenda else None),
            'code': commit,
        }
