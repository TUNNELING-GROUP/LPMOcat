import json
from types import SimpleNamespace

from django.test import SimpleTestCase

from .views import _site_cells

# SmAA10A as pipeline steps 8c-8e describe it
SITE = {'motif': 'E…N-X-F', 'short': 'E-N-F',
        'residues': [{'role': 'gln_glu', 'name': 'GLU', 'residue': 60},
                     {'role': 'middle', 'name': 'ASN', 'residue': 185},
                     {'role': 'aromatic', 'name': 'PHE', 'residue': 187}],
        'distal_axial': {'name': 'ALA', 'residue': 112, 'atom': 'CB', 'distance': 3.8, 'angle': 24}}
BRACE = {'n_his': {'motif': 'H-G-X4-P-X2-R', 'sequence': 'HGYVESPASR',
                   'residues': [{'role': 'gly', 'name': 'GLY', 'residue': 29, 'offset': 1},
                                {'role': 'pro', 'name': 'PRO', 'residue': 34, 'offset': 6},
                                {'role': 'arg', 'name': 'ARG', 'residue': 37, 'offset': 9}]},
         'internal_his': {'motif': None, 'sequence': 'HSTTS', 'residues': []}}


def _protein(site=None, brace=None):
    return SimpleNamespace(second_sphere_motif=json.dumps(site) if site is not None else None,
                           brace_motifs=json.dumps(brace) if brace is not None else None)


class SiteCellsTest(SimpleTestCase):
    def test_cells(self):
        cells = _site_cells(_protein(SITE, BRACE))
        self.assertEqual(cells['second_sphere'], {'text': 'E…N-X-F', 'title': 'Glu60 · Asn185 · Phe187'})
        self.assertEqual(cells['shielding']['text'], 'Ala112')
        self.assertIn('CB 3.8 Å from the copper', cells['shielding']['title'])
        self.assertEqual(cells['n_his']['text'], 'H-G-X4-P-X2-R')
        self.assertEqual(cells['n_his']['title'], 'H-G-Y-V · motif residues Gly29, Pro34, Arg37')
        # a His with none of the motifs after it
        self.assertEqual(cells['internal_his'], {'text': 'Other', 'title': 'H-S-T-T'})

    def test_open_position(self):
        cells = _site_cells(_protein({**SITE, 'distal_axial': {'name': None}}, BRACE))
        self.assertEqual(cells['shielding']['text'], 'None (open)')

    def test_nothing_found(self):
        cells = _site_cells(_protein())
        self.assertEqual({key: cell['text'] for key, cell in cells.items()},
                         {'second_sphere': '', 'shielding': '', 'n_his': '', 'internal_his': ''})
