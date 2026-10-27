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


class PropertyCellsTest(SimpleTestCase):
    def _cells(self, **values):
        from unittest import mock
        from . import views
        base = {'net_charge': -7.0, 'surface_charge': -5.2, 'pi': 4.81, 'gravy': -0.4123,
                'disulfides': ['Cys14–Cys33 (L2–β9)'], 'metal': 'Cu', 'bacterial': False,
                'sequons': [{'residue': 45, 'motif': 'N-K-T', 'rsa': 42, 'region': 'L2'}],
                'n_his_methylated': 'yes', 'experimental': True,
                'modified_residues': [{'residue': 1, 'name': 'HIC', 'parent': 'HIS'}]}
        with mock.patch.object(views, '_property_values', return_value={**base, **values}):
            return views._property_cells(None)

    def test_cells(self):
        cells = self._cells()
        self.assertEqual({key: cell['text'] for key, cell in cells.items()},
                         {'net_charge': '-7.0', 'surface_charge': '-5.2', 'pi': '4.81', 'gravy': '-0.41',
                          'disulfides': '1', 'sequons': '1', 'n_his_methylated': 'Yes', 'metal': 'Cu'})
        self.assertIn('HIC1, shown here as His', cells['n_his_methylated']['title'])
        self.assertEqual(cells['sequons']['title'], 'Asn45 (N-K-T; L2, RSA 42 %)')
        self.assertEqual(cells['disulfides']['title'], 'Cys14–Cys33 (L2–β9)')

    def test_not_computed_and_none(self):
        cells = self._cells(net_charge=None, pi=None, disulfides=[], metal='')
        self.assertEqual(cells['net_charge']['text'], '')
        self.assertEqual(cells['pi']['text'], '')
        # no bond is 0, not unknown; no metal is a dash
        self.assertEqual(cells['disulfides']['text'], '0')
        self.assertEqual(cells['metal']['text'], '–')
        self.assertEqual(self._cells(disulfides=None)['disulfides']['text'], '')

    def test_sequons_of_a_bacterial_protein_are_potential_only(self):
        cells = self._cells(bacterial=True)
        self.assertIn('bacteria do not N-glycosylate', cells['sequons']['title'])
        self.assertEqual(self._cells(sequons=[])['sequons'],
                         {'text': '0', 'title': 'No N-X-S/T sequon in the catalytic domain'})

    def test_methylation_of_a_model_and_a_plain_his(self):
        self.assertEqual(self._cells(n_his_methylated=None, experimental=False, modified_residues=[])
                         ['n_his_methylated']['text'], '–')
        self.assertEqual(self._cells(n_his_methylated='no', modified_residues=[])['n_his_methylated']['text'], 'No')


class NHisMethylationTest(SimpleTestCase):
    def _protein(self, struct_id, found, brace='(1, 86)'):
        return SimpleNamespace(selected_structure=struct_id, histidine_brace=brace,
                               modified_residues=json.dumps(found) if found is not None else None)

    def test_states(self):
        from .overview import n_his_methylation
        hic = [{'residue': 1, 'name': 'HIC', 'parent': 'HIS'}]
        self.assertEqual(n_his_methylation(self._protein('4EIR', hic)), 'yes')
        self.assertEqual(n_his_methylation(self._protein('8RRY', [])), 'no')
        # another modified residue is not the N-His
        self.assertEqual(n_his_methylation(self._protein('7NTL', [{'residue': 10, 'name': 'SNN', 'parent': 'ASN'}])), 'no')
        # a model, or nothing recorded: unknown
        self.assertIsNone(n_his_methylation(self._protein('AF-Q7SHI8', [])))
        self.assertIsNone(n_his_methylation(self._protein('4EIR', None)))
