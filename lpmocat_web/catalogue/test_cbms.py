import json
from types import SimpleNamespace

from django.test import SimpleTestCase

from .cbms import annotated_cbms, cbm_family, structure_cbms, uniprot_cbms
from .management.commands.map_domain_ranges import Command as MapDomainRanges


def _domain(source, accession, name, start, end):
    return {'source': source, 'accession': accession, 'name': name, 'start': start, 'end': end}


# ScAA10C-like entry: AA10 catalytic domain, then a CBM2 that InterPro also files under the
# CBM2/CBM3 superfamily
CBM2_ENTRY = {'domains': [
    _domain('interpro', 'IPR004302', 'Cellulose/chitin-binding protein, N-terminal', 30, 220),
    _domain('uniprot', '', 'Chitin-binding type-4', 30, 220),
    _domain('interpro', 'IPR008965', 'CBM2/CBM3, carbohydrate-binding domain superfamily', 255, 362),
    _domain('interpro', 'IPR012291', 'CBM2, carbohydrate-binding domain superfamily', 256, 360),
    _domain('interpro', 'IPR001919', 'Carbohydrate-binding type-2 domain', 258, 364),
    _domain('uniprot', '', 'CBM2', 258, 364),
]}
# CjAA10A-like entry: two chitin-binding modules that only the old names call "type-3"
CBM5_ENTRY = {'domains': [
    _domain('uniprot', '', 'Chitin-binding type-3', 256, 302),
    _domain('interpro', 'IPR003610', 'Carbohydrate binding module family 5/12', 256, 302),
    _domain('interpro', 'IPR003610', 'Carbohydrate binding module family 5/12', 342, 389),
    _domain('pfam', 'PF02839', 'Carbohydrate-binding module family 5/12', 257, 300),
]}


class CbmFamilyTests(SimpleTestCase):
    def test_family_names(self):
        self.assertEqual(cbm_family({'source': 'interpro', 'accession': 'IPR000254'}), 'CBM1')
        self.assertEqual(cbm_family({'source': 'uniprot', 'name': 'Chitin-binding type-3'}), 'CBM5/12')
        # the AA10 catalytic domain (former CBM33) and the mixed superfamily name no family
        self.assertIsNone(cbm_family({'source': 'uniprot', 'name': 'Chitin-binding type-4'}))
        self.assertIsNone(cbm_family({'source': 'interpro', 'accession': 'IPR008965'}))

    def test_one_module_per_region(self):
        self.assertEqual(annotated_cbms(CBM2_ENTRY), [{
            'family': 'CBM2', 'name': 'Carbohydrate-binding type-2 domain', 'accession': 'IPR001919',
            'start': 256, 'end': 364}])
        self.assertEqual([(c['family'], c['start'], c['end']) for c in annotated_cbms(CBM5_ENTRY)],
                         [('CBM5/12', 256, 302), ('CBM5/12', 342, 389)])
        self.assertEqual(annotated_cbms(None), [])

    def test_uniprot_cbms_of_all_entries(self):
        annotations = {'A': CBM5_ENTRY, 'B': CBM2_ENTRY, 'C': {'domains': []}}
        self.assertEqual(uniprot_cbms(['C', 'A', 'B', 'missing'], annotations), ['CBM2', 'CBM5/12'])

    def test_structure_cbms(self):
        ranges = {'X': [{'type': 'catalytic'}, {'type': 'CBM', 'family': 'CBM20'},
                        {'type': 'CBM', 'accession': 'IPR000254'}]}
        self.assertEqual(structure_cbms('X', ranges), ['CBM1', 'CBM20'])
        self.assertEqual(structure_cbms('Y', ranges), [])


class PlacedByAlignmentTests(SimpleTestCase):
    """map_domain_ranges places a CBM of another UniProt entry of the protein on the structure."""

    def setUp(self):
        self.command = MapDomainRanges()
        self.command.stderr = SimpleNamespace(write=lambda text: None)
        # the entry: 10 residues, a CBM1 at 21-40 (1-based), 10 more residues
        self.cbm = 'QCGGIGWSGPTNCAPGSTCS'
        self.entry = 'MKLSAVLAAV' + 'HGFVQNIVID' + self.cbm + 'KWLEHHHHHH'
        self.annotations = {'P1': {'domains': [_domain('uniprot', '', 'CBM1', 21, 40)]}}

    def place(self, structure_sequence, first=11):
        protein = SimpleNamespace(sequence=structure_sequence, selected_structure='AF-P0',
                                  uniprot_all=json.dumps({'P0': {}, 'P1': {}}))
        residues = list(range(first, first + len(structure_sequence)))
        return self.command._from_uniprot_entries(protein, residues, self.annotations, {'P1': self.entry})

    def test_module_in_the_structure(self):
        [cbm] = self.place('HGFVQNIVID' + self.cbm + 'KWLE')
        self.assertEqual((cbm['family'], cbm['start_residue'], cbm['end_residue']), ('CBM1', 21, 40))
        self.assertEqual(cbm['source'], 'UniProt P1, aligned')

    def test_module_mostly_missing(self):
        self.assertEqual(self.place('HGFVQNIVID' + self.cbm[:8]), [])

    def test_other_protein(self):
        self.assertEqual(self.place('WWWWWWWWWWPPPPPPPPPPDDDDDDDDDD'), [])
