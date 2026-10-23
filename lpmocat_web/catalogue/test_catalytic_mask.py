from unittest import mock

from django.test import SimpleTestCase

from . import makePlots
from .makePlots import catalytic_masked_sequence


def _pdb(residue_numbers):
    """Minimal PDB text: one CA atom per residue with the given auth numbers."""
    return '\n'.join(f'ATOM  {i:5d}  CA  ALA A{num:4d}      0.000   0.000   0.000'
                     for i, num in enumerate(residue_numbers, 1))


class CatalyticMaskTests(SimpleTestCase):
    def setUp(self):
        patcher = mock.patch.dict(makePlots.CATALYTIC_BOUNDARIES,
                                  {'1ABC': {'end_residue': 24}}, clear=True)
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_masks_residues_past_end_in_auth_numbering(self):
        # structure numbered from 22 (like PcAA14A): residues 22-24 are catalytic
        seq = 'HGYVKL'
        masked = catalytic_masked_sequence('HG-YV--KL', seq, _pdb(range(22, 28)), '1ABC')
        self.assertEqual(masked, 'HG-Y-----')

    def test_unknown_structure_is_unchanged(self):
        self.assertEqual(catalytic_masked_sequence('HG-YV', 'HGYV', _pdb(range(1, 5)), '9XYZ'), 'HG-YV')

    def test_falls_back_to_sequence_numbering_when_pdb_does_not_match(self):
        seq = 'A' * 30
        masked = catalytic_masked_sequence(seq, seq, _pdb([1, 2]), ' 1ABC ')
        self.assertEqual(masked, 'A' * 24 + '-' * 6)
