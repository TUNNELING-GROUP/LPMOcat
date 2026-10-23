from django.test import SimpleTestCase

from .protein_mass import polypeptide_kda, residues_kda


class PolypeptideMassTests(SimpleTestCase):
    def test_single_glycine_is_the_free_amino_acid(self):
        self.assertAlmostEqual(polypeptide_kda('G'), 0.07507, places=4)

    def test_peptide_bond_releases_one_water(self):
        self.assertAlmostEqual(polypeptide_kda('GG'), 2 * 0.07507 - 0.018015, places=4)

    def test_gaps_case_and_unknown_letters_are_ignored(self):
        self.assertEqual(polypeptide_kda('g-G\nX'), polypeptide_kda('GG'))

    def test_no_residues(self):
        self.assertIsNone(polypeptide_kda('--'))
        self.assertIsNone(polypeptide_kda(None))


class ResidueSetMassTests(SimpleTestCase):
    """Residues that are not one chain (surface, loops, binding surface): no terminal water."""
    def test_is_the_sum_of_residue_masses(self):
        self.assertAlmostEqual(residues_kda('G'), 0.07507 - 0.018015, places=4)
        self.assertAlmostEqual(residues_kda('GG'), 2 * (0.07507 - 0.018015), places=4)

    def test_one_water_below_the_polypeptide(self):
        self.assertAlmostEqual(polypeptide_kda('HGDKLE') - residues_kda('HGDKLE'), 0.018015, places=5)

    def test_no_residues(self):
        self.assertIsNone(residues_kda(''))
        self.assertIsNone(residues_kda(None))
