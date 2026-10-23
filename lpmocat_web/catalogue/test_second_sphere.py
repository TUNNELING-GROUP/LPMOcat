import json
import sys
from pathlib import Path
from types import SimpleNamespace

from django.test import SimpleTestCase

from .management.commands.populate_db import second_sphere_json
from .overview import catalogue_overview
from .second_sphere import NOT_DETERMINED, OPEN_POSITION, filter_options, motif_residue_text, second_sphere

# the pipeline's step-8 module and pipeline.py (lpmocat_pipeline/ is not a package)
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lpmocat_pipeline' / 'steps'))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lpmocat_pipeline'))
from brace_interactions import distal_axial_residue, second_sphere_motif  # noqa: E402
from pipeline import _copper_position, _pseudo_cu  # noqa: E402

# A copper at the origin held in the z = 0 plane: His-1 (1) N and ND1, His-2 (30) NE2; the axial
# aromatic lies below the plane (z < 0), so the solvent-facing axial position is at z > 0.
HIS1 = (1, 'HIS', {'N': (2.0, 0.0, 0.0), 'ND1': (0.0, 2.0, 0.0), 'NE2': (0.0, 4.0, 0.5)})
HIS2 = (30, 'HIS', {'ND1': (-4.0, 0.0, 0.5), 'NE2': (-2.0, 0.0, 0.0)})


def _chain(names, start=2, **atoms):
    """Residues numbered from `start` with the given names; atoms={number: {name: xyz}}."""
    return [(start + i, name, dict(atoms.get(str(start + i), {})) or {'CA': (50.0 + i, 50.0, 50.0)})
            for i, name in enumerate(names)]


def _residues(names, **atoms):
    return [HIS1] + _chain(names, **atoms) + [HIS2]


class SecondSphereMotifTests(SimpleTestCase):
    def test_his_stacked_on_his2_holds_the_gln(self):
        # H at 3, eight residues, Q at 12, X, Y at 14
        residues = _residues(['ALA', 'HIS'] + ['GLY'] * 8 + ['GLN', 'SER', 'TYR'])
        index = next(i for i, r in enumerate(residues) if r[0] == 14)
        motif = second_sphere_motif(residues, index, stacked_his={3}, cu=(0, 0, 0))
        self.assertEqual(motif['motif'], 'H-X8-Q-X-Y')
        self.assertEqual(motif['short'], 'H-Q-Y')
        self.assertEqual([(r['role'], r['residue']) for r in motif['residues']],
                         [('first', 3), ('gln_glu', 12), ('aromatic', 14)])
        self.assertEqual(motif['residues'][0]['found_by'], 'stacked on His-2')

    def test_published_motif_needs_its_aromatic(self):
        # R-X4-E-X-F (AA10); the same Arg before an E-X-Y with an Asn (AA11) is N-X-E-X-Y
        aa10 = _residues(['ARG', 'GLY', 'GLY', 'GLY', 'GLY', 'GLU', 'SER', 'PHE'])
        motif = second_sphere_motif(aa10, len(aa10) - 2, set(), (0, 0, 0))
        self.assertEqual(motif['motif'], 'R-X4-E-X-F')
        self.assertEqual(motif['residues'][0]['found_by'], 'sequence motif')
        aa11 = _residues(['ARG', 'GLY', 'GLY', 'ASN', 'ARG', 'GLU', 'MET', 'TYR'])
        self.assertEqual(second_sphere_motif(aa11, len(aa11) - 2, set(), (0, 0, 0))['motif'],
                         'N-X-E-X-Y')

    def test_aa9_his_without_stacking_by_the_published_motif(self):
        # MtAA9 (5UFV): His161 ... Gln167-X-Tyr169, His161 not stacked on His-2
        residues = _residues(['HIS', 'GLN', 'ALA', 'ASN', 'ASN', 'PRO', 'GLN', 'PHE', 'TYR'])
        motif = second_sphere_motif(residues, len(residues) - 2, set(), (0, 0, 0))
        self.assertEqual((motif['motif'], motif['residues'][0]['found_by']),
                         ('H-X5-Q-X-Y', 'sequence motif'))
        # an H-X9 His is too far for the published AA9 motif
        far = _residues(['HIS'] + ['GLY'] * 9 + ['GLN', 'PHE', 'TYR'])
        self.assertEqual(second_sphere_motif(far, len(far) - 2, set(), (0, 0, 0))['motif'], 'Q-X-Y')

    def test_hydrogen_bond_partner_else_the_gln_alone(self):
        names = ['LYS', 'GLY', 'GLY', 'GLN', 'SER', 'TYR']
        bonded = _residues(names, **{'2': {'NZ': (10.0, 0.0, 0.0)}, '5': {'OE1': (10.0, 3.0, 0.0)}})
        motif = second_sphere_motif(bonded, len(bonded) - 2, set(), (0, 0, 0))
        self.assertEqual(motif['motif'], 'K-X2-Q-X-Y')
        self.assertEqual(motif['residues'][0]['found_by'], 'hydrogen bond')
        apart = _residues(names, **{'2': {'NZ': (10.0, 0.0, 0.0)}, '5': {'OE1': (10.0, 5.0, 0.0)}})
        motif = second_sphere_motif(apart, len(apart) - 2, set(), (0, 0, 0))
        self.assertEqual((motif['motif'], motif['short']), ('Q-X-Y', 'Q-Y'))

    def test_without_gln_glu_the_one_nearest_the_copper_heads_it(self):
        names = ['GLU', 'GLY', 'GLY', 'ASN', 'SER', 'PHE']
        near = _residues(names, **{'2': {'OE1': (0.0, -5.0, 0.0), 'OE2': (0.0, -6.0, 0.0)}})
        motif = second_sphere_motif(near, len(near) - 2, set(), (0, 0, 0))
        self.assertEqual((motif['motif'], motif['short']), ('E…N-X-F', 'E-N-F'))
        self.assertEqual(motif['residues'][0], {'role': 'gln_glu', 'name': 'GLU', 'residue': 2,
                                                'cu_distance': 5.0})
        self.assertEqual(motif['residues'][1]['role'], 'middle')
        far = _residues(names, **{'2': {'OE1': (0.0, -9.0, 0.0)}})
        self.assertEqual(second_sphere_motif(far, len(far) - 2, set(), (0, 0, 0))['motif'], 'N-X-F')
        self.assertEqual(second_sphere_motif(near, len(near) - 2, set(), None)['motif'], 'N-X-F')

    def test_no_residues_before_the_aromatic(self):
        self.assertIsNone(second_sphere_motif(_residues(['TYR']), 1, set(), (0, 0, 0)))
        self.assertIsNone(second_sphere_motif(_residues(['TYR']), None, set(), (0, 0, 0)))


class CopperPositionTests(SimpleTestCase):
    """Step 6f's copper position, which steps 7c, 8 and 12 use."""

    def test_metal_held_by_both_rings(self):
        position, source = _copper_position((0.1, 0.0, 0.0), HIS1[2], HIS2[2], [9.0, 9.0, 9.0])
        self.assertEqual(source, 'metal')
        self.assertAlmostEqual(position[0], 0.1)

    def test_amine_not_needed(self):
        # 7VFC: the N-His ring is flipped, so the Cu sits on its NE2, far from the amine N
        his1 = {'N': (6.0, 0.0, 0.0), 'ND1': (0.0, 4.0, 0.0), 'NE2': (0.0, 1.9, 0.0)}
        self.assertEqual(_copper_position((0.0, 0.0, 0.0), his1, HIS2[2], None)[1], 'metal')

    def test_metal_off_the_rings_or_missing_gives_the_pseudo_cu(self):
        centroid = _pseudo_cu(HIS1[2], HIS2[2])          # as step 6f stores it
        self.assertAlmostEqual(centroid[0], (2.0 + 0.0 + 0.0 - 4.0 - 2.0) / 5)
        for metal in ((5.0, 5.0, 5.0), None):
            position, source = _copper_position(metal, HIS1[2], HIS2[2], centroid)
            self.assertEqual(source, 'pseudo-Cu')
            self.assertEqual(position, centroid)
        # one ring alone does not hold it
        self.assertEqual(_copper_position((0.0, 2.0, 0.0), HIS1[2], {}, centroid)[1], 'pseudo-Cu')
        self.assertEqual(_copper_position(None, {}, {}, None), (None, None))

    def test_pseudo_cu_needs_two_brace_nitrogens(self):
        self.assertIsNotNone(_pseudo_cu(HIS1[2], {}))
        self.assertIsNone(_pseudo_cu({'N': (0.0, 0.0, 0.0)}, {}))
        self.assertIsNone(_pseudo_cu({}, {}))


class DistalAxialTests(SimpleTestCase):
    aromatic = (40, 'TYR', {'OH': (0.0, 0.0, -3.0), 'CZ': (0.0, 0.0, -4.4)})

    def distal(self, *others, cu=(0.0, 0.0, 0.0)):
        return distal_axial_residue([HIS1, HIS2, self.aromatic, *others], 1, 30, 40, cu)

    def test_alanine_over_the_copper(self):
        alanine = (112, 'ALA', {'CA': (0.0, 0.5, 4.0), 'CB': (0.0, 0.3, 3.4), 'HB1': (0.0, 0.0, 2.0)})
        self.assertEqual(self.distal(alanine),
                         {'name': 'ALA', 'residue': 112, 'atom': 'CB', 'distance': 3.4, 'angle': 5})

    def test_glycine_counts_by_its_ca(self):
        glycine = (7, 'GLY', {'CA': (0.0, 0.0, 3.2), 'N': (0.0, 0.0, 2.5)})
        self.assertEqual(self.distal(glycine)['name'], 'GLY')

    def test_open_below_equatorial_far_or_backbone(self):
        others = [
            (5, 'LEU', {'CD1': (3.5, 0.0, 0.5)}),                    # nearly in the plane
            (6, 'ALA', {'CB': (0.0, 0.0, 5.0)}),                      # too far
            (8, 'SER', {'CA': (0.0, 0.0, 3.0), 'N': (0.0, 0.3, 3.0)}),  # backbone only
            (9, 'PHE', {'CZ': (0.0, 0.3, -3.5)}),                     # the aromatic's side
        ]
        self.assertEqual(self.distal(*others), {'name': None})

    def test_without_copper_or_aromatic(self):
        self.assertIsNone(self.distal(cu=None))
        self.assertIsNone(distal_axial_residue([HIS1, HIS2], 1, 30, 40, (0.0, 0.0, 0.0)))


class SecondSphereFieldTests(SimpleTestCase):
    site = {
        'second_sphere_motif': {'motif': 'H-X8-Q-X-Y', 'short': 'H-Q-Y', 'residues': [
            {'role': 'first', 'name': 'HIS', 'residue': 157, 'found_by': 'stacked on His-2'},
            {'role': 'gln_glu', 'name': 'GLN', 'residue': 166},
            {'role': 'aromatic', 'name': 'TYR', 'residue': 168}]},
        'distal_axial': {'name': 'ALA', 'residue': 112, 'atom': 'CB', 'distance': 3.4, 'angle': 5},
        'copper': 'metal',
    }

    def test_loader_and_reader(self):
        text = second_sphere_json(self.site)
        site = second_sphere(SimpleNamespace(second_sphere_motif=text))
        self.assertEqual((site['motif'], site['short'], site['distal']), ('H-X8-Q-X-Y', 'H-Q-Y', 'Ala'))
        self.assertEqual(site['distal_residue']['residue'], 112)
        self.assertEqual(len(site['residues']), 3)

    def test_open_position_and_missing_data(self):
        open_site = second_sphere(second_sphere_json({**self.site, 'distal_axial': {'name': None}}))
        self.assertEqual((open_site['distal'], open_site['distal_residue']), (OPEN_POSITION, None))
        self.assertIsNone(second_sphere_json({'second_sphere_motif': None, 'distal_axial': None}))
        for value in (None, '', 'not json', SimpleNamespace()):
            self.assertEqual(second_sphere(value), {'motif': '', 'short': '', 'residues': [],
                                                    'distal': '', 'distal_residue': None})

    def test_hover_text(self):
        first, gln = self.site['second_sphere_motif']['residues'][:2]
        self.assertEqual(motif_residue_text('H-X8-Q-X-Y', first),
                         'Second-sphere motif H-X8-Q-X-Y: holds the Gln/Glu (stacked on His-2)')
        self.assertEqual(motif_residue_text('E…N-X-F', {**gln, 'cu_distance': 5.2}),
                         'Second-sphere motif E…N-X-F: Gln/Glu near the copper (5.2 Å from the copper)')

    def test_filter_options_most_common_first_open_and_undetermined_last(self):
        values = ['E-N-F', NOT_DETERMINED, 'H-Q-Y', OPEN_POSITION, '', 'H-Q-Y', 'R-E-F', 'E-N-F', 'H-Q-Y']
        self.assertEqual(filter_options(values),
                         [('H-Q-Y', 3), ('E-N-F', 2), ('R-E-F', 1), (OPEN_POSITION, 1), (NOT_DETERMINED, 1)])

    def test_overview_record(self):
        def protein(code, text):
            return SimpleNamespace(
                protein_code=code, clade=SimpleNamespace(name='AA9'),
                organism=SimpleNamespace(name='Neurospora crassa', kingdom='Eukaryota', lineage=None),
                substrate=SimpleNamespace(all=lambda: []), regioselectivity=SimpleNamespace(all=lambda: []),
                associated_CBMs=None, selected_structure=code,
                cat_domain_length=None, sequence='HA', loop_regions=None, second_sphere_motif=text)

        data = catalogue_overview([protein('A', second_sphere_json(self.site)), protein('B', None)],
                                  lambda p: '')
        records = {r['code']: r for r in data['proteins']}
        self.assertEqual((records['A']['motif'], records['A']['distal']), (['H-Q-Y', 'H-X8-Q-X-Y'], 'Ala'))
        self.assertEqual((records['B']['motif'], records['B']['distal']),
                         ([NOT_DETERMINED, ''], NOT_DETERMINED))
        self.assertEqual((data['not_determined'], data['open_position']), (NOT_DETERMINED, OPEN_POSITION))
