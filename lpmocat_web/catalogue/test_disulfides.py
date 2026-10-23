import json

from django.test import SimpleTestCase

from .makePlots import makePlots
from .views import _with_disulfide_conect


def _atom(serial, name, resname, number):
    return (f'ATOM  {serial:5d} {name:<4} {resname:>3} A{number:4d}      '
            f'0.000   0.000   0.000  1.00  0.00           {name[0]}')


PDB = '\n'.join([
    _atom(1, 'N', 'HIS', 1),
    _atom(2, 'SG', 'CYS', 20),
    _atom(3, 'SG', 'CYS', 45),
    _atom(10732, 'SG', 'CYS', 90),
    'HETATM11000 CU    CU A 201       0.000   0.000   0.000  1.00  0.00          CU',
    'CONECT    111000',
    'CONECT11000    1',
    'END',
    '',
])


class DisulfideConectTests(SimpleTestCase):
    def test_bonds_merged_with_own_records_both_ways_sorted_before_end(self):
        lines = _with_disulfide_conect(PDB, json.dumps([[20, 45, 2.72], [45, 90, 2.05]])).split('\n')
        conect = [line for line in lines if line.startswith('CONECT')]
        self.assertEqual(conect, ['CONECT    111000',
                                  'CONECT    2    3',
                                  'CONECT    3    210732',
                                  'CONECT10732    3',
                                  'CONECT11000    1'])
        self.assertEqual(lines[-2:], ['END', ''])
        self.assertEqual(lines.index('END'), lines.index('CONECT11000    1') + 1)

    def test_unchanged_without_bonds_or_matching_cysteines(self):
        self.assertEqual(_with_disulfide_conect(PDB, None), PDB)
        self.assertEqual(_with_disulfide_conect(PDB, '[]'), PDB)
        self.assertEqual(_with_disulfide_conect(PDB, json.dumps([[20, 99, 2.0]])), PDB)


class DisulfideBondsTests(SimpleTestCase):
    def _plotter(self, aligned, numbers, bonds):
        plotter = makePlots.__new__(makePlots)
        plotter.oneLetterSequence = [list(aligned)]
        plotter.residueNumbers = [numbers]
        plotter.disulfideList = [bonds]
        return plotter

    def test_columns_and_positions_in_structure_numbering(self):
        # numbered from 25; the Cys are residues 27 and 29 = positions 3 and 5
        plotter = self._plotter('HG-C-AC', [25, 26, 27, 28, 29], [[27, 29, 2.03]])
        self.assertEqual(plotter.disulfideBonds(0), [{
            'a': 27, 'b': 29, 'distance': 2.03, 'group': 'Cys27–Cys29',
            'columns': (4, 7), 'positions': (3, 5)}])

    def test_bond_with_a_masked_cysteine_is_left_out(self):
        # the second Cys lies past the catalytic domain, masked as a gap
        plotter = self._plotter('HGC---', [1, 2, 3], [[3, 6, 2.03]])
        self.assertEqual(plotter.disulfideBonds(0), [])
