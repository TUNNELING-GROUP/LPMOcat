import sys
from pathlib import Path
from types import SimpleNamespace

from Bio.SeqUtils import seq3
from django.test import SimpleTestCase

from .brace_motifs import OTHER, brace_motifs, filter_options, motif_residue_text
from .management.commands.populate_db import brace_motifs_json
from .overview import catalogue_overview
from .second_sphere import NOT_DETERMINED

# the pipeline's step-10 module (lpmocat_pipeline/steps/ is not a package)
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lpmocat_pipeline' / 'steps'))
from brace_interactions import brace_motifs as find_brace_motifs  # noqa: E402


def _chain(sequence, start=1, breaks=()):
    """Residues numbered from `start` with backbone atoms 1.3 A apart along x, so that consecutive
    residues are bonded -- except before the numbers in `breaks`, where the chain jumps 10 A."""
    residues, x = [], 0.0
    for i, letter in enumerate(sequence):
        number = start + i
        if number in breaks:
            x += 10.0
        residues.append((number, seq3(letter).upper(), {'N': (x, 0.0, 0.0), 'CA': (x + 1.2, 0.0, 0.0),
                                                        'C': (x + 2.4, 0.0, 0.0)}))
        x += 3.7
    return residues


def _motifs(n_his, internal_his, breaks=()):
    """The step-8e result of a chain N-His + filler + Internal-His (His-2 at residue 101)."""
    sequence = n_his + 'A' * (100 - len(n_his)) + internal_his
    return find_brace_motifs(_chain(sequence, breaks=breaks), 1, 101)


class BraceMotifPipelineTests(SimpleTestCase):
    def test_n_his_motifs(self):
        cases = {'HGYVESPASR': 'H-G-X4-P-X2-R', 'HGYIAKPAPS': 'H-G-X4-P', 'HGFVQNIVID': 'H-G',
                 'HYTFPDLTYD': None}
        for sequence, motif in cases.items():
            with self.subTest(sequence):
                result = _motifs(sequence, 'HAAAA')['n_his']
                self.assertEqual((result['motif'], result['sequence']), (motif, sequence))

    def test_n_his_motif_residues(self):
        residues = _motifs('HGYVESPASR', 'HAAAA')['n_his']['residues']
        self.assertEqual([(r['role'], r['name'], r['residue'], r['offset']) for r in residues],
                         [('gly', 'GLY', 2, 1), ('pro', 'PRO', 7, 6), ('arg', 'ARG', 10, 9)])

    def test_internal_his_motifs(self):
        cases = {'HHGPV': 'H-X-G-P', 'HTIGP': 'H-X2-G-P', 'HGGMF': 'H-X-G', 'HSTTS': None}
        for sequence, motif in cases.items():
            with self.subTest(sequence):
                result = _motifs('HAAAAAAAAA', sequence)['internal_his']
                self.assertEqual((result['motif'], result['sequence']), (motif, sequence))
        residues = _motifs('HAAAAAAAAA', 'HTIGP')['internal_his']['residues']
        self.assertEqual([(r['role'], r['residue'], r['offset']) for r in residues],
                         [('gly', 104, 3), ('pro', 105, 4)])

    def test_positions_after_a_chain_break_are_missing(self):
        # the Arg of H-G-X4-P-X2-R follows a break: only H-G-X4-P is certain
        result = _motifs('HGYVESPASR', 'HHGPV', breaks={10})['n_his']
        self.assertEqual((result['motif'], result['sequence']), ('H-G-X4-P', 'HGYVESPAS'))

    def test_without_the_his(self):
        self.assertEqual(find_brace_motifs(_chain('AGYVESPASR'), 1, 5), {'n_his': None, 'internal_his': None})


class BraceMotifFieldTests(SimpleTestCase):
    site = {'brace_motifs': {
        'n_his': {'motif': 'H-G-X4-P-X2-R', 'sequence': 'HGYVESPASR', 'residues': [
            {'role': 'gly', 'name': 'GLY', 'residue': 29, 'offset': 1},
            {'role': 'pro', 'name': 'PRO', 'residue': 34, 'offset': 6},
            {'role': 'arg', 'name': 'ARG', 'residue': 37, 'offset': 9}]},
        'internal_his': {'motif': None, 'sequence': 'HSTTS', 'residues': []},
    }}

    def test_loader_and_reader(self):
        motifs = brace_motifs(SimpleNamespace(brace_motifs=brace_motifs_json(self.site)))
        self.assertEqual((motifs['n_his']['motif'], motifs['n_his']['residues']), ('H-G-X4-P-X2-R', 'H-G-Y-V'))
        self.assertEqual([r['residue'] for r in motifs['n_his']['motif_residues']], [29, 34, 37])
        self.assertEqual((motifs['internal_his']['motif'], motifs['internal_his']['residues']), (OTHER, 'H-S-T-T'))
        self.assertIsNone(brace_motifs_json({'brace_motifs': None}))

    def test_residues_of_an_insertion_variant(self):
        text = brace_motifs_json({'brace_motifs': {'internal_his': {
            'motif': 'H-X2-G-P', 'sequence': 'HTIGP', 'residues': []}}})
        motifs = brace_motifs(text)
        self.assertEqual(motifs['internal_his']['residues'], 'H-T-I-G-P')
        self.assertEqual(motifs['n_his']['motif'], NOT_DETERMINED)

    def test_missing_data(self):
        for value in (None, '', 'not json', SimpleNamespace()):
            motifs = brace_motifs(value)
            self.assertEqual({m['motif'] for m in motifs.values()}, {NOT_DETERMINED})
            self.assertEqual({m['residues'] for m in motifs.values()}, {''})

    def test_hover_text(self):
        gly, pro, arg = self.site['brace_motifs']['n_his']['residues']
        self.assertEqual(motif_residue_text('n_his', 'H-G-X4-P-X2-R', arg),
                         'N-His motif H-G-X4-P-X2-R: Arg (N-His +9)')
        self.assertEqual(motif_residue_text('internal_his', 'H-X-G-P', {**gly, 'offset': 2}),
                         'Internal-His motif H-X-G-P: Gly (Internal-His +2)')

    def test_filter_options_most_common_first_other_and_undetermined_last(self):
        values = ['H-X-G', OTHER, 'H-X-G-P', NOT_DETERMINED, '', 'H-X-G-P', OTHER, 'H-X2-G-P']
        self.assertEqual(filter_options(values),
                         [('H-X-G-P', 2), ('H-X-G', 1), ('H-X2-G-P', 1), (OTHER, 2), (NOT_DETERMINED, 1)])

    def test_overview_record(self):
        def protein(code, text):
            return SimpleNamespace(
                protein_code=code, clade=SimpleNamespace(name='AA10'),
                organism=SimpleNamespace(name='Serratia marcescens', kingdom='Bacteria', lineage=None),
                substrate=SimpleNamespace(all=lambda: []), regioselectivity=SimpleNamespace(all=lambda: []),
                associated_CBMs=None, selected_structure=code,
                cat_domain_length=None, sequence='HA', loop_regions=None, brace_motifs=text)

        data = catalogue_overview([protein('A', brace_motifs_json(self.site)), protein('B', None)],
                                  lambda p: '')
        records = {r['code']: r for r in data['proteins']}
        self.assertEqual(records['A']['brace'], {'n_his': ['H-G-X4-P-X2-R', 'H-G-Y-V'],
                                                 'internal_his': [OTHER, 'H-S-T-T']})
        self.assertEqual(records['B']['brace']['n_his'], [NOT_DETERMINED, ''])
        self.assertEqual(data['other_motif'], OTHER)
