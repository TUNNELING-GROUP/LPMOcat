import json
import random

from django.test import SimpleTestCase, TestCase
from django.urls import reverse

from .models import Clade, Electrostatics, Organism, Protein, SurfaceAccesability
from .sequence_search import (MAX_QUERY_LENGTH, N_WEAK_SHOWN, QueryError, align_search,
                              motif_search, parse_query, search)

AMINO = 'ACDEFGHIKLMNPQRSTVWY'


def _random_seq(rng, n):
    return ''.join(rng.choice(AMINO) for _ in range(n))


def _mutate(rng, seq, fraction):
    seq = list(seq)
    for i in rng.sample(range(len(seq)), int(fraction * len(seq))):
        seq[i] = rng.choice(AMINO.replace(seq[i], ''))
    return ''.join(seq)


class ParseQueryTests(SimpleTestCase):
    def test_plain_sequence(self):
        self.assertEqual(parse_query('  hgyvesv\n'), ('', 'HGYVESV'))

    def test_fasta_uses_first_record_only(self):
        text = '>sp|Q7Z9M7|LPMO test\nHGYVES\nVTNGG\n>second\nWWWWWW\n'
        self.assertEqual(parse_query(text), ('sp|Q7Z9M7|LPMO test', 'HGYVESVTNGG'))

    def test_genbank_numbering_gaps_and_stop_are_dropped(self):
        self.assertEqual(parse_query('1 hgyv-esv.tn 11 ggkw*')[1], 'HGYVESVTNGGKW')

    def test_invalid_characters(self):
        with self.assertRaisesMessage(QueryError, 'unexpected character(s): # @'):
            parse_query('HGYV@ES#V')

    def test_nucleotide_sequence_rejected(self):
        with self.assertRaisesMessage(QueryError, 'nucleotide'):
            parse_query('ATGGCGTTAACGGATCCGTAGCTAGCTAGGCT')

    def test_empty_and_too_long(self):
        with self.assertRaises(QueryError):
            parse_query('>header only\n')
        with self.assertRaises(QueryError):
            parse_query('A' * (MAX_QUERY_LENGTH + 1))


class SearchTests(SimpleTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        rng = random.Random(7)
        cls.rng = rng
        cls.targets = [(f'P{i}', _random_seq(rng, 200)) for i in range(20)]

    def test_self_hit_is_best_with_full_identity(self):
        code, seq = self.targets[3]
        hits = align_search(seq, self.targets)
        self.assertEqual(hits[0]['code'], code)
        self.assertEqual(hits[0]['identity'], 100.0)
        self.assertEqual(hits[0]['query_coverage'], 100.0)
        self.assertTrue(hits[0]['significant'])

    def test_homologue_found_despite_signal_peptide(self):
        code, seq = self.targets[5]
        query = 'MKSFTLLALAAAVSAQA' + _mutate(self.rng, seq, 0.4)
        hit = align_search(query, self.targets)[0]
        self.assertEqual(hit['code'], code)
        self.assertTrue(hit['significant'])
        self.assertGreater(hit['query_start'], 1)        # signal peptide left unaligned
        self.assertLess(hit['identity'], 100.0)

    def test_unrelated_query_returns_few_weak_matches(self):
        hits = align_search(_random_seq(random.Random(99), 200), self.targets)
        self.assertLessEqual(len(hits), N_WEAK_SHOWN)
        self.assertTrue(all(not h['significant'] for h in hits))

    def test_motif_wildcards_and_overlaps(self):
        targets = [('A', 'HGYVAHGFV'), ('B', 'HGWVHGW'), ('C', 'KKKK'), ('D', 'NNNDN')]
        hits = {h['code']: h['positions'] for h in motif_search('HGXV', targets)}
        self.assertEqual(hits, {'A': [1, 6], 'B': [1]})
        self.assertEqual(motif_search('NBN', targets)[0]['positions'], [1, 3])   # B = N, then D

    def test_dispatch_by_length(self):
        self.assertEqual(search('HGYV', self.targets)['mode'], 'motif')
        self.assertEqual(search(self.targets[0][1], self.targets)['mode'], 'alignment')
        with self.assertRaises(QueryError):
            search('HG', self.targets)


class SequenceSearchViewTests(TestCase):
    @classmethod
    def setUpTestData(cls):
        rng = random.Random(3)
        org = Organism.objects.create(name='Test organism')
        cls.seqs = {}
        for i, clade in enumerate(['AA9', 'AA10']):
            code = f'Tt{clade}A'
            cls.seqs[code] = _random_seq(rng, 150)
            Protein.objects.create(
                protein_code=code, protein_name=code, other_names=f'Alias{i}',
                organism=org, clade=Clade.objects.create(name=clade), pdb_code='',
                selected_structure='', sequence=cls.seqs[code], histidine_brace='',
                surface_accesability=SurfaceAccesability.objects.create(prot_code=code, sasa_file=''),
                electrostatics=Electrostatics.objects.create(prot_code=code, pqr_file=''))

    def _post(self, query):
        return self.client.post(reverse('sequence_search'), data=json.dumps({'query': query}),
                                content_type='application/json')

    def test_alignment_hit_carries_display_info(self):
        resp = self._post('>my query\n' + self.seqs['TtAA10A'])
        data = resp.json()
        self.assertEqual(resp.status_code, 200)
        self.assertEqual((data['mode'], data['header'], data['query_length']),
                         ('alignment', 'my query', 150))
        self.assertEqual(data['hits'][0]['code'], 'TtAA10A')
        self.assertEqual((data['hits'][0]['clade'], data['hits'][0]['other_names']),
                         ('AA10', 'Alias1'))

    def test_motif_mode(self):
        motif = self.seqs['TtAA9A'][10:16]
        data = self._post(motif).json()
        self.assertEqual((data['mode'], data['motif']), ('motif', motif))
        self.assertIn('TtAA9A', [h['code'] for h in data['hits']])

    def test_errors(self):
        resp = self._post('not a sequence!')
        self.assertEqual(resp.status_code, 400)
        self.assertFalse(resp.json()['ok'])
        self.assertEqual(self.client.get(reverse('sequence_search')).status_code, 405)
