from unittest import mock

from django.test import TestCase
from django.urls import reverse

from . import compare_limits
from .models import Clade, Electrostatics, Organism, Protein, SurfaceAccesability


class CompareLimitTests(TestCase):
    """A comparison holds a limited number of proteins, or one whole family."""

    @classmethod
    def setUpTestData(cls):
        organism = Organism.objects.create(name='Test organism')
        families = {'AA9': Clade.objects.create(name='AA9'), 'AA10': Clade.objects.create(name='AA10')}
        # AA9: 6 proteins, AA10: 3 -- with the limit patched to 4, AA9 is a family above it
        for family, n in (('AA9', 6), ('AA10', 3)):
            for i in range(n):
                code = f'Xx{family}{chr(70 - i)}'      # created in reverse alphabetical order
                Protein.objects.create(
                    protein_code=code, protein_name=code, sequence='HGAA', organism=organism,
                    clade=families[family], pdb_code='', selected_structure=code,
                    surface_accesability=SurfaceAccesability.objects.create(prot_code=code, sasa_file=''),
                    electrostatics=Electrostatics.objects.create(prot_code=code, pqr_file=''))
        cls.aa9 = sorted(Protein.objects.filter(clade__name='AA9').values_list('protein_code', flat=True))
        cls.aa10 = sorted(Protein.objects.filter(clade__name='AA10').values_list('protein_code', flat=True))

    def setUp(self):
        patcher = mock.patch.object(compare_limits, 'MAX_COMPARE_PROTEINS', 4)
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_small_selection_is_fine(self):
        self.assertFalse(compare_limits.too_many(self.aa9[:4]))

    def test_mixed_selection_above_the_limit_is_too_many(self):
        self.assertTrue(compare_limits.too_many(self.aa9[:3] + self.aa10[:2]))
        self.assertTrue(compare_limits.too_many(self.aa9[:5]))

    def test_repeats_do_not_count(self):
        self.assertFalse(compare_limits.too_many(self.aa9[:3] * 3))

    def test_whole_family_is_allowed_in_any_order(self):
        shuffled = list(reversed(self.aa9))
        self.assertEqual(compare_limits.whole_family(shuffled), 'AA9')
        self.assertFalse(compare_limits.too_many(shuffled))
        self.assertEqual(compare_limits.canonical_order(shuffled), self.aa9)

    def test_family_plus_one_is_not_a_family(self):
        codes = self.aa9 + self.aa10[:1]
        self.assertIsNone(compare_limits.whole_family(codes))
        self.assertTrue(compare_limits.too_many(codes))

    def test_other_selections_keep_their_order(self):
        codes = [self.aa9[2], self.aa10[0], self.aa9[0]]
        self.assertEqual(compare_limits.canonical_order(codes), codes)

    def test_compare_page_sends_too_many_back_to_the_selection(self):
        codes = self.aa9[:3] + self.aa10[:2]
        response = self.client.get(reverse('compare_proteins', args=[','.join(codes)]))
        self.assertEqual(response.status_code, 302)
        self.assertTrue(response.url.startswith(reverse('index') + '?select='))
        for code in codes:
            self.assertIn(code, response.url)

    def test_compare_page_puts_a_whole_family_in_its_fixed_order(self):
        response = self.client.get(reverse('compare_proteins', args=[','.join(reversed(self.aa9))]))
        self.assertEqual(response.status_code, 302)
        self.assertEqual(response.url, reverse('compare_proteins', args=[','.join(self.aa9)]))

    def test_data_endpoints_refuse_too_many(self):
        ids = ','.join(self.aa9[:3] + self.aa10[:2])
        for name, params in (('compare_plots', {'ids': ids, 'ptype': 'Sequence'}),
                             ('compare_analysis', {'ids': ids}),
                             ('get_pdb_structures', {'ids': ids})):
            response = self.client.get(reverse(name), params)
            self.assertEqual(response.status_code, 400, name)
            self.assertIn('up to 4 proteins', response.json()['error'], name)

    def test_whole_family_has_its_own_cache(self):
        from django.core.cache import caches
        from .views import _plots_cache
        self.assertIs(_plots_cache(self.aa9), caches['family_plots'])
        self.assertIs(_plots_cache(list(reversed(self.aa9))), caches['plots'])
        self.assertIs(_plots_cache(self.aa9[:3]), caches['plots'])

    def test_large_comparison_opens_with_plots_closed(self):
        self.assertFalse(compare_limits.plots_closed(self.aa9[:3]))
        with self.settings(COMPARE_PLOTS_CLOSED_ABOVE=2):
            self.assertTrue(compare_limits.plots_closed(self.aa9[:3]))
            self.assertFalse(compare_limits.plots_closed(self.aa9[:2]))

    def test_coloring_endpoint_checks_its_input(self):
        ids = ','.join(self.aa9[:2])
        self.assertEqual(self.client.get(reverse('compare_coloring'), {'ids': ids, 'kind': 'nonsense', 'name': 'x'}).status_code, 400)
        self.assertEqual(self.client.get(reverse('compare_coloring'), {'kind': 'hydrophobicityScales', 'name': 'x'}).status_code, 400)
        too_many = ','.join(self.aa9[:3] + self.aa10[:2])
        self.assertEqual(self.client.get(reverse('compare_coloring'),
                                         {'ids': too_many, 'kind': 'hydrophobicityScales', 'name': 'x'}).status_code, 400)
