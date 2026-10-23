"""The binding surface of pipeline step 7c: a region holds several short runs (the plane
touches scattered residues), where the former templates gave one range per region."""
from types import SimpleNamespace

from django.test import SimpleTestCase

from .views import _all_ranges, _loop_flexibility

# SmAA10A as step 7c finds it: 17 residues over five loops plus the N-His under "other"
SURFACE = {'other': ['1-1'],
           'L2': ['24-24', '27-28', '30-31', '33-33', '55-55'],
           'L3': ['84-88'], 'LS': ['126-126'], 'L8': ['155-156'], 'β9': ['158-158']}
LOOPS = {'L2': ['7-69'], 'L3': ['84-88'], 'LS': ['122-129'], 'L8': ['154-157']}


def _protein(surface=SURFACE, loops=LOOPS, rmsf=None):
    flexibility = SimpleNamespace(
        all=lambda: [SimpleNamespace(method='cabsflex', rmsf=rmsf)] if rmsf else [])
    return SimpleNamespace(binding_surface=str(surface), loop_regions=str(loops),
                           flexibility=flexibility)


class AllRangesTest(SimpleTestCase):
    def test_every_run_is_read(self):
        self.assertEqual(_all_ranges(SURFACE['L2']),
                         [(24, 24), (27, 28), (30, 31), (33, 33), (55, 55)])

    def test_bad_input(self):
        self.assertEqual(_all_ranges(None), [])
        self.assertEqual(_all_ranges(['not-a-range', '7-9']), [(7, 9)])


class LoopOverlapTest(SimpleTestCase):
    def test_all_runs_of_a_loop_are_kept(self):
        """The old code took the first run only, so L2 looked like a single residue."""
        _, overlap, _ = _loop_flexibility(_protein())
        self.assertEqual(overlap['L2'], '24-24,27-28,30-31,33-33,55-55')
        self.assertEqual(overlap['L3'], '84-88')

    def test_runs_are_clipped_to_the_loop(self):
        surface = dict(SURFACE, L8=['150-156'])          # starts before the loop (154-157)
        _, overlap, _ = _loop_flexibility(_protein(surface=surface))
        self.assertEqual(overlap['L8'], '154-156')

    def test_a_loop_without_surface_is_absent(self):
        _, overlap, _ = _loop_flexibility(_protein(surface={'L3': ['84-88']}))
        self.assertNotIn('L2', overlap)

    def test_mean_flexibility_covers_every_run(self):
        # RMSF 1.0 at positions 24-33, 0.0 elsewhere: L2's surface is 24, 27, 28, 30, 31, 33
        # (all 1.0) and 55 (0.0) -> 6/7
        values = [0.0] * 200
        for position in range(24, 34):
            values[position - 1] = 1.0
        _, _, averages = _loop_flexibility(_protein(rmsf='{"values": %s}' % values))
        self.assertEqual(averages['L2'], round(6 / 7, 2))
