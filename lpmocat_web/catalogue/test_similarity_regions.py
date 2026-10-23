"""Region numbering in the Further-analysis tracks (catalogue.similarity).

loop_regions, beta_core and binding_surface hold 1-based SEQUENCE POSITIONS, while
col2auth, the histidine brace and hole hopping hold the structure's own RESIDUE NUMBERS.
49 of 79 AA9 structures do not start at residue 1, so mixing the two shifts every loop by
the offset (AgAA9A: L2 positions 17-51 are residues 36-70, not 20-51).
"""
import json

from django.test import SimpleTestCase

from .similarity import _in_region, _ranges_by_name


class RangesByNameTests(SimpleTestCase):
    def test_reads_json_and_python_repr(self):
        """populate_db writes JSON; older rows hold a Python dict repr."""
        self.assertEqual(_ranges_by_name('{"L2": ["3-5"]}'), {"L2": {3, 4, 5}})
        self.assertEqual(_ranges_by_name("{'L2': ['3-5']}"), {"L2": {3, 4, 5}})

    def test_several_runs_and_empty(self):
        self.assertEqual(_ranges_by_name('{"L3": ["1-2", "7-7"]}'), {"L3": {1, 2, 7}})
        self.assertEqual(_ranges_by_name(None), {})
        self.assertEqual(_ranges_by_name(""), {})
        self.assertEqual(_ranges_by_name("not a dict"), {})

    def test_skips_a_malformed_range_without_losing_the_rest(self):
        self.assertEqual(_ranges_by_name('{"L2": ["3-5", "oops"]}'), {"L2": {3, 4, 5}})


def _protein(first_residue, n=6):
    """A stub of _extract's per-protein entry: sequence positions 1..n at columns 0..n-1,
    numbered from `first_residue` in the structure. L2 = positions 2-3."""
    return {
        "col2auth": {c: first_residue + c for c in range(n)},
        "col2pos": {c: c + 1 for c in range(n)},
        "cat_end": first_residue + n,
        "loop_pos": {"L2": {2, 3}},
        "core_pos": {"b1": {5}},
        "bind_pos": {"L2": {3}},
        "site_auth": {},
        "hh_auth": {first_residue + 4},
        "props": {"rsa": [50.0] * n},
    }


class InRegionNumberingTests(SimpleTestCase):
    def test_loop_follows_positions_when_numbering_starts_at_one(self):
        d = _protein(1)
        self.assertEqual([c for c in range(6) if _in_region(d, c, "L2")], [1, 2])

    def test_loop_follows_positions_when_numbering_is_offset(self):
        """The regression: with residue numbers 20.. the loop must stay at positions 2-3
        (columns 1-2), not move to the columns whose residue numbers are 2 and 3."""
        d = _protein(20)
        self.assertEqual([c for c in range(6) if _in_region(d, c, "L2")], [1, 2])

    def test_loops_and_core_are_complementary(self):
        d = _protein(20)
        loops = [c for c in range(6) if _in_region(d, c, "loops")]
        core = [c for c in range(6) if _in_region(d, c, "core")]
        self.assertEqual(loops, [1, 2])
        self.assertEqual(sorted(loops + core), list(range(6)))

    def test_binding_surface_region(self):
        self.assertEqual([c for c in range(6) if _in_region(_protein(20), c, "binding")], [2])

    def test_hole_hopping_still_uses_residue_numbers(self):
        """hh_auth comes from the structure, so it is matched against col2auth."""
        self.assertEqual([c for c in range(6) if _in_region(_protein(20), c, "hh")], [4])

    def test_past_the_catalytic_end_is_excluded(self):
        d = _protein(20)
        d["cat_end"] = 21                      # residues 20, 21 only
        self.assertEqual([c for c in range(6) if _in_region(d, c, "whole")], [0, 1])
