from django.test import SimpleTestCase

from .dendrogram import average_linkage_tree, prune_tree


def _leaves(node):
    return [node['name']] if 'children' not in node else [x for c in node['children'] for x in _leaves(c)]


# Five proteins: A and B close, C near them, D and E a separate pair
LABELS = ['A', 'B', 'C', 'D', 'E']
MATRIX = [
    [0.0, 0.1, 0.3, 0.9, 0.9],
    [0.1, 0.0, 0.35, 0.9, 0.9],
    [0.3, 0.35, 0.0, 0.8, 0.85],
    [0.9, 0.9, 0.8, 0.0, 0.2],
    [0.9, 0.9, 0.85, 0.2, 0.0],
]


class AverageLinkageTreeTest(SimpleTestCase):
    def test_merges_and_heights(self):
        tree = average_linkage_tree(MATRIX, LABELS)
        # root: {A, B, C} + {D, E} at the mean of their six distances
        self.assertAlmostEqual(tree['distance'], (0.9 * 4 + 0.8 + 0.85) / 6)
        self.assertEqual(sorted(_leaves(tree)), LABELS)
        # optimizeTree: the smaller subtree first ({D, E}), then the lower merge first
        self.assertEqual(_leaves(tree), ['D', 'E', 'C', 'A', 'B'])
        abc = tree['children'][1]
        self.assertAlmostEqual(abc['distance'], (0.3 + 0.35) / 2)
        self.assertAlmostEqual(abc['children'][1]['distance'], 0.1)

    def test_ties_break_on_the_first_pair(self):
        # every distance equal: the JS merges the first pair it finds, A and B
        tree = average_linkage_tree([[0, 1, 1], [1, 0, 1], [1, 1, 0]], ['A', 'B', 'C'])
        pair = [c for c in tree['children'] if 'children' in c][0]
        self.assertEqual(sorted(_leaves(pair)), ['A', 'B'])

    def test_single_protein(self):
        self.assertEqual(average_linkage_tree([[0.0]], ['A']), {'name': 'A'})
        self.assertIsNone(average_linkage_tree([], []))


class PruneTreeTest(SimpleTestCase):
    def setUp(self):
        self.tree = average_linkage_tree(MATRIX, LABELS)

    def test_keeps_topology_and_merge_distances(self):
        pruned = prune_tree(self.tree, {'A', 'C', 'D'})
        self.assertEqual(_leaves(pruned), ['D', 'C', 'A'])
        # D's partner E is gone, so D hangs from the root; A joins C where {A, B} joined C
        self.assertEqual(pruned['distance'], self.tree['distance'])
        self.assertEqual(pruned['children'][0], {'name': 'D'})
        self.assertAlmostEqual(pruned['children'][1]['distance'], (0.3 + 0.35) / 2)

    def test_one_or_no_leaf(self):
        self.assertEqual(prune_tree(self.tree, {'B'}), {'name': 'B'})
        self.assertIsNone(prune_tree(self.tree, set()))
        self.assertIsNone(prune_tree(self.tree, {'unknown'}))

    def test_does_not_change_the_tree(self):
        before = repr(self.tree)
        prune_tree(self.tree, {'A', 'E'})
        self.assertEqual(repr(self.tree), before)
