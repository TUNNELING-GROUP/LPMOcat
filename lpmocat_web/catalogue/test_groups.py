import json

from django.test import SimpleTestCase

from . import groups


class GroupParsingTests(SimpleTestCase):
    codes = ['A', 'B', 'C', 'D', 'E']

    def test_cleans_names_codes_and_duplicates(self):
        raw = json.dumps([{'name': '  Cellulose   only ', 'codes': ['B', 'A', 'X']},
                          {'name': '', 'codes': ['A', 'C']},
                          {'name': 'Empty', 'codes': ['X']}])
        parsed = groups.parse(raw, self.codes)
        self.assertEqual([(g['name'], g['codes']) for g in parsed],
                         [('Cellulose only', ['B', 'A']), ('Group 2', ['C'])])
        self.assertEqual([g['color'] for g in parsed], groups.GROUP_COLORS[:2])

    def test_at_most_three_groups(self):
        raw = json.dumps([{'name': n, 'codes': [n]} for n in 'ABCD'])
        self.assertEqual(len(groups.parse(raw, self.codes)), groups.MAX_GROUPS)

    def test_bad_input_gives_no_groups(self):
        for raw in (None, '', 'not json', '{"a": 1}', '[1, 2]'):
            self.assertEqual(groups.parse(raw, self.codes), [])

    def test_order_and_ungrouped(self):
        parsed = groups.parse(json.dumps([{'name': 'x', 'codes': ['D', 'B']}]), self.codes)
        self.assertEqual(groups.ordered_codes(parsed, self.codes), ['B', 'D', 'A', 'C', 'E'])
        everything = groups.with_ungrouped(parsed, self.codes)
        self.assertEqual(everything[-1]['name'], groups.UNGROUPED_NAME)
        self.assertEqual(everything[-1]['codes'], ['A', 'C', 'E'])
        self.assertEqual(json.loads(groups.to_param(everything)), [{'name': 'x', 'codes': ['D', 'B']}])
        self.assertEqual(groups.with_ungrouped([], self.codes), [])
