"""Groups of compared proteins: named sets the user builds on the browse page ("Add as group")
and the compare page shows (bands in the table and the plots, 3D colouring, ...).

They travel in the compare page's address as ?groups=<JSON>: [{"name": ..., "codes": [...]}].
Every compared protein outside them belongs to "Ungrouped".

At most MAX_GROUPS named groups: their colours are the first three slots of the data-viz
reference palette, the only ones that stay apart for every pair (normal vision and
colour-blindness, `validate_palette.py --pairs all`); a fourth would clash with one of them.
Each colour is always shown with the group's name."""
import json

MAX_GROUPS = 3
MAX_NAME = 40
GROUP_COLORS = ['#eb6834', '#1baf7a', '#2a78d6']   # orange, aqua, blue
UNGROUPED_NAME = 'Ungrouped'
UNGROUPED_COLOR = '#a3a3a3'


def parse(raw, codes):
    """The named groups of `raw` (the ?groups= JSON) restricted to `codes`, cleaned: at most
    MAX_GROUPS, names trimmed (empty -> "Group N"), a protein in its first group only, empty
    groups dropped; [] when there are none. Each: {'name', 'codes', 'color'}."""
    try:
        data = json.loads(raw) if raw else []
    except (TypeError, ValueError):
        return []
    if not isinstance(data, list):
        return []
    known, taken, groups = set(codes), set(), []
    for item in data:
        if len(groups) == MAX_GROUPS:
            break
        if not isinstance(item, dict):
            continue
        members = [c for c in item.get('codes') or [] if isinstance(c, str)
                   and c in known and c not in taken]
        if not members:
            continue
        taken.update(members)
        name = ' '.join(str(item.get('name') or '').split())[:MAX_NAME]
        groups.append({'name': name or f'Group {len(groups) + 1}', 'codes': members,
                       'color': GROUP_COLORS[len(groups)]})
    return groups


def with_ungrouped(groups, codes):
    """The named groups plus "Ungrouped" (the other proteins of `codes`), if any are left."""
    if not groups:
        return []
    grouped = {c for g in groups for c in g['codes']}
    rest = [c for c in codes if c not in grouped]
    extra = [{'name': UNGROUPED_NAME, 'codes': rest, 'color': UNGROUPED_COLOR,
              'ungrouped': True}] if rest else []
    return groups + extra


def ordered_codes(groups, codes):
    """`codes` in group order: each group's proteins in their order in `codes`, then the rest."""
    position = {c: i for i, c in enumerate(codes)}
    grouped = [c for g in groups for c in sorted(g['codes'], key=position.get)]
    return grouped + [c for c in codes if c not in set(grouped)]


def to_param(groups):
    """The ?groups= value of the named groups (compact JSON)."""
    return json.dumps([{'name': g['name'], 'codes': g['codes']}
                       for g in groups if not g.get('ungrouped')], separators=(',', ':'))
