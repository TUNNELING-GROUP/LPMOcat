"""How many proteins one comparison may hold.

The compare page loads every structure into the 3D viewer and draws every plot once per
protein, so the page grows with the selection (about 0.7 MB per protein). A comparison is
therefore limited to MAX_COMPARE_PROTEINS. One exception: a whole family. Its comparison is
always the same set, so it is kept in the canonical order below and can be built ahead on the
server (`manage.py warm_compare_cache`), which leaves only the browser's share of the work."""
from django.conf import settings

from .models import Protein

MAX_COMPARE_PROTEINS = 50
# Above this many proteins the compare page opens with every plot section closed (the Sequence
# plot is otherwise open and part of the page); a section is fetched when it is opened.
PLOTS_CLOSED_ABOVE = 30


def unique(codes):
    """The codes without repeats, in their first order."""
    return list(dict.fromkeys(c for c in codes if c))


def family_codes():
    """{family: its protein codes, sorted} -- the order a whole-family comparison is shown and cached in."""
    families = {}
    for code, family in Protein.objects.values_list('protein_code', 'clade__name'):
        families.setdefault(family or '', []).append(code)
    return {family: sorted(codes) for family, codes in families.items() if family}


def whole_family(codes):
    """The family whose proteins are exactly `codes` (in any order), else None."""
    chosen = set(codes)
    for family, members in family_codes().items():
        if chosen == set(members):
            return family
    return None


def canonical_order(codes):
    """A whole family in its fixed order (so its cache entries are found), any other set as given."""
    codes = unique(codes)
    return sorted(codes) if whole_family(codes) else codes


def too_many(codes):
    """More proteins than one comparison may hold, and not a whole family."""
    codes = unique(codes)
    return len(codes) > MAX_COMPARE_PROTEINS and not whole_family(codes)


def plots_closed(codes):
    """The comparison is large enough to open with its plot sections closed."""
    return len(unique(codes)) > plots_closed_above()


def plots_closed_above():
    """PLOTS_CLOSED_ABOVE, or settings.COMPARE_PLOTS_CLOSED_ABOVE (to try the large-comparison
    page on a few proteins)."""
    return getattr(settings, 'COMPARE_PLOTS_CLOSED_ABOVE', PLOTS_CLOSED_ABOVE)


def limit_message(n):
    return (f'A comparison holds up to {MAX_COMPARE_PROTEINS} proteins ({n} were selected). '
            'A whole family can be compared as one selection.')
