"""Column-subset a full-alignment plot figure down to a selection's columns.

Plots are generated once per protein against the *full* MSA coordinate system
(makePlots(..., no_trim=True)). For a given selection we drop the alignment
columns that are all-gap across the selected set and renumber the rest to
1..K — exactly what per-combination generation used to do, but as a cheap
array remap instead of rebuilding every Plotly figure.

Key fact that makes this exact: a column is dropped only when *every* selected
protein gaps there, so every residue always survives. Subsetting therefore only
removes empty columns and shifts x-positions — colours, letters and the
Sequence heatmap's plot-relative colorscale are untouched.
"""
import gzip
import json
import re
from pathlib import Path

_ALIGN_POS_RE = re.compile(r'Alignment position: (\d+)')

# Where `manage.py pregenerate_plots` writes full-alignment figures (keep in sync with it).
PREGEN_DIR = Path(__file__).resolve().parents[2] / 'data' / 'django_outputs' / 'pregenerated_plots'

# Plot types that are pre-generated + subset (combination-independent). The rest are
# regenerated on demand.
PREGEN_MARKER_TYPES = ('Sequence', 'Charge', 'Hole Hopping', 'Hydrophobicity')


def pregen_type_slug(ptype, scale=None):
    """File slug for a (plot type, hydrophobicity scale). Matches pregenerate_plots.py."""
    if ptype == 'Hydrophobicity' and scale:
        return 'hydrophobicity__' + scale.lower().replace(' ', '_').replace('/', '_')
    return ptype.lower().replace(' ', '_').replace('/', '_')


def load_pregen_figure(protein_code, ptype, scale=None):
    """Load a pre-generated full-alignment figure dict, or None if absent."""
    path = PREGEN_DIR / protein_code / (pregen_type_slug(ptype, scale) + '.json.gz')
    if not path.exists():
        return None
    with gzip.open(path, 'rt', encoding='utf-8') as fh:
        return json.loads(fh.read())


def assemble_from_pregen(codes, ptype, scale, keep_indices):
    """Column-subset the pre-generated figures for a selection into per-protein plot
    JSON strings (in `codes` order), matching per-combination output.

    `keep_indices` are the full-alignment 0-indexed columns to keep (all-gap-trimmed for
    the full-length view, or the catalytic-domain columns). Returns None if any selected
    protein lacks a pre-generated file (caller should then regenerate live).
    """
    n = len(codes)
    out = []
    for i, code in enumerate(codes):
        full = load_pregen_figure(code, ptype, scale)
        if full is None:
            return None
        out.append(json.dumps(subset_figure(full, keep_indices, is_last=(i == n - 1))))
    return out


def compute_keep_indices(one_letter_sequences):
    """0-indexed full-alignment columns that survive for this selection: columns
    where at least one sequence has a residue (i.e. not all-gap)."""
    if not one_letter_sequences:
        return []
    return [i for i, col in enumerate(zip(*one_letter_sequences))
            if not all(res == '-' for res in col)]


def build_remap(keep_indices):
    """Map full-alignment column (1-indexed) -> trimmed position (1-indexed)."""
    return {i + 1: rank + 1 for rank, i in enumerate(keep_indices)}


def _rewrite_align_pos(value, old_col, remap):
    """Renumber a hover string's 'Alignment position: <full>' to the trimmed
    position. Passes non-strings through unchanged."""
    if not isinstance(value, str):
        return value
    return _ALIGN_POS_RE.sub(
        lambda m: ('Alignment position: %d' % remap[old_col]
                   if int(m.group(1)) == old_col and old_col in remap else m.group(0)),
        value)


def _subset_trace(tr, remap, keepset, renumber):
    tr = {k: (list(v) if isinstance(v, (list, tuple)) else v) for k, v in tr.items()}
    x = tr.get('x')
    if not isinstance(x, list) or not x:
        return tr

    # Domain-band polygon (Sequence): a filled rectangle whose x are the span
    # corners (start-0.5 .. end+0.5). Remap the endpoints, keep the 5-point shape.
    if tr.get('fill') == 'toself':
        xs = [v for v in x if v is not None]
        if xs:
            sc = int(round(min(xs) + 0.5))
            ec = int(round(max(xs) - 0.5))
            if sc in remap and ec in remap:
                nx0, nx1 = remap[sc] - 0.5, remap[ec] + 0.5
                tr['x'] = [nx0, nx1, nx1, nx0, nx0]
        return tr

    # Generic per-column trace: keep None separators, drop points at all-gap
    # columns (only the Sequence heatmap actually has any — data traces reference
    # residue columns, which all survive), remap the rest, renumber hover text.
    #
    # NOTE: this is exact only when trimming removes no residues (the all-gap /
    # full-length view). The catalytic-domain view removes residues, which changes
    # plot-relative encodings (Sequence heatmap colorscale, hover-label colour arrays),
    # so that view is REGENERATED live rather than subset — see views._get_section_plots.
    n = len(x)
    cols = [None if xi is None else round(xi) for xi in x]
    mask = [(c is None) or (c in keepset) for c in cols]
    for key, arr in list(tr.items()):
        if not (isinstance(arr, list) and len(arr) == n):
            continue
        out = []
        for j in range(n):
            if not mask[j]:
                continue
            v = arr[j]
            if key == 'x':
                v = None if cols[j] is None else remap[cols[j]]
            elif renumber and cols[j] is not None:
                if isinstance(v, str):
                    v = _rewrite_align_pos(v, cols[j], remap)
                elif isinstance(v, list):
                    v = [_rewrite_align_pos(e, cols[j], remap) for e in v]
            out.append(v)
        tr[key] = out
    return tr


def subset_figure(fig, keep_indices, renumber=True, is_last=None):
    """Return a column-subset copy of a full-alignment Plotly figure dict, matching
    what per-combination generation produced for these `keep_indices`.

    `fig` is a plain dict (json-loaded Plotly figure). `renumber=True` rewrites the
    hover 'Alignment position' labels to the trimmed numbering. `is_last` — whether
    this protein is the last row in the selection; controls x-axis visibility, the
    only position-dependent layout (only the bottom plot shows the axis). Leave None
    to keep the pre-generated value.
    """
    remap = build_remap(keep_indices)
    keepset = set(remap)
    K = len(keep_indices)
    data = [_subset_trace(t, remap, keepset, renumber) for t in fig.get('data', [])]
    layout = dict(fig.get('layout', {}))
    # Width and x-range are sized to the column count; the full figure carries the
    # 4749-column values, so resize them to the trimmed K (matches make_plot_for_protein).
    if 'width' in layout:
        layout['width'] = max(K * 10, 1200)
    if isinstance(layout.get('xaxis'), dict):
        xaxis = dict(layout['xaxis'])
        if 'range' in xaxis:
            xaxis['range'] = [0.5, K + 1]
        if is_last is not None:
            xaxis['visible'] = bool(is_last)
        layout['xaxis'] = xaxis
    anns = layout.get('annotations')
    if anns:
        new_anns = []
        for a in anns:
            col = round(a['x']) if isinstance(a.get('x'), (int, float)) else None
            if col in keepset:
                a = dict(a, x=remap[col])
                if renumber and isinstance(a.get('text'), str):
                    a['text'] = _rewrite_align_pos(a['text'], col, remap)
                new_anns.append(a)
        layout['annotations'] = new_anns
    return {'data': data, 'layout': layout}
