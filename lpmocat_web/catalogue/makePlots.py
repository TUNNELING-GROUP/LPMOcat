import plotly.graph_objs
from plotly.subplots import make_subplots
from scipy.signal import residue
import itertools
from .models import Protein, Substrate, Clade, Organism, Flexibility
import numpy as np
import textwrap
import plotly.graph_objects as go
from .sec_struct_plot import sec_struct_plot
from .surface import SURFACE_RSA, read_sasa_file
from .second_sphere import motif_residue_text, second_sphere
from .brace_motifs import brace_motifs, motif_residue_text as brace_motif_residue_text
import json
import colorsys
import math
import re
from functools import lru_cache
import csv
import os
from pathlib import Path

# Catalytic domain + CBM ranges per structure, in the structure's own (auth) residue
# numbering. Produced by the `map_domain_ranges` management command.
_DOMAIN_RANGES_PATH = Path(__file__).resolve().parents[2] / 'data' / 'django_outputs' / 'domain_ranges_auth.json'
try:
    with open(_DOMAIN_RANGES_PATH) as _f:
        DOMAIN_RANGES = json.load(_f)
except (FileNotFoundError, ValueError):
    DOMAIN_RANGES = {}

# Refined catalytic-domain end residues (pipeline step 6): THE catalytic-domain
# assignment used across the app (3D view, plots, domain bar, dendrogram, lengths, masses,
# charges). The annotation-derived catalytic range in DOMAIN_RANGES can overshoot (a fallback
# to the last residue, a loose annotation) or stop short, so it is overridden by this end.
_CAT_BOUNDARIES_PATH = Path(__file__).resolve().parents[2] / 'data' / '06_catalytic_domain' / 'catalytic_domain_boundaries.json'
try:
    with open(_CAT_BOUNDARIES_PATH) as _f:
        CATALYTIC_BOUNDARIES = json.load(_f)
except (FileNotFoundError, ValueError):
    CATALYTIC_BOUNDARIES = {}


def apply_catalytic_boundary(ranges, cat_end):
    """Domain ranges with the catalytic segment ending exactly at `cat_end` (the refined
    catalytic-domain assignment), shortening or extending the annotated end. Other domains
    (CBM, …) keep their annotated range, clipped to start after the catalytic domain; one left
    empty is dropped. Unchanged when cat_end is None."""
    if cat_end is None:
        return list(ranges)
    result = []
    for dom in ranges:
        if dom.get('type') == 'catalytic':
            dom = {**dom, 'end_residue': cat_end}
        elif dom.get('start_residue') is not None and dom['start_residue'] <= cat_end:
            if dom.get('end_residue') is not None and dom['end_residue'] <= cat_end:
                continue
            dom = {**dom, 'start_residue': cat_end + 1}
        result.append(dom)
    return result


def _capped_domain_ranges(struct_id):
    """Domain ranges for a structure with the catalytic domain as assigned in
    CATALYTIC_BOUNDARIES (see apply_catalytic_boundary)."""
    key = (struct_id or '').strip()
    return apply_catalytic_boundary(DOMAIN_RANGES.get(key, []),
                                    CATALYTIC_BOUNDARIES.get(key, {}).get('end_residue'))

# Domain bar: a slim band drawn above the sequence row. Solid (opaque) colours — not rgba:
# Plotly splits an rgba fill into fill + fill-opacity, and the frontend's per-hover
# style.fill alpha would then multiply with that lingering fill-opacity, giving the wrong
# shades (dark-by-default, lighter-after-unhover). Solid colours avoid that — light by
# default, saturated on hover (the hover colour is also the tooltip bg, so they match).
DOMAIN_BAR_COLORS = {
    'catalytic': 'rgb(187, 214, 200)',   # light site green
    'CBM': 'rgb(239, 207, 172)',         # light amber
}
DOMAIN_LABEL_COLORS = {
    'catalytic': '#3f7359',
    'CBM': '#a2661f',
}
# Saturated hues — applied by the frontend while the bar is hovered (and used as the
# tooltip background), so the darkened band and its tooltip are exactly the same colour.
DOMAIN_BAR_HOVER_COLORS = {
    'catalytic': 'rgb(86, 153, 118)',
    'CBM': 'rgb(217, 141, 58)',
}
# Short on-plot labels — the raw InterPro names are far too long to draw
# (e.g. "Cellulose/chitin-binding protein, N-terminal"). Keyed by InterPro accession.
DOMAIN_SHORT_LABELS = {
    'IPR000254': 'CBM1',    'IPR035971': 'CBM1',
    'IPR001919': 'CBM2',    'IPR012291': 'CBM2',
    'IPR008965': 'CBM2/3',
    'IPR001956': 'CBM3',    'IPR036966': 'CBM3',
    'IPR003610': 'CBM5/12', 'IPR036573': 'CBM5/12',
    'IPR009031': 'CBM10',   'IPR002883': 'CBM10',  'IPR036601': 'CBM10',
    'IPR002044': 'CBM20',   'IPR034836': 'CBM20',
}
# The sequence heatmap is a single row at y=0, spanning -0.5..0.5. The domain bar is a
# slim band sitting entirely above it, with a little clearance.
DOMAIN_BAR_Y0 = 0.72
DOMAIN_BAR_Y1 = 1.00
DOMAIN_BAR_YMAX = 1.08   # top of the y-axis range, so the bar isn't clipped

# Height of the residue/data row — the coloured band the user actually looks at.
# Every non-sequence plot is exactly this tall. The Sequence plot is taller by the
# domain bar's share, but its residue row still renders at DATA_ROW_HEIGHT so all the
# data rows line up visually.
#
# Sequence y-range is [-0.5, DOMAIN_BAR_YMAX]; the residue row spans -0.5..0.5 (1.0 of
# the 1.58 total), so SEQUENCE_PLOT_HEIGHT = DATA_ROW_HEIGHT * (1.58 / 1.0).
# These MUST match `.plotWrapper` / `.yaxis-label` in compare_protein.css — the CSS
# forces the plot div height, so a mismatch squashes the SVG and drifts the labels.
DATA_ROW_HEIGHT = 36
SEQUENCE_PLOT_HEIGHT = round(DATA_ROW_HEIGHT * (0.5 + 1.08) / 1.0)   # -> 57

# Hover tooltip font size, applied at the layout level to every plot.
HOVER_FONT_SIZE = 10

# The Secondary Structure and Loops plot: the secondary structure row at y=0 (-0.5..0.5, like
# the sequence row; glyphs scaled up to the fixed range), the core-strand bar above it, the
# loop bar on top and a dot on each disulfide-bonded Cys below it. Its height must match
# compare_protein.css.
SS_LOOPS_PLOT_TYPE = 'Secondary Structure and Loops'
SS_GLYPH_SCALE = 4.5
CORE_BAR_Y0, CORE_BAR_Y1 = 0.62, 0.90
LOOP_BAR_Y0, LOOP_BAR_Y1 = 0.98, 1.26
SS_BARS_YMAX = 1.32
BINDING_TICK_Y = -0.62          # binding-surface ticks, right below the secondary structure row
DISULFIDE_DOT_Y = -1.00
SS_LANE_YMIN = -1.20
SS_LOOPS_PLOT_HEIGHT = round(DATA_ROW_HEIGHT * (SS_BARS_YMAX - SS_LANE_YMIN) / 1.0)   # -> 91
CORE_STRAND_COLOR = (251, 192, 45)   # the strand colour of the secondary structure row
# Disulfide marks and sticks: dark gold, sulfur-like but far darker than the strand yellow
# (3.3:1 on white)
DISULFIDE_COLOR = (184, 134, 11)
# Binding-surface ticks (pipeline step 7c): deep teal, the one hue the plot does not already
# use -- helix red, strand and core-strand yellow, coil grey, the loop colours (orange, violet,
# blue, green) and the disulfide gold (4.8:1 on white)
BINDING_SURFACE_COLOR = (0, 121, 107)

# Sequence plot: the second-sphere motif residues (pipeline 8c; the aromatic keeps the axial
# colour) and the shielding residue over the copper's solvent side (8d). Teal and brown:
# OKLab dE >= 9.7 from each other and every other sequence colour under protanopia and
# deuteranopia (>= 17.1 in normal vision); black letters on them 7.0:1 and 4.5:1.
SEQ_MOTIF = 'Second-sphere motif'
SEQ_DISTAL = 'Shielding residue'
SEQ_MOTIF_COLOR = (38, 166, 154)
SEQ_DISTAL_COLOR = (141, 110, 99)
# Sequence plot: the Gly, Pro and Arg of the motifs after the brace histidines (pipeline 8e,
# H-G-X4-P-X2-R after the N-terminal His, H-X-G-P after the internal one). Blue: OKLab dE >= 10.8
# from every other sequence colour under protanopia and deuteranopia (>= 18.4 in normal
# vision); black letters on it 6.2:1.
SEQ_BRACE_MOTIF = 'Brace histidine motif'
SEQ_BRACE_MOTIF_COLOR = (68, 136, 252)


def core_strand_label(name):
    """Display name of a core strand: 'b3' -> 'β3'."""
    return 'β' + name[1:] if name.startswith('b') else name

# Secondary-structure residue colours for Mol* — match the SS plot (helix red, sheet
# yellow, undefined/coil grey). Keyed by the codes stored in StructList.
SS_MOLSTAR_COLORS = {
    'H': (198, 40, 40),    # helix  (#C62828)
    'E': (251, 192, 45),   # sheet  (#FBC02D)
    'U': (97, 97, 97),     # undefined / coil (#616161)
}


# Relative SASA is thresholded at 25%: below = buried, above = exposed. Muted (low
# intensity) blue/red — used for BOTH the RSA plot and the Mol* coloring so they match.
RSA_THRESHOLD = SURFACE_RSA
RSA_BURIED_RGB = (95, 140, 200)     # blue
RSA_EXPOSED_RGB = (208, 100, 88)    # red


def _sasa_to_rgb(rsa):
    """Relative SASA (%) → colour matching the RSA plot: muted blue when buried
    (below the threshold), muted red when exposed."""
    if rsa is None:
        return (255, 255, 255)
    return RSA_EXPOSED_RGB if float(rsa) >= RSA_THRESHOLD else RSA_BURIED_RGB


# Hole hopping: every protein in a comparison gets its own colour, in comparison order; within a
# protein, surface-exposed donor residues use a clearly darker shade than bridging residues.
# The palette is generated in OKLCH (tune the look with the HH_* parameters): soft chroma, hues
# with magenta / pink / rose skipped (which also keeps clear of the Mol* hover pink) and orange /
# brown / tan / brownish olive skipped (reserved for the orange copper dot of the path overlay),
# and two lightness tiers with the same bridge->donor step, giving 64 (bridge, donor) pairs. They
# are ordered greedily, each next protein taking the pair most different from all earlier ones,
# so the first proteins are the most distinct and every protein still gets its own colour (with
# many proteins, neighbouring colours become similar; past 64 the list repeats).
HH_HUES = [24] + list(range(114, 295, 6))    # OKLCH hue, degrees: salmon, then yellow-green .. violet
HH_CHROMA = (0.085, 0.10)                    # (bridge, donor)
HH_TIERS = ((0.81, 0.56), (0.75, 0.50))      # (bridge L, donor L) per lightness tier
HH_FIRST_HUE = 252                           # the first protein is blue
HOLE_HOPPING_NONE_HOVER = 'rgba(176,190,197,0.7)'


def _hh_srgb_channel(c):
    return 12.92 * c if c <= 0.0031308 else 1.055 * c ** (1 / 2.4) - 0.055


def _hh_linear_channel(c):
    return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4


def _hh_oklch_to_rgb(L, C, hue):
    """OKLCH -> 8-bit sRGB, lowering chroma until the colour fits the sRGB gamut."""
    while True:
        a, b = C * math.cos(math.radians(hue)), C * math.sin(math.radians(hue))
        l_ = (L + 0.3963377774 * a + 0.2158037573 * b) ** 3
        m_ = (L - 0.1055613458 * a - 0.0638541728 * b) ** 3
        s_ = (L - 0.0894841775 * a - 1.2914855480 * b) ** 3
        lin = (4.0767416621 * l_ - 3.3077115913 * m_ + 0.2309699292 * s_,
               -1.2684380046 * l_ + 2.6097574011 * m_ - 0.3413193965 * s_,
               -0.0041960863 * l_ - 0.7034186147 * m_ + 1.7076147010 * s_)
        if all(-1e-4 <= v <= 1 + 1e-4 for v in lin) or C <= 0:
            return tuple(round(max(0.0, min(1.0, _hh_srgb_channel(max(0.0, v)))) * 255) for v in lin)
        C -= 0.002


def _hh_rgb_to_oklab(rgb):
    r, g, b = (_hh_linear_channel(v / 255) for v in rgb)
    l_ = (0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b) ** (1 / 3)
    m_ = (0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b) ** (1 / 3)
    s_ = (0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b) ** (1 / 3)
    return (0.2104542553 * l_ + 0.7936177850 * m_ - 0.0040720468 * s_,
            1.9779984951 * l_ - 2.4285922050 * m_ + 0.4505937099 * s_,
            0.0259040371 * l_ + 0.7827717662 * m_ - 0.8086757660 * s_)


@lru_cache(maxsize=1)
def hole_hopping_palette():
    """All (bridge RGB, donor RGB) pairs in protein order: each next pair is the candidate most
    different from every earlier one (the smallest OKLab distance between any of their colours)."""
    candidates = []
    for tier, (bridge_L, donor_L) in enumerate(HH_TIERS):
        for hue in HH_HUES:
            pair = (_hh_oklch_to_rgb(bridge_L, HH_CHROMA[0], hue), _hh_oklch_to_rgb(donor_L, HH_CHROMA[1], hue))
            candidates.append((tier, hue, pair, tuple(_hh_rgb_to_oklab(c) for c in pair)))

    def distance(p, q):
        return min(math.dist(p[3][i], q[3][j]) for i in (0, 1) for j in (0, 1))

    first = min(candidates, key=lambda c: (c[0], abs(c[1] - HH_FIRST_HUE)))
    order, rest = [first], [c for c in candidates if c is not first]
    nearest = {id(c): distance(c, first) for c in rest}
    while rest:
        nxt = max(rest, key=lambda c: nearest[id(c)])
        rest.remove(nxt)
        order.append(nxt)
        for c in rest:
            nearest[id(c)] = min(nearest[id(c)], distance(c, nxt))
    return tuple(c[2] for c in order)


def hole_hopping_colors(position):
    """{'bridge': RGB, 'donor': RGB} for the protein at `position` (0-based) in a comparison."""
    palette = hole_hopping_palette()
    bridge, donor = palette[position % len(palette)]
    return {'bridge': bridge, 'donor': donor}


def _hh_hover_rgba(rgb):
    return f'rgba({rgb[0]},{rgb[1]},{rgb[2]},0.9)'


def recolor_hole_hopping_figure(fig, position):
    """Re-apply the colour pair for comparison `position` to a Hole Hopping figure dict.
    Pre-generated figures are built one protein at a time (position 0), so the selection's
    order is applied when they are assembled. Roles are read from the heatmap z (1 bridge,
    2 donor) and the hover text."""
    colors = hole_hopping_colors(position)
    bridge = 'rgb({},{},{})'.format(*colors['bridge'])
    donor = 'rgb({},{},{})'.format(*colors['donor'])
    for trace in fig.get('data', []):
        if trace.get('type') == 'heatmap' and trace.get('name') == 'hole_hopping':
            trace['colorscale'] = [[0.0, bridge], [0.5, bridge], [1.0, donor]]
            trace['opacity'] = 1   # matches makeHoleHoppingPlots (older pre-generated files used 0.7)
        elif trace.get('type') == 'scatter' and isinstance((trace.get('hoverlabel') or {}).get('bgcolor'), list):
            trace['hoverlabel']['bgcolor'] = [
                _hh_hover_rgba(colors['donor']) if 'Donor' in str(text)
                else _hh_hover_rgba(colors['bridge']) if 'Bridge' in str(text)
                else HOLE_HOPPING_NONE_HOVER
                for text in (trace.get('customdata') or [])]
    return fig


def parse_pdb_residue_numbers(pdb_text):
    """Return the author residue numbers (auth_seq_id) of a PDB structure, in
    chain/file order, one entry per residue. These are the original numbers from
    the source file — they do not necessarily start at 1 and may contain gaps."""
    numbers = []
    seen = set()
    for line in (pdb_text or '').split('\n'):
        if not line.startswith('ATOM'):
            continue
        # auth_seq_id occupies cols 23-26, insertion code col 27 (1-indexed)
        key = (line[21], line[22:27])
        if key in seen:
            continue
        seen.add(key)
        try:
            numbers.append(int(line[22:26]))
        except ValueError:
            continue
    return numbers


def catalytic_masked_sequence(aligned, sequence, pdb_structure, struct_id):
    """`aligned` with every residue *beyond* the structure's catalytic-domain end_residue
    (CATALYTIC_BOUNDARIES, auth numbering) replaced by a gap ('-'), so CBMs, linkers and
    disordered tails drop out. Unchanged when no boundary is known. Residues are numbered
    from the PDB, falling back to 1..len(sequence) when the PDB doesn't match the sequence.
    Auth numbering is monotonic, so only trailing residues are masked."""
    cat_end = CATALYTIC_BOUNDARIES.get((struct_id or '').strip(), {}).get('end_residue')
    if cat_end is None:
        return aligned
    res_nums = parse_pdb_residue_numbers(pdb_structure)
    seq_len = len(sequence or '')
    if len(res_nums) != seq_len:
        res_nums = list(range(1, seq_len + 1))

    out, seq_pos = [], 0
    for ch in aligned:
        if ch == '-':
            out.append('-')
            continue
        auth = res_nums[seq_pos] if seq_pos < len(res_nums) else seq_pos + 1
        seq_pos += 1
        out.append(ch if auth <= cat_end else '-')
    return ''.join(out)


# Flexibility methods the viewer can show, in display order, with per-method hover labels and
# units. Values live in the Flexibility model (one row per method per protein); a protein may
# have only some methods (e.g. B-factor for crystals vs pLDDT for AlphaFold models). Each
# method has its own value range/units, so colour gradients and y-axes are scaled per method.
FLEX_METHOD_ORDER = ['cabsflex', 'gnm', 'md_rmsf', 'bfactor', 'plddt']
# `scheme='alphafold'`: colour by the official AlphaFold pLDDT confidence bands instead of the
# red flexibility gradient (pLDDT is confidence, not a fluctuation magnitude).
FLEX_METHOD_META = {
    'cabsflex': {'label': 'CABSflex RMSF', 'hover': 'RMSF',      'unit': ' Å'},
    'gnm':      {'label': 'GNM',           'hover': 'GNM',       'unit': ''},
    'md_rmsf':  {'label': 'MD RMSF',       'hover': 'RMSF',      'unit': ' Å'},
    'bfactor':  {'label': 'B-factor',      'hover': 'B-factor',  'unit': ' Å²'},
    'plddt':    {'label': 'pLDDT',         'hover': 'pLDDT',     'unit': '', 'scheme': 'alphafold'},
}


def _plddt_rgb(v):
    """Official AlphaFold pLDDT confidence-band colour for a value (0-100)."""
    if v >= 90:   return (0, 83, 214)      # Very high   (dark blue)
    if v >= 70:   return (101, 203, 243)   # Confident   (light blue)
    if v >= 50:   return (255, 219, 19)    # Low         (yellow)
    return (255, 125, 69)                  # Very low    (orange)


class makePlots:

    # Color per region family — sub-regions share parent color
    MOLSTAR_COLOR_BY_GROUP = {
        'L2a': (255, 69, 0),
        'L2': (255, 69, 0),
        'L2b': (255, 69, 0),
        'L3': (255, 140, 0),
        'LS': (138, 43, 226),
        'L8': (30, 144, 255),
        'LC': (34, 139, 34),
        'LCb': (34, 139, 34),
    }

    # This class is used to create plots for the LPMO catalogue.
    def __init__(self, protObjList:list, msa:list=None, no_trim:bool=False, keep_columns=None):
        # no_trim=True keeps every alignment column (the full MSA width) instead of
        # dropping columns that are all-gap across the selected set. Used to
        # pre-generate per-protein plots against the fixed full-alignment coordinate
        # system, which can then be column-subset per selection.
        # keep_columns=<0-indexed column list> trims to exactly those alignment columns
        # (used to restrict the plots to the catalytic domain); overrides no_trim.
        self._no_trim = no_trim
        self._keep_columns = keep_columns

        # Load all hydrophobicity scales from CSV
        self.allHydrophobicityScales = {}
        self.hydrophobicityScaleNames = []
        csv_path = os.path.join(Path(__file__).parent, 'data', 'hydrophobicity_scales.csv')
        try:
            with open(csv_path, 'r') as f:
                reader = csv.DictReader(f)
                scale_columns = [col for col in reader.fieldnames if col not in ('aa', 'aa full', 'aa classes')]
                self.hydrophobicityScaleNames = scale_columns
                for scale_name in scale_columns:
                    self.allHydrophobicityScales[scale_name] = {'-': None}
                for row in reader:
                    aa = row['aa']
                    for scale_name in scale_columns:
                        self.allHydrophobicityScales[scale_name][aa] = float(row[scale_name])
        except Exception as e:
            print(f"Warning: Could not load hydrophobicity scales CSV: {e}")
            self.allHydrophobicityScales = {'Kyte-Doolittle': {
                "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8, "G": -0.4, "H": -3.2,
                "I": 4.5, "K": -3.9, "L": 3.8, "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5,
                "R": -4.5, "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3, "-": None
            }}
            self.hydrophobicityScaleNames = ['Kyte-Doolittle']
        self.defaultHydrophobScale = self.hydrophobicityScaleNames[0]
        self.aaHydrophobicity = self.allHydrophobicityScales[self.defaultHydrophobScale]
        three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        '---': '-',
        }

        one_to_three = {
            'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
            'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
            'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
            'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
            '-': '---'
        }



        iterator = 0
        self.oneLetterSequence = []
        self.threeLetterSequence = []
        self.flexibilityList = []
        # Flexibility for every method the selected proteins have, so the user can switch which
        # measure is shown (self.flexibilityList stays the default method for legacy callers).
        _flex_present = set()
        for _p in protObjList:
            _flex_present.update(_p.flexibility.values_list('method', flat=True))
        self.flexibilityMethods = [m for m in FLEX_METHOD_ORDER if m in _flex_present]
        self.defaultFlexMethod = ('cabsflex' if 'cabsflex' in self.flexibilityMethods
                                  else (self.flexibilityMethods[0] if self.flexibilityMethods else 'cabsflex'))
        self.allFlexibilityLists = {m: [] for m in self.flexibilityMethods}
        # flexMin/MaxByMethod follow the built sequence (per-view: catalytic when the sequence is
        # catalytic-masked, full otherwise) and drive the PLOT. flexMin/MaxByMethodCat is always
        # the catalytic-domain range (residues <= end_residue), used to normalise the default 3D
        # colouring so the catalytic view matches; full-range 3D colouring uses flexMin/MaxByMethod.
        self.flexMinByMethod = {m: None for m in self.flexibilityMethods}
        self.flexMaxByMethod = {m: None for m in self.flexibilityMethods}
        self.flexMinByMethodCat = {m: None for m in self.flexibilityMethods}
        self.flexMaxByMethodCat = {m: None for m in self.flexibilityMethods}
        self.relSASAList = []
        self.hydrophobicityList = []
        self.allHydrophobicityLists = {name: [] for name in self.hydrophobicityScaleNames}
        self.electrostaticsList = []
        self.atomChargesList = []  # Store atom charges for molstar coloring
        self.StructList = []
        self.bindingSurfaceYN = []
        self.loopRegionsYN = []
        self.histidineBraceTuple = []
        self.bindingSurfaceDict = []
        self.loopRegionsDict = []
        self.betaCoreDict = []      # per protein {'β1': [[start, end], ...]} (pieces, sequence positions)
        self.disulfideList = []     # per protein [[residue, residue, S-S distance], ...] (residue numbers)
        self.residueNumList = [1]
        # Original PDB residue numbers per protein (one per non-gap residue, in
        # sequence order). Falls back to 1..N when no structure is stored.
        self.residueNumbers = []
        self.protId = [id.protein_code for id in protObjList]
        self.protName = [id.protein_name for id in protObjList]
        self.sequenceGaps = []
        self.braceInteractionResidues = []
        self.metalResidueNumbers = []
        # Catalytic domain / CBM ranges (auth numbering) per protein, for the domain bar
        self.domainRanges = [
            _capped_domain_ranges(p.selected_structure) for p in protObjList
        ]
        # Catalytic-domain end residue (auth numbering) per protein, None when unknown
        self.catalyticEnds = [
            CATALYTIC_BOUNDARIES.get((p.selected_structure or '').strip(), {}).get('end_residue')
            for p in protObjList
        ]
        # Total proteins in this comparison (kept when subsetting to single-protein plots),
        # so hovers can show the alignment column only when more than one is displayed.
        self.totalProteinCount = len(protObjList)
        self.holeHoppingList = []
        self.holeHoppingPaths = []
        self.maxFlexVal = None
        self.minFlexVal = None
        self.maxElectrostaticsVal = None
        self.minElectrostaticsVal = None
        # self.secondSphereResidues = [eval(prot.second_sphere_aromatic) for prot in protObjList]
        def safe_eval_second_sphere(value):
            if value is None or value == '' or str(value).lower() == 'nan':
                return None
            try:
                result = eval(str(value))
                if isinstance(result, (tuple, list)):
                    return result
                return (result,)
            except (NameError, SyntaxError):
                return None
        self.secondSphereResidues = [safe_eval_second_sphere(prot.second_sphere_aromatic) for prot in protObjList]
        # Second-sphere motif and shielding residue per protein (second_sphere())
        self.secondSphereMotifs = [second_sphere(prot) for prot in protObjList]
        # Motifs after the N-terminal and the internal His of the brace (brace_motifs())
        self.braceMotifs = [brace_motifs(prot) for prot in protObjList]

        
        for protein in protObjList:

            # Original residue numbering from the structure file. The stored
            # sequence corresponds 1:1 to the ATOM residues, so a positional
            # lookup maps sequence index -> original residue number.
            resNums = parse_pdb_residue_numbers(protein.pdb_structure)
            if len(resNums) != len(protein.sequence or ''):
                resNums = list(range(1, len(protein.sequence or '') + 1))
            self.residueNumbers.append(resNums)

            if msa is not None and protein.protein_code in msa['name']:
                idx = msa['name'].index(protein.protein_code)
                aligned_seq = msa['seq'][idx]
                self.oneLetterSequence.append(list(aligned_seq))
                self.threeLetterSequence.append([one_to_three.get(aa, '---') for aa in aligned_seq])
            else:
                # protein.sequence is 1-letter codes (e.g. "HGQVQN...")
                one_letter = list(protein.sequence)
                self.oneLetterSequence.append(one_letter)
                self.threeLetterSequence.append(
                    [one_to_three.get(aa, '---') for aa in one_letter]
                )

            # Build a per-residue list (aligned to the sequence, None at gaps) for EVERY
            # flexibility method this protein has; track each method's global min/max for its
            # colour gradient. A method the protein lacks -> all None (empty plot for it).
            cat_end = CATALYTIC_BOUNDARIES.get((protein.selected_structure or '').strip(), {}).get('end_residue')
            flex_by_method = {fx.method: fx for fx in protein.flexibility.all()}
            for method in self.flexibilityMethods:
                fx = flex_by_method.get(method)
                try:
                    method_values = json.loads(fx.rmsf)["values"] if fx else None
                except (ValueError, KeyError, TypeError):
                    method_values = None
                # Catalytic-domain range: min/max over residues with auth number <= end_residue
                # (view-independent). method_values is in residue order, so index k -> resNums[k].
                if method_values:
                    for k, x in enumerate(method_values):
                        try:
                            fv = float(x)
                        except (ValueError, TypeError):
                            continue
                        auth = resNums[k] if k < len(resNums) else k + 1
                        if cat_end is None or auth <= cat_end:
                            if self.flexMinByMethodCat[method] is None or fv < self.flexMinByMethodCat[method]:
                                self.flexMinByMethodCat[method] = fv
                            if self.flexMaxByMethodCat[method] is None or fv > self.flexMaxByMethodCat[method]:
                                self.flexMaxByMethodCat[method] = fv
                # Per-residue list aligned to the (possibly masked) sequence; None at gaps/masked.
                # flexMin/MaxByMethod tracks the values actually in the sequence -> per-view range.
                seq_vals = []
                vi = 0
                for residue in self.oneLetterSequence[iterator]:
                    if residue == '-':
                        seq_vals.append(None)
                        continue
                    val = None
                    if method_values is not None and vi < len(method_values):
                        try:
                            val = float(method_values[vi])
                        except (ValueError, TypeError):
                            val = None
                    vi += 1
                    seq_vals.append(val)
                    if val is not None:
                        if self.flexMinByMethod[method] is None or val < self.flexMinByMethod[method]:
                            self.flexMinByMethod[method] = val
                        if self.flexMaxByMethod[method] is None or val > self.flexMaxByMethod[method]:
                            self.flexMaxByMethod[method] = val
                self.allFlexibilityLists[method].append(seq_vals)

            # Legacy single-method attributes = the default method (keeps export / 3D query /
            # subsetting / the plot's default working unchanged).
            if self.defaultFlexMethod in self.allFlexibilityLists:
                self.flexibilityList.append(self.allFlexibilityLists[self.defaultFlexMethod][iterator])
                self.minFlexVal = self.flexMinByMethod.get(self.defaultFlexMethod)
                self.maxFlexVal = self.flexMaxByMethod.get(self.defaultFlexMethod)
            else:
                self.flexibilityList.append([None] * len(self.oneLetterSequence[iterator]))

            # relative SASA (%) per residue of the structure sequence (NaN where unknown)
            sasaRows = read_sasa_file(protein.surface_accesability.sasa_file)
            self.relSASAList.append([])
            residue = 0
            for res in self.oneLetterSequence[iterator]:
                if res != '-':
                    self.relSASAList[iterator].append(
                        sasaRows[residue][3] if residue < len(sasaRows) else math.nan)
                    residue += 1


            self.hydrophobicityList.append(list((map(lambda x: self.aaHydrophobicity[x], self.oneLetterSequence[iterator]))))
            for scale_name in self.hydrophobicityScaleNames:
                scale_dict = self.allHydrophobicityScales[scale_name]
                self.allHydrophobicityLists[scale_name].append(
                    list(map(lambda x: scale_dict.get(x), self.oneLetterSequence[iterator]))
                )

            # Extract charge data from PQR file
            self.electrostaticsList.append([])
            self.atomChargesList.append([])  # Store atom charges: list of dicts with residue_number, atom_name, charge
            total_charge = None
            residue_charges = {}  # residue_number -> total charge
            
            try:
                # Get PQR file content from electrostatics model
                pqr_content = protein.electrostatics.pqr_file if hasattr(protein, 'electrostatics') and protein.electrostatics else None
                
                if pqr_content:
                    pqr_lines = pqr_content.splitlines()
                    
                    # Extract total charge from REMARK line
                    import re
                    for line in pqr_lines:
                        if 'REMARK' in line and 'Total charge' in line:
                            try:
                                # Extract charge value using regex: "Total charge on this protein: -8.0000 e"
                                # Look for number (with optional sign and decimal) followed by 'e'
                                match = re.search(r'([-+]?\d+\.?\d*)\s+e', line)
                                if match:
                                    total_charge = float(match.group(1))
                                else:
                                    # Fallback: try splitting by colon
                                    if ':' in line:
                                        charge_part = line.split(':')[1].strip()
                                        charge_str = charge_part.split('e')[0].strip()
                                        total_charge = float(charge_str)
                            except (ValueError, IndexError, AttributeError):
                                pass
                            break
                    
                    # Extract per-residue charges and atom charges from ATOM lines
                    current_residue = None
                    current_residue_charge = 0.0
                    
                    for line in pqr_lines:
                        if line.startswith('ATOM'):
                            parts = line.split()
                            if len(parts) >= 10:
                                try:
                                    # The chain ID column is optional in PQR, so index from the
                                    # end: ... resSeq x y z charge radius
                                    residue_num = int(parts[-6])
                                    charge = float(parts[-2])
                                    atom_name = parts[2] if len(parts) > 2 else ''  # Atom name
                                    
                                    # Store atom charge for molstar coloring
                                    self.atomChargesList[iterator].append({
                                        'residue_number': residue_num,
                                        'atom_name': atom_name,
                                        'charge': charge
                                    })
                                    
                                    if current_residue is None:
                                        current_residue = residue_num
                                        current_residue_charge = 0.0
                                    
                                    if residue_num != current_residue:
                                        # Save previous residue charge
                                        residue_charges[current_residue] = current_residue_charge
                                        current_residue = residue_num
                                        current_residue_charge = 0.0
                                    
                                    current_residue_charge += charge
                                except (ValueError, IndexError):
                                    continue
                    
                    # Don't forget the last residue
                    if current_residue is not None:
                        residue_charges[current_residue] = current_residue_charge
                
                # Map charges to sequence positions (handling gaps). residue_charges is
                # keyed by the original residue numbers from the PQR/PDB, so look up by
                # those rather than by a 1-based position.
                seqPos = 0
                protResNums = self.residueNumbers[iterator]
                for res in self.oneLetterSequence[iterator]:
                    if res != '-':
                        resNum = protResNums[seqPos] if seqPos < len(protResNums) else seqPos + 1
                        seqPos += 1
                        # Get charge for this residue number, default to 0 if not found
                        charge_val = residue_charges.get(resNum, 0.0)
                        self.electrostaticsList[iterator].append(charge_val)
                        
                        # Track min/max for normalization
                        if self.minElectrostaticsVal is None or charge_val < self.minElectrostaticsVal:
                            self.minElectrostaticsVal = charge_val
                        if self.maxElectrostaticsVal is None or charge_val > self.maxElectrostaticsVal:
                            self.maxElectrostaticsVal = charge_val
                    else:
                        self.electrostaticsList[iterator].append(None)
                
                # Store total charge for this protein (for filtering)
                if not hasattr(self, 'totalChargeList'):
                    self.totalChargeList = []
                self.totalChargeList.append(total_charge)
                
            except Exception as e:
                print(f"Warning: Could not extract charge data for {protein.protein_code}: {e}")
                # Fill with None values
                for res in self.oneLetterSequence[iterator]:
                    self.electrostaticsList[iterator].append(None)
                if not hasattr(self, 'totalChargeList'):
                    self.totalChargeList = []
                self.totalChargeList.append(None)
                # Ensure atomChargesList has an entry even if empty
                if len(self.atomChargesList) <= iterator:
                    self.atomChargesList.append([])

            try:
                binding_dict = eval(protein.binding_surface) if protein.binding_surface else {}
            except Exception:
                binding_dict = {}
            self.bindingSurfaceDict.append(binding_dict if isinstance(binding_dict, dict) else {})

            loop_regions_raw = getattr(protein, 'loop_regions', None)
            loop_regions_dict = {}
            if loop_regions_raw and str(loop_regions_raw).strip().lower() not in ('nan', 'none', ''):
                try:
                    parsed_loop_regions = eval(loop_regions_raw)
                    if isinstance(parsed_loop_regions, dict):
                        loop_regions_dict = parsed_loop_regions
                except Exception as e:
                    print(f"WARNING: failed to evaluate loop regions for {protein.protein_code}: {e}")
                    loop_regions_dict = {}
            self.loopRegionsDict.append(loop_regions_dict)

            try:
                core = json.loads(protein.beta_core) if getattr(protein, 'beta_core', None) else {}
            except (TypeError, ValueError):
                core = {}
            self.betaCoreDict.append({core_strand_label(name): [piece.split('-') for piece in pieces]
                                      for name, pieces in core.items()})

            try:
                bonds = json.loads(protein.disulfides) if getattr(protein, 'disulfides', None) else []
            except (TypeError, ValueError):
                bonds = []
            self.disulfideList.append(bonds)

            # Secondary structure is stored on the cabsflex Flexibility row (DSSP simplified).
            _ss_obj = flex_by_method.get('cabsflex') or flex_by_method.get(self.defaultFlexMethod)
            ss_raw = _ss_obj.sec_struct if _ss_obj else None
            if ss_raw:
                ss_list = json.loads(ss_raw)["ss"]
            else:
                ss_list = ['-'] * len(resNums)

            self.StructList.append('')
            residueIndex = 0
            for res in self.oneLetterSequence[iterator]:
                if res != '-':
                    ss = ss_list[residueIndex] if residueIndex < len(ss_list) else 'U'
                    # DSSP coil is '-'; remap to 'U' so it's distinct from alignment gaps
                    if ss == '-':
                        ss = 'U'
                    self.StructList[iterator] += ss
                    residueIndex += 1

                else:
                    self.StructList[iterator] += '-'


            # {region: [[start, end], ...]}, sequence positions. Every range is kept: the plane
            # of pipeline step 7c touches a few residues per loop, so a region holds several
            # short runs, where the former templates gave one range per region.
            surface = {}
            for key, ranges in self.bindingSurfaceDict[iterator].items():
                pieces = ranges if isinstance(ranges, (list, tuple)) else [ranges]
                segments = []
                for piece in pieces:
                    if isinstance(piece, str):
                        parts = [value.strip() for value in piece.split('-')]
                    elif isinstance(piece, (list, tuple)) and len(piece) >= 2:
                        parts = [str(piece[0]).strip(), str(piece[1]).strip()]
                    else:
                        continue
                    try:
                        segments.append([int(parts[0]), int(parts[1])])
                    except (ValueError, IndexError):
                        continue
                if segments:
                    surface[key] = sorted(segments)
            self.bindingSurfaceDict[iterator] = surface

            for key, ranges in list(self.loopRegionsDict[iterator].items()):
                try:
                    if isinstance(ranges, (list, tuple)) and len(ranges) == 0:
                        del self.loopRegionsDict[iterator][key]
                        continue
                    first_segment = ranges[0] if isinstance(ranges, (list, tuple)) and ranges else ranges
                    if isinstance(first_segment, str):
                        start_end = [val.strip() for val in first_segment.split('-')]
                    elif isinstance(first_segment, (list, tuple)) and len(first_segment) >= 2:
                        start_end = [str(first_segment[0]).strip(), str(first_segment[1]).strip()]
                    else:
                        del self.loopRegionsDict[iterator][key]
                        continue
                    if len(start_end) == 2:
                        self.loopRegionsDict[iterator][key] = start_end
                    else:
                        del self.loopRegionsDict[iterator][key]
                except Exception:
                    del self.loopRegionsDict[iterator][key]
                    continue


            on_surface = {}
            for key, segments in self.bindingSurfaceDict[iterator].items():
                for start_val, end_val in segments:
                    for position in range(start_val, end_val + 1):
                        on_surface[position] = key
            self.bindingSurfaceYN.append(
                [f'Yes [{on_surface[position]}]' if position in on_surface else 'No'
                 for position in range(1, len(self.oneLetterSequence[iterator]) + 1)])

            self.loopRegionsYN.append([])
            for key in self.loopRegionsDict[iterator].keys():
                try:
                    start_val = int(self.loopRegionsDict[iterator][key][0])
                    end_val = int(self.loopRegionsDict[iterator][key][1])
                except (ValueError, TypeError, IndexError):
                    continue
                self.loopRegionsYN[iterator] += ['No'] * (start_val - len(self.loopRegionsYN[iterator]) - 1)
                self.loopRegionsYN[iterator] += [f'Yes [{key}]'] * (end_val - start_val + 1)
            self.loopRegionsYN[iterator] += ['No'] * (len(self.oneLetterSequence[iterator]) - len(self.loopRegionsYN[iterator]))

            # histidine_brace holds residue numbers in the structure's own numbering.
            # It may be missing, empty or 'nan' for uncharacterised proteins.
            try:
                brace = eval(protein.histidine_brace) if (protein.histidine_brace or '').strip() not in ('', 'nan') else []
            except (SyntaxError, NameError, TypeError):
                brace = []
            self.histidineBraceTuple.append(brace)

            # Find auth_seq_id of the active site metal ion from the stored PDB structure
            metal_res_num = None
            if protein.metal and protein.pdb_structure:
                metal_name = protein.metal.upper()
                for _line in protein.pdb_structure.splitlines():
                    if not _line.startswith('HETATM'):
                        continue
                    rname = _line[17:20].strip().upper().rstrip('0123456789')
                    if rname == metal_name:
                        try:
                            metal_res_num = int(_line[22:26])
                        except (ValueError, IndexError):
                            pass
                        break
            self.metalResidueNumbers.append(metal_res_num)

            interactions_raw = getattr(protein, 'histidine_brace_interactions', None)
            interaction_residues = []
            if interactions_raw and str(interactions_raw).lower() not in ['nan', 'none', 'null', '']:
                try:
                    interactions_json = json.loads(interactions_raw)
                    residues_set = set()
                    for hb in (interactions_json.get('hydrogen_bonds') or []):
                        res_id = hb.get('other_res_id')
                        if isinstance(res_id, int):
                            residues_set.add(res_id)
                    for pi in (interactions_json.get('pi_pi_stacking') or []):
                        res_id = pi.get('other_res_id')
                        if isinstance(res_id, int):
                            residues_set.add(res_id)
                    interaction_residues = sorted(residues_set)
                except (json.JSONDecodeError, TypeError, ValueError):
                    interaction_residues = []
            self.braceInteractionResidues.append(interaction_residues)

            # Parse hole hopping data
            hh_roles = {}  # residue_number -> 'donor' or 'bridge'
            prot_hh_paths = []  # list of paths with residue numbers for 3D visualization
            try:
                if hasattr(protein, 'hole_hopping') and protein.hole_hopping:
                    paths_data = json.loads(protein.hole_hopping.paths_json)
                    # A path ends at the active-site metal: e.g. CU301, CU999 (the pseudo-Cu the
                    # pipeline adds when a structure has no copper) or ZN/NI/MG when that metal
                    # stood in for it. Node names are upper case.
                    metal_node = re.compile(r'(CU|ZN|NI|MG|CO|MN|FE|MO)\d+$', re.IGNORECASE)
                    donor_residue_nums = set()
                    # First pass: identify all donor residue numbers
                    for path_entry in paths_data:
                        path = path_entry.get('path', [])
                        if path:
                            m = re.match(r'[A-Za-z]+(\d+)', path[0])
                            if m and not metal_node.match(path[0]):
                                donor_residue_nums.add(int(m.group(1)))
                    # Second pass: classify all residues and build path graph data
                    for path_entry in paths_data:
                        path_residue_nums = []
                        for residue_name in path_entry.get('path', []):
                            m = re.match(r'[A-Za-z]+(\d+)', residue_name)
                            if not m:
                                continue
                            res_num = int(m.group(1))
                            is_cu = bool(metal_node.match(residue_name))
                            path_residue_nums.append({'residue_number': res_num, 'name': residue_name, 'is_cu': is_cu})
                            if not is_cu:
                                if res_num in donor_residue_nums:
                                    hh_roles[res_num] = 'donor'
                                elif res_num not in hh_roles:
                                    hh_roles[res_num] = 'bridge'
                        if path_residue_nums:
                            prot_hh_paths.append({
                                'rank': path_entry.get('rank', 0),
                                'path': path_residue_nums,
                                'mean_residence_time': path_entry.get('mean_residence_time', 0)
                            })
            except Exception as e:
                print(f"Warning: Could not parse hole hopping for {protein.protein_code}: {e}")
            self.holeHoppingPaths.append(prot_hh_paths)

            self.holeHoppingList.append([])
            # hh_roles is keyed by the structure's own residue numbering
            protResNums = self.residueNumbers[iterator]
            seqPos = 0
            for res in self.oneLetterSequence[iterator]:
                if res != '-':
                    resNum = protResNums[seqPos] if seqPos < len(protResNums) else seqPos + 1
                    seqPos += 1
                    self.holeHoppingList[iterator].append(hh_roles.get(resNum))
                else:
                    self.holeHoppingList[iterator].append(None)

            if len(self.threeLetterSequence[iterator]) > len(self.residueNumList):
                self.residueNumList =  np.linspace(1, len(self.oneLetterSequence[iterator]), len(self.oneLetterSequence[iterator]),
                                dtype=int).astype(
                        dtype='str')

            iterator += 1

        self.secondSphereCombinedResidues = []
        for idx in range(len(self.secondSphereResidues)):
            combined = set()
            second = self.secondSphereResidues[idx]
            if isinstance(second, (list, tuple, set)):
                for item in second:
                    if isinstance(item, int):
                        combined.add(item)
            elif isinstance(second, int):
                combined.add(second)
            combined.update([res for res in self.braceInteractionResidues[idx] if isinstance(res, int)])
            if idx < len(self.histidineBraceTuple):
                combined.difference_update([res for res in self.histidineBraceTuple[idx] if isinstance(res, int)])
            self.secondSphereCombinedResidues.append(sorted(combined))
        
        # Extract categorization of second-sphere residues (position numbers by category)
        self.axialResiduePositions = []  # Position numbers of axial residues
        self.piStackingResiduePositions = []  # Position numbers of pi-stacking residues
        self.cationPiResiduePositions = []  # Position numbers of cation-pi residues
        self.hbondResiduePositions = []  # Position numbers of H-bond residues
        
        for idx, protein in enumerate(protObjList):
            axial_positions = set()
            pi_stacking_positions = set()
            cation_pi_positions = set()
            hbond_positions = set()
            
            try:
                # Get axial residue position from second_sphere_aromatic
                second = self.secondSphereResidues[idx]
                if isinstance(second, (list, tuple, set)):
                    for item in second:
                        if isinstance(item, int):
                            axial_positions.add(item)
                elif isinstance(second, int):
                    axial_positions.add(second)
                
                # Get categorized residues from histidine_brace_interactions
                interactions_raw = getattr(protein, 'histidine_brace_interactions', None)
                if interactions_raw and str(interactions_raw).lower() not in ['nan', 'none', 'null', '']:
                    interactions_json = json.loads(interactions_raw)
                    
                    # Debug for first protein (commented out)
                    # if idx == 0:
                    #     print(f"DEBUG EXTRACTION - Protein {idx}:")
                    #     print(f"  Found interactions JSON with keys: {interactions_json.keys()}")
                    #     print(f"  Pi-pi stacking entries: {interactions_json.get('pi_pi_stacking', [])}")
                    #     print(f"  H-bond entries: {interactions_json.get('hydrogen_bonds', [])}")
                    
                    # Extract pi-pi stacking positions (exclude axial)
                    for pi in (interactions_json.get('pi_pi_stacking') or []):
                        res_num = pi.get('other_res_id')  # Fixed: was 'other_res_num'
                        if res_num and isinstance(res_num, int) and res_num not in axial_positions:
                            pi_stacking_positions.add(res_num)
                    
                    # Extract cation-pi positions (exclude axial and pi-stacking)
                    for cp in (interactions_json.get('cation_pi') or []):
                        res_num = cp.get('other_res_id')
                        if res_num and isinstance(res_num, int) and res_num not in axial_positions and res_num not in pi_stacking_positions:
                            cation_pi_positions.add(res_num)
                    
                    # Extract H-bond positions (exclude axial, pi-stacking, and cation-pi)
                    for hb in (interactions_json.get('hydrogen_bonds') or []):
                        res_num = hb.get('other_res_id')  # Fixed: was 'other_res_num'
                        if res_num and isinstance(res_num, int) and res_num not in axial_positions and res_num not in pi_stacking_positions and res_num not in cation_pi_positions:
                            hbond_positions.add(res_num)
            except (json.JSONDecodeError, TypeError, ValueError, AttributeError):
                pass
            
            # Exclude histidine brace positions
            if idx < len(self.histidineBraceTuple):
                brace_positions = set([res for res in self.histidineBraceTuple[idx] if isinstance(res, int)])
                axial_positions.difference_update(brace_positions)
                pi_stacking_positions.difference_update(brace_positions)
                cation_pi_positions.difference_update(brace_positions)
                hbond_positions.difference_update(brace_positions)
            
            self.axialResiduePositions.append(sorted(axial_positions))
            self.piStackingResiduePositions.append(sorted(pi_stacking_positions))
            self.cationPiResiduePositions.append(sorted(cation_pi_positions))
            self.hbondResiduePositions.append(sorted(hbond_positions))

        # After filling self.oneLetterSequence, self.threeLetterSequence, and all property lists:
        if len(self.oneLetterSequence) > 0:
            # Transpose to get columns
            columns = list(zip(*self.oneLetterSequence))
            if self._keep_columns is not None:
                # Trim to an explicit column set (used to restrict plots to the catalytic
                # domain); keep only in-range, sorted, de-duplicated columns.
                keep_indices = sorted(i for i in set(self._keep_columns) if 0 <= i < len(columns))
            else:
                # Indices of columns that are NOT all gaps
                keep_indices = [i for i, col in enumerate(columns) if not all(res == '-' for res in col)]
            # Helper to filter with bounds check
            def filter_by_indices(data_list, indices):
                return [
                    [row[i] for i in indices if i < len(row)]
                    for row in data_list
                ]

            # For StructList (strings), filter as string
            def filter_struct_list(struct_list, indices):
                return [
                    ''.join([row[i] for i in indices if i < len(row)])
                    for row in struct_list
                ]

            # Trim when a keep_columns set is given (catalytic-domain view) or by the normal
            # all-gap rule; no_trim keeps the full alignment (pre-generation).
            if self._keep_columns is not None or not self._no_trim:
                self.oneLetterSequence = filter_by_indices(self.oneLetterSequence, keep_indices)
                self.threeLetterSequence = filter_by_indices(self.threeLetterSequence, keep_indices)
                self.flexibilityList = filter_by_indices(self.flexibilityList, keep_indices)
                for method in self.flexibilityMethods:
                    self.allFlexibilityLists[method] = filter_by_indices(self.allFlexibilityLists[method], keep_indices)
                self.hydrophobicityList = filter_by_indices(self.hydrophobicityList, keep_indices)
                for scale_name in self.hydrophobicityScaleNames:
                    self.allHydrophobicityLists[scale_name] = filter_by_indices(self.allHydrophobicityLists[scale_name], keep_indices)
                if hasattr(self, 'electrostaticsList'):
                    self.electrostaticsList = filter_by_indices(self.electrostaticsList, keep_indices)
                if hasattr(self, 'atomChargesList'):
                    # Filter atomChargesList (list of lists of dicts)
                    self.atomChargesList = [self.atomChargesList[i] for i in keep_indices if i < len(self.atomChargesList)]
                self.bindingSurfaceYN = filter_by_indices(self.bindingSurfaceYN, keep_indices)
                self.loopRegionsYN = filter_by_indices(self.loopRegionsYN, keep_indices)
                self.holeHoppingList = filter_by_indices(self.holeHoppingList, keep_indices)
                self.StructList = filter_struct_list(self.StructList, keep_indices)
            self.residueNumList = [str(i+1) for i in range(len(max(self.oneLetterSequence)))]
            self.sequenceGaps = [(''.join(seq)).split('-') for seq in self.oneLetterSequence]

            for seq in self.StructList:
                if len(seq) < len(self.residueNumList):
                    seq.join('-')


    def alignedResidueNumbers(self, iterator: int) -> list:
        """Original residue number for each alignment column of protein `iterator`,
        or None where the alignment has a gap."""
        protResNums = self.residueNumbers[iterator] if iterator < len(self.residueNumbers) else []
        aligned = []
        seqPos = 0
        for res in self.oneLetterSequence[iterator]:
            if res == '-':
                aligned.append(None)
                continue
            aligned.append(protResNums[seqPos] if seqPos < len(protResNums) else seqPos + 1)
            seqPos += 1
        return aligned

    def _align_pos_template(self):
        """Hovertemplate fragment showing the MSA column (the point's x), shown only
        when comparing several proteins. Non-sequence plots plot on the alignment
        x-axis, so %{x} is the alignment position."""
        return '<br>Alignment position: %{x:.0f}' if getattr(self, 'totalProteinCount', 1) > 1 else ''

    def _domainBarSegments(self, iterator: int, alignedResNums: list) -> list:
        """Translate this protein's domain ranges (auth residue numbering) into
        alignment-column spans for the domain bar.

        Domain boundaries are stored in the structure's own numbering, which is not
        the alignment coordinate, so we locate the columns whose residue falls inside
        each range and take the first/last such column.
        """
        domains = self.domainRanges[iterator] if iterator < len(self.domainRanges) else []
        segments = []
        for dom in domains:
            start_res, end_res = dom.get('start_residue'), dom.get('end_residue')
            if start_res is None or end_res is None:
                continue
            cols = [i + 1 for i, rn in enumerate(alignedResNums)
                    if rn is not None and start_res <= rn <= end_res]
            if not cols:
                continue  # domain not represented in this alignment
            dtype = dom.get('type', 'catalytic')
            name = dom.get('name') or dtype
            segments.append({
                'type': dtype,
                'name': name,
                'label': dom.get('family') or self._domainShortLabel(dtype, dom.get('accession'), name),
                'source': dom.get('source', ''),
                'start_res': start_res,
                'end_res': end_res,
                'start_col': cols[0],
                'end_col': cols[-1],
            })
        return segments

    @staticmethod
    def _domainShortLabel(dtype: str, accession, name: str) -> str:
        """Short label to draw on the sequence overlay."""
        if dtype == 'catalytic':
            return 'Catalytic domain'
        short = DOMAIN_SHORT_LABELS.get(accession)
        if short:
            return short
        m = re.search(r'CBM\s*(\d+(?:/\d+)?)', name or '', re.I)
        return f'CBM{m.group(1)}' if m else 'CBM'

    def makeHydrophobColorScale(self, hydrophobDict: dict, type: str) -> dict or list:

        normalizedVals = []
        hydrophobRGBlist = []
        notNoneDict = {k:hydrophobDict[k] for k in hydrophobDict.keys() if k != '-'}

        sortedDict = {k: v for k, v in sorted(notNoneDict.items(), key=lambda item: item[1])}


        maxVal = max(notNoneDict.values())
        minVal = min(notNoneDict.values())

        for key in sortedDict.keys():
            normalizedVal = round((hydrophobDict[key] - minVal) / (maxVal - minVal), 2)


            if type == 'colorscale' and normalizedVal in normalizedVals:
                continue

            normalizedVals.append(normalizedVal)

            if normalizedVal <= 0.5:
                ratio = normalizedVal * 2
                hydrophobRGB = [int(255 * ratio), int(255 * ratio), 255]

            else:

                ratio = (normalizedVal - 0.5) / 0.5
                hydrophobRGB = [255, int(255 * (1 - ratio)), int(255 * (1 - ratio))]

            hydrophobRGBlist.append(hydrophobRGB)

        if type == 'colorscale':
            return [[n, f'rgb({rgb[0]}, {rgb[1]}, {rgb[2]})'] for n, rgb in zip(normalizedVals, hydrophobRGBlist)]

        elif type == 'aaRGB':
            return {key:rgb for key, rgb in zip(sortedDict.keys(), hydrophobRGBlist)}


    def makeHydrophobicityPlots(self, figure:plotly.graph_objs.Figure, row:int=None, col:int=None, scale_name:str=None) -> None:
        '''Plots hydrophobicity as a heatmap for each protein.'''

        iterator = 0
        annotationList = []

        if scale_name and scale_name in self.allHydrophobicityScales:
            hydrophobDict = self.allHydrophobicityScales[scale_name]
            hydrophobicityData = self.allHydrophobicityLists[scale_name]
        else:
            hydrophobDict = self.aaHydrophobicity
            hydrophobicityData = self.hydrophobicityList

        colorScale = self.makeHydrophobColorScale(hydrophobDict, 'colorscale')

        # Calculate zmin/zmax from the scale values (excluding None/gaps)
        scale_values = [v for v in hydrophobDict.values() if v is not None]
        zmin_val = min(scale_values)
        zmax_val = max(scale_values)

        for protein in hydrophobicityData:


            customdata = []
            hover_colors = []
            alignedResNums = self.alignedResidueNumbers(iterator)
            for idx, letter in enumerate(self.threeLetterSequence[iterator]):
                if letter != '---':
                    residue_number = alignedResNums[idx]
                    customdata.append(f'{letter}{residue_number}<br>Hydrophobicity: {hydrophobicityData[iterator][idx]:.2f} ')
                    
                    # Get hydrophobicity value
                    hydro_val = hydrophobicityData[iterator][idx]
                    
                    if hydro_val is not None:
                        # Normalize to 0-1 range
                        normalized = max(0, min(1, (hydro_val - zmin_val) / (zmax_val - zmin_val))) if zmax_val != zmin_val else 0.5
                        
                        # Find color in colorScale
                        color_found = False
                        for i in range(len(colorScale) - 1):
                            pos1, color1 = colorScale[i]
                            pos2, color2 = colorScale[i + 1]
                            
                            if pos1 <= normalized <= pos2:
                                # Interpolate
                                ratio = (normalized - pos1) / (pos2 - pos1) if pos2 != pos1 else 0
                                
                                # Parse RGB - handle different formats
                                try:
                                    rgb1 = [int(x.strip()) for x in color1.replace('rgb(', '').replace(')', '').split(',')]
                                    rgb2 = [int(x.strip()) for x in color2.replace('rgb(', '').replace(')', '').split(',')]
                                    
                                    r = int(rgb1[0] + (rgb2[0] - rgb1[0]) * ratio)
                                    g = int(rgb1[1] + (rgb2[1] - rgb1[1]) * ratio)
                                    b = int(rgb1[2] + (rgb2[2] - rgb1[2]) * ratio)
                                    
                                    hover_colors.append(f'rgba({r},{g},{b},0.7)')
                                    color_found = True
                                    break
                                except Exception as e:
                                    print(f"Error parsing colors: {e}, color1={color1}, color2={color2}")
                                    hover_colors.append('rgba(176,190,197,0.7)')
                                    color_found = True
                                    break
                        
                        if not color_found:
                            # Use last color
                            try:
                                rgb = [int(x.strip()) for x in colorScale[-1][1].replace('rgb(', '').replace(')', '').split(',')]
                                hover_colors.append(f'rgba({rgb[0]},{rgb[1]},{rgb[2]},0.7)')
                            except:
                                hover_colors.append('rgba(176,190,197,0.7)')
                    else:
                        hover_colors.append('rgba(176,190,197,0.7)')
                else:
                    customdata.append('')
                    hover_colors.append('rgba(176,190,197,0.7)')


            figure.add_trace(go.Heatmap(
                x=[int(x) for x in self.residueNumList],
                y=np.linspace(0, 0, len(self.threeLetterSequence[iterator]), dtype='int').tolist(),
                z=hydrophobicityData[iterator],
                showscale=False,
                colorscale=colorScale,
                zmid=(zmin_val + zmax_val) / 2,
                zmax=zmax_val,
                zmin=zmin_val,
                opacity=0.6,
                hoverinfo='skip',  # Disable hover on heatmap, we'll use scatter instead
                name='hydrophobicity',
                showlegend=False,
                colorbar={'xref': 'container', 'title': 'hydrophobicity', 'yref': 'container', 'orientation': 'v',
                          'y': 1000}
            ), col=col, row=row+iterator)
            
            # Add invisible scatter trace on top for hover with custom colors
            # Filter out gap positions (empty customdata)
            x_no_gaps = []
            y_no_gaps = []
            customdata_no_gaps = []
            hover_colors_no_gaps = []
            for idx, cd in enumerate(customdata):
                if cd != '':  # Skip gaps
                    x_no_gaps.append(int(self.residueNumList[idx]))
                    y_no_gaps.append(0)
                    customdata_no_gaps.append(cd)
                    hover_colors_no_gaps.append(hover_colors[idx])
            
            figure.add_trace(go.Scatter(
                x=x_no_gaps,
                y=y_no_gaps,
                mode='markers',
                marker=dict(size=10, opacity=0),  # Invisible markers
                customdata=customdata_no_gaps,
                hovertemplate=' %{customdata}' + self._align_pos_template() + '<extra></extra>',
                hoverlabel=dict(bgcolor=hover_colors_no_gaps),
                showlegend=False,
                hoverinfo='text'
            ), col=col, row=row+iterator)



            for num in self.residueNumList:
                xref_axis = 'x'
                y_axis_idx = (row + iterator) if (row is not None) else (iterator + 1)
                y_axis_name = 'y' if y_axis_idx == 1 else f'y{y_axis_idx}'
                yref_str = f'{y_axis_name} domain'
                
                annotationList.append({'font': {'color': 'black', 'size': 11, 'weight': 'bold'},
                                'font': {'color': 'black', 'size': 11, 'weight': 'bold'},
                                'showarrow': False,
                                'text': self.oneLetterSequence[iterator][int(num) - 1],
                                'x': int(num) ,
                                'xanchor': 'center',
                                'xref': xref_axis,
                                'y': 0.5,
                                'yanchor': 'middle',
                                'yref': yref_str
                                })

            iterator += 1

        figure.update_layout(annotations= list(figure['layout']['annotations']) + annotationList)

    def makeChargePlots(self, figure:plotly.graph_objs.Figure, row:int=None, col:int=None) -> None:
        '''Plots charge as a heatmap for each protein. Uses fixed scale from -2 to 2, where 0 is white.'''
        
        iterator = 0
        annotationList = []
        
        # Create color scale for charge: Red (negative) -> White (neutral) -> Blue (positive)
        if not hasattr(self, 'electrostaticsList') or not self.electrostaticsList:
            return
        
        # Use fixed scale from -2 to 2, where 0 corresponds to white
        charge_min = -2.0
        charge_max = 2.0
        charge_mid = 0.0
        
        # Create color scale: Red -> White -> Blue
        chargeColorScale = [
            [0.0, 'rgb(255, 0, 0)'],      # Red for charge = -2
            [0.5, 'rgb(255, 255, 255)'],  # White for charge = 0
            [1.0, 'rgb(0, 0, 255)']       # Blue for charge = 2
        ]
        
        # Create charge color function for hover colors
        chargeColorFunc = self.makeChargeColorScale()
        
        for protein_charges in self.electrostaticsList:
            customdata = []
            hover_colors = []
            alignedResNums = self.alignedResidueNumbers(iterator)

            for idx, letter in enumerate(self.threeLetterSequence[iterator]):
                if letter != '---':
                    residue_number = alignedResNums[idx]
                    charge_val = protein_charges[idx] if idx < len(protein_charges) else None
                    if charge_val is not None:
                        # Per-residue formal charge is integer-valued — show it without decimals
                        customdata.append(f'{letter}{residue_number}<br>Charge: {int(round(charge_val))} e')

                        # Use the charge color function to get RGB
                        chargeRGB = chargeColorFunc(charge_val)
                        hover_colors.append(f'rgba({chargeRGB[0]},{chargeRGB[1]},{chargeRGB[2]},0.7)')
                    else:
                        customdata.append(f'{letter}{residue_number}<br>Charge: N/A')
                        hover_colors.append('rgba(176,190,197,0.7)')
                else:
                    customdata.append('')
                    hover_colors.append('rgba(176,190,197,0.7)')
            
            # Normalize charges for heatmap display using fixed scale -2 to 2
            normalized_charges = []
            for idx_c, charge_val in enumerate(protein_charges):
                # Use None for gaps (transparent)
                if idx_c < len(self.threeLetterSequence[iterator]) and self.threeLetterSequence[iterator][idx_c] == '---':
                    normalized_charges.append(None)
                elif charge_val is not None:
                    # Clamp charge value to the scale range
                    charge_val = max(charge_min, min(charge_max, charge_val))

                    # Normalize to 0-1 range where 0.5 corresponds to charge=0
                    if charge_val <= charge_mid:
                        # Negative charges: map from -2 (0.0) to 0 (0.5)
                        normalized = 0.5 * (charge_val - charge_min) / (charge_mid - charge_min)
                    else:
                        # Positive charges: map from 0 (0.5) to 2 (1.0)
                        normalized = 0.5 + 0.5 * (charge_val - charge_mid) / (charge_max - charge_mid)

                    normalized_charges.append(max(0.0, min(1.0, normalized)))
                else:
                    normalized_charges.append(None)  # Transparent for missing data
            
            figure.add_trace(go.Heatmap(
                x=[int(x) for x in self.residueNumList],
                y=np.linspace(0, 0, len(self.threeLetterSequence[iterator]), dtype='int').tolist(),
                z=normalized_charges,
                showscale=False,
                colorscale=chargeColorScale,
                zmid=0.5,
                zmax=1.0,
                zmin=0.0,
                opacity=0.6,
                hoverinfo='skip',
                name='charge',
                showlegend=False,
            ), col=col, row=row+iterator)
            
            # Add invisible scatter trace for hover
            x_no_gaps = []
            y_no_gaps = []
            customdata_no_gaps = []
            hover_colors_no_gaps = []
            for idx, cd in enumerate(customdata):
                if cd != '':
                    x_no_gaps.append(int(self.residueNumList[idx]))
                    y_no_gaps.append(0)
                    customdata_no_gaps.append(cd)
                    hover_colors_no_gaps.append(hover_colors[idx])
            
            figure.add_trace(go.Scatter(
                x=x_no_gaps,
                y=y_no_gaps,
                mode='markers',
                marker=dict(size=10, opacity=0),
                customdata=customdata_no_gaps,
                hovertemplate=' %{customdata}' + self._align_pos_template() + '<extra></extra>',
                hoverlabel=dict(bgcolor=hover_colors_no_gaps),
                showlegend=False,
                hoverinfo='text'
            ), col=col, row=row+iterator)
            
            # Add residue letter annotations
            for num in self.residueNumList:
                xref_axis = 'x'
                y_axis_idx = (row + iterator) if (row is not None) else (iterator + 1)
                y_axis_name = 'y' if y_axis_idx == 1 else f'y{y_axis_idx}'
                yref_str = f'{y_axis_name} domain'
                
                annotationList.append({
                    'font': {'color': 'black', 'size': 11, 'weight': 'bold'},
                    'showarrow': False,
                    'text': self.oneLetterSequence[iterator][int(num) - 1],
                    'x': int(num),
                    'xanchor': 'center',
                    'xref': xref_axis,
                    'y': 0.5,
                    'yanchor': 'middle',
                    'yref': yref_str
                })
            
            iterator += 1
        
        figure.update_layout(annotations=list(figure['layout']['annotations']) + annotationList if figure['layout']['annotations'] else annotationList)

    def makeSequencePlots(self, figure:plotly.graph_objs.Figure, row:int=None, col:int=None) -> None:
        '''Creates a heatmap for the aligned sequence and adds residue annotations.'''
        
        iterator = 0
        annotationList = []
        self.sequenceDict = {
            'positive': [], 
            'negative': [], 
            'aromatic': [], 
            'sulphur': [], 
            'brace': [],
            'Axial position (to Cu)': [],  # Dark orange
            'Pi-pi stacking': [],  # Medium orange
            'Cation-pi': [],  # Light golden orange
            'H-bond': [],  # Light orange
            SEQ_MOTIF: [],
            SEQ_DISTAL: [],
            SEQ_BRACE_MOTIF: [],
        }
        residueIndex = []

        protein_ids = getattr(self, 'protId', None)
        if not protein_ids:
            protein_ids = [str(i) for i in range(len(self.threeLetterSequence))]

        for seq in self.threeLetterSequence:
            # Build colors and hover data for each residue
            colors_rgb = []  # Will store RGB tuples for each position
            hover_data = []
            hover_colors = []
            
            alignedResNums = self.alignedResidueNumbers(iterator)
            residueIndex.append([])

            for val in self.sequenceDict.values():
                val.append([])

            combined_second_sphere = set()
            if hasattr(self, 'secondSphereCombinedResidues') and len(self.secondSphereCombinedResidues) > iterator:
                combined_second_sphere = set(self.secondSphereCombinedResidues[iterator])
            
            # Get categorization sets for labeling
            axial_set = set()
            pi_stacking_set = set()
            cation_pi_set = set()
            hbond_set = set()
            if hasattr(self, 'axialResiduePositions') and len(self.axialResiduePositions) > iterator:
                axial_set = set(self.axialResiduePositions[iterator])
            if hasattr(self, 'piStackingResiduePositions') and len(self.piStackingResiduePositions) > iterator:
                pi_stacking_set = set(self.piStackingResiduePositions[iterator])
            if hasattr(self, 'cationPiResiduePositions') and len(self.cationPiResiduePositions) > iterator:
                cation_pi_set = set(self.cationPiResiduePositions[iterator])
            if hasattr(self, 'hbondResiduePositions') and len(self.hbondResiduePositions) > iterator:
                hbond_set = set(self.hbondResiduePositions[iterator])
            # motif residues other than the aromatic: residue number -> hover text
            site = (self.secondSphereMotifs[iterator]
                    if iterator < len(getattr(self, 'secondSphereMotifs', [])) else second_sphere(None))
            motif_text = {r['residue']: motif_residue_text(site['motif'], r) for r in site['residues']}
            motif_aromatic = {r['residue'] for r in site['residues'] if r.get('role') == 'aromatic'}
            distal = site['distal_residue'] or {}
            # Gly/Pro/Arg of the motifs after the brace histidines: residue number -> hover text
            brace = (self.braceMotifs[iterator]
                     if iterator < len(getattr(self, 'braceMotifs', [])) else brace_motifs(None))
            brace_text = {r['residue']: brace_motif_residue_text(which, m['motif'], r)
                          for which, m in brace.items() for r in m['motif_residues']}
            
            # Debug: print categorization (commented out)
            # if iterator == 0:  # Only print for first protein
            #     print(f"DEBUG - Protein {iterator}:")
            #     print(f"  Axial positions: {axial_set}")
            #     print(f"  Pi-stacking positions: {pi_stacking_set}")
            #     print(f"  H-bond positions: {hbond_set}")
            #     print(f"  Combined second sphere: {combined_second_sphere}")
            #     print(f"  Uncategorized residues: {combined_second_sphere - axial_set - pi_stacking_set - hbond_set}")

            for num in self.residueNumList:
                residue_letter = self.oneLetterSequence[iterator][int(num) - 1]
                # Original residue number from the structure file (None at gaps)
                residue = alignedResNums[int(num) - 1]
                residueIndex[iterator].append(residue if residue_letter != '-' else 'NaN')

                if residue_letter == '-':
                    text = ''
                    rgb_color = None  # Transparent for gaps
                    rgba_hover = 'rgba(176,190,197,0.7)'
                else:
                    text = f'{self.threeLetterSequence[iterator][int(num) -1]}{residue}'
                    rgb_color = None  # Transparent for non-highlighted residues
                    rgba_hover = 'rgba(176,190,197,0.7)'

                    if residue in self.histidineBraceTuple[iterator]:
                        text += '<br>Histidine brace '
                        rgb_color = (186, 3, 186)
                        rgba_hover = 'rgba(186,3,186,0.7)'
                        self.sequenceDict['brace'][iterator].append({'residue_number': residue, 'color': {'r':186, 'g': 3, 'b': 186}})

                    elif residue in axial_set:
                        text += '<br>Second sphere: Axial position (to Cu) '
                        if residue in motif_aromatic:
                            text += f"<br>Second-sphere motif {site['motif']}"
                        rgb_color = (255, 140, 0)  # Dark Orange
                        rgba_hover = 'rgba(255,140,0,0.7)'
                        self.sequenceDict['Axial position (to Cu)'][iterator].append({'residue_number': residue, 'color': {'r': 255, 'g': 140, 'b': 0}})

                    elif residue in motif_text:
                        text += f'<br>{motif_text[residue]} '
                        rgb_color = SEQ_MOTIF_COLOR
                        rgba_hover = 'rgba({},{},{},0.85)'.format(*SEQ_MOTIF_COLOR)
                        self.sequenceDict[SEQ_MOTIF][iterator].append({'residue_number': residue, 'color': dict(zip('rgb', SEQ_MOTIF_COLOR))})

                    elif residue == distal.get('residue'):
                        text += (f"<br>Second sphere: shielding residue on the copper's solvent side "
                                 f"({distal.get('atom')} {distal.get('distance')} Å from the copper) ")
                        rgb_color = SEQ_DISTAL_COLOR
                        rgba_hover = 'rgba({},{},{},0.85)'.format(*SEQ_DISTAL_COLOR)
                        self.sequenceDict[SEQ_DISTAL][iterator].append({'residue_number': residue, 'color': dict(zip('rgb', SEQ_DISTAL_COLOR))})

                    elif residue in brace_text:
                        text += f'<br>{brace_text[residue]} '
                        rgb_color = SEQ_BRACE_MOTIF_COLOR
                        rgba_hover = 'rgba({},{},{},0.85)'.format(*SEQ_BRACE_MOTIF_COLOR)
                        self.sequenceDict[SEQ_BRACE_MOTIF][iterator].append({'residue_number': residue, 'color': dict(zip('rgb', SEQ_BRACE_MOTIF_COLOR))})
                    
                    elif residue in pi_stacking_set:
                        text += '<br>Second sphere: Pi-pi stacking '
                        rgb_color = (255, 165, 0)  # Medium Orange
                        rgba_hover = 'rgba(255,165,0,0.7)'
                        self.sequenceDict['Pi-pi stacking'][iterator].append({'residue_number': residue, 'color': {'r': 255, 'g': 165, 'b': 0}})
                    
                    elif residue in cation_pi_set:
                        text += '<br>Second sphere: Cation-pi interaction '
                        rgb_color = (255, 175, 50)  # Light Golden Orange
                        rgba_hover = 'rgba(255,175,50,0.7)'
                        self.sequenceDict['Cation-pi'][iterator].append({'residue_number': residue, 'color': {'r': 255, 'g': 175, 'b': 50}})
                    
                    elif residue in hbond_set:
                        text += '<br>Second sphere: H-bond interaction '
                        rgb_color = (255, 200, 100)  # Light Orange
                        rgba_hover = 'rgba(255,200,100,0.7)'
                        self.sequenceDict['H-bond'][iterator].append({'residue_number': residue, 'color': {'r': 255, 'g': 200, 'b': 100}})

                    elif residue_letter in ['R', 'H', 'K'] and residue not in self.histidineBraceTuple[iterator] and residue not in axial_set and residue not in pi_stacking_set and residue not in cation_pi_set and residue not in hbond_set:
                        text += '<br>Positively charged '
                        rgb_color = (255, 180, 180)  # Pastel red/pink
                        rgba_hover = 'rgba(255,180,180,0.7)'
                        self.sequenceDict['positive'][iterator].append({'residue_number': residue, 'color': {'r':255, 'g': 180, 'b': 180}})

                    elif residue_letter in ['D', 'E'] and residue not in axial_set and residue not in pi_stacking_set and residue not in cation_pi_set and residue not in hbond_set:
                        text += '<br>Negatively charged '
                        rgb_color = (180, 180, 255)  # Pastel blue
                        rgba_hover = 'rgba(180,180,255,0.7)'
                        self.sequenceDict['negative'][iterator].append({'residue_number': residue, 'color': {'r':180, 'g': 180, 'b': 255}})

                    elif residue_letter in ['F', 'Y', 'W'] and residue not in axial_set and residue not in pi_stacking_set and residue not in cation_pi_set and residue not in hbond_set:
                        text += '<br>Aromatic '
                        rgb_color = (180, 255, 180)  # Pastel green
                        rgba_hover = 'rgba(180,255,180,0.7)'
                        self.sequenceDict['aromatic'][iterator].append({'residue_number': residue, 'color': {'r':180, 'g': 255, 'b': 180}})

                    elif residue_letter in ['M', 'C'] and residue not in axial_set and residue not in pi_stacking_set and residue not in cation_pi_set and residue not in hbond_set:
                        text += '<br>Contains sulphur '
                        rgb_color = (255, 255, 180)  # Pastel yellow
                        rgba_hover = 'rgba(255,255,180,0.7)'
                        self.sequenceDict['sulphur'][iterator].append({'residue_number': residue, 'color': {'r':255, 'g': 255, 'b': 180}})

                    # Alignment column — only informative when comparing several proteins,
                    # where the same column lines up equivalent positions across structures.
                    if getattr(self, 'totalProteinCount', 1) > 1:
                        text += f'<br>Alignment position: {int(num)}'

                colors_rgb.append(rgb_color)
                hover_data.append(text)
                hover_colors.append(rgba_hover)
                
                # Still add text annotation for the letter
                xref_axis = 'x'
                y_axis_idx = (row + iterator) if (row is not None) else (iterator + 1)
                y_axis_name = 'y' if y_axis_idx == 1 else f'y{y_axis_idx}'

                annotationList.append({
                    'font': {'color': 'black', 'size': 11, 'weight': 'bold'},
                    'showarrow': False,
                    'text': self.oneLetterSequence[iterator][int(num) - 1],
                    'x': int(num),
                    'xanchor': 'center',
                    'xref': xref_axis,
                    'y': 0.,
                    'yanchor': 'middle',
                    'yref': y_axis_name,
                })

            # Create colorscale from RGB colors, using None for transparent positions
            # Build z-values: sequential index for colored positions, None for transparent
            z_values = []
            colored_rgb = []  # Only non-None colors
            for rgb in colors_rgb:
                if rgb is None:
                    z_values.append(None)
                else:
                    z_values.append(len(colored_rgb))
                    colored_rgb.append(rgb)

            if len(colored_rgb) >= 2:
                colorscale_seq = [[i/(len(colored_rgb)-1), f'rgb({rgb[0]},{rgb[1]},{rgb[2]})']
                                 for i, rgb in enumerate(colored_rgb)]
            elif len(colored_rgb) == 1:
                c = colored_rgb[0]
                colorscale_seq = [[0.0, f'rgb({c[0]},{c[1]},{c[2]})'], [1.0, f'rgb({c[0]},{c[1]},{c[2]})']]
            else:
                colorscale_seq = [[0.0, 'rgb(240,240,240)'], [1.0, 'rgb(240,240,240)']]

            # Add colored heatmap
            n_colors = max(len(colored_rgb) - 1, 1)
            figure.add_trace(go.Heatmap(
                x=[int(x) for x in self.residueNumList],
                y=np.linspace(0, 0, len(seq), dtype='int').tolist(),
                z=z_values,
                showscale=False,
                colorscale=colorscale_seq,
                zmin=0,
                zmax=n_colors,
                hoverinfo='skip',
                showlegend=False,
                name=protein_ids[iterator],
                meta={'protein_id': protein_ids[iterator]},
            ), col=col, row=row+iterator)

            # ── Domain bar (catalytic domain / CBMs) above the sequence ──────────
            # Domain ranges are in the structure's own auth numbering; alignedResNums
            # maps each alignment column to its auth residue number, so we translate
            # by scanning for the columns whose residue falls inside each range.
            #
            # Each segment is two traces: a filled band (the visual) and an invisible
            # marker row (the hover target). Splitting them is what makes the tooltip
            # work — a single filled trace routes hover through fill logic and swallows
            # it. The band is a *trace* (not a layout shape) so the frontend can darken
            # its `.js-fill` path directly on hover, without a Plotly redraw that would
            # wipe the tooltip.
            seg_index = 0
            for seg in self._domainBarSegments(iterator, alignedResNums):
                x0, x1 = seg['start_col'] - 0.5, seg['end_col'] + 0.5
                base_color = DOMAIN_BAR_COLORS.get(seg['type'], 'rgba(150,150,150,0.40)')
                hover_color = DOMAIN_BAR_HOVER_COLORS.get(seg['type'], 'rgba(150,150,150,0.80)')

                # 1. visible band — hoverinfo skip so it never competes for the hover
                figure.add_trace(go.Scatter(
                    x=[x0, x1, x1, x0, x0],
                    y=[DOMAIN_BAR_Y0, DOMAIN_BAR_Y0, DOMAIN_BAR_Y1, DOMAIN_BAR_Y1, DOMAIN_BAR_Y0],
                    mode='lines', fill='toself', fillcolor=base_color,
                    line=dict(width=0), hoverinfo='skip', name='', showlegend=False,
                ), col=col, row=row + iterator)

                # 2. invisible hover targets — exact recipe as the residue hover scatter:
                # markers, customdata + hovertemplate with a leading space, NO hoverinfo
                # (hoverinfo='text' alongside hovertemplate renders an empty label).
                mid_y = (DOMAIN_BAR_Y0 + DOMAIN_BAR_Y1) / 2
                hover_x = [float(c) for c in range(seg['start_col'], seg['end_col'] + 1)]
                tooltip = (f"<b>{seg['label']}</b><br>{seg['name']}"
                           f"<br>Residues {seg['start_res']}–{seg['end_res']}"
                           f"<br>Source: {seg['source']}")
                figure.add_trace(go.Scatter(
                    x=hover_x,
                    y=[mid_y] * len(hover_x),
                    mode='markers',
                    marker=dict(size=14, opacity=0),
                    customdata=[tooltip] * len(hover_x),
                    hovertemplate=' %{customdata}<extra></extra>',
                    hoverlabel=dict(bgcolor=hover_color, font=dict(color='white', size=HOVER_FONT_SIZE)),
                    # Tagged so the frontend can (a) skip the residue-hover handler and
                    # the cursor and (b) darken this segment's band by index. segIndex is
                    # the position among this plot's `.js-fill` band paths (bands are the
                    # only filled traces, added in this same order).
                    meta={'domainBar': True, 'segIndex': seg_index,
                          'baseColor': base_color, 'hoverColor': hover_color,
                          'startRes': seg['start_res'], 'endRes': seg['end_res']},
                    showlegend=False,
                ), col=col, row=row + iterator)
                seg_index += 1

            # Residue row occupies the bottom of the plot (data -0.5..0.5); the domain
            # bar sits above it. The range top leaves headroom for the bar.
            figure.update_yaxes(range=[-0.5, DOMAIN_BAR_YMAX], col=col, row=row + iterator)

            # Add invisible scatter for hover with custom colors
            # Filter out gap positions (empty hover_data)
            x_no_gaps = []
            y_no_gaps = []
            customdata_no_gaps = []
            hover_colors_no_gaps = []
            for idx, hd in enumerate(hover_data):
                if hd != '':  # Skip gaps
                    x_no_gaps.append(int(self.residueNumList[idx]))
                    y_no_gaps.append(0)
                    customdata_no_gaps.append(hd)
                    hover_colors_no_gaps.append(hover_colors[idx])
            
            figure.add_trace(go.Scatter(
                x=x_no_gaps,
                y=y_no_gaps,
                mode='markers',
                marker=dict(size=10, opacity=0),
                customdata=customdata_no_gaps,
                hovertemplate='%{customdata}<extra></extra>',
                hoverlabel=dict(bgcolor=hover_colors_no_gaps),
                showlegend=False,
                hoverinfo='text'
            ), col=col, row=row+iterator)

            iterator += 1

        existing_annotations = figure.layout.annotations if figure.layout.annotations else []
        figure.update_layout(
            annotations=list(existing_annotations) + annotationList,
            autosize=True
        )


    def makeFlexibilityPlot(self, figure:go.Figure, proteinIndex:int=0, row:int=1, col:int=1, showLegend:bool=True, method_name:str=None) -> None:
        '''Plots flexibility for each protein as a line plot. `method_name` selects which
        flexibility measure to show (cabsflex/gnm/md_rmsf/bfactor/plddt); defaults to the
        default method. Colour gradient + y-axis + hover use that method's own range/units.'''

        # Resolve the method: the requested one if this comparison has it, else the default.
        method = method_name if (method_name in self.allFlexibilityLists) else self.defaultFlexMethod
        method_lists = self.allFlexibilityLists.get(method)
        flex_min = self.flexMinByMethod.get(method, self.minFlexVal)
        flex_max = self.flexMaxByMethod.get(method, self.maxFlexVal)
        meta = FLEX_METHOD_META.get(method, {'hover': 'Flexibility', 'unit': ''})
        hov, unit = meta['hover'], meta['unit']
        invert = meta.get('invert', False)
        use_alphafold = meta.get('scheme') == 'alphafold'

        def _flex_color(value):
            # pLDDT: official AlphaFold confidence bands; everything else: red gradient
            # scaled to this method's own range (red = more flexible).
            if use_alphafold:
                r, g, b = _plddt_rgb(value)
                return f'rgb({r}, {g}, {b})'
            normalized = (value - flex_min) / ((flex_max - flex_min) + 1e-9)
            cval = (1 - normalized) if invert else normalized
            c = int(200 * (1 - cval))
            return f'rgb(255, {c}, {c})'

        customData = []
        alignedResNums = self.alignedResidueNumbers(proteinIndex)
        for idx, residue in enumerate(self.threeLetterSequence[proteinIndex]):
            if residue != '---':
                customData.append(residue + str(alignedResNums[idx]))
            else:
                customData.append(None)



        # flexibilityList already has None values for gaps, aligned to sequence positions!
        # Just like sec_struct_plot.py uses the sequence directly
        y_values = (method_lists[proteinIndex] if method_lists and proteinIndex < len(method_lists)
                    else self.flexibilityList[proteinIndex])

        # Residue numbering = the structure's original (auth) residue numbers, so the hover
        # matches the Sequence/Charge plots and the 3D viewer (not a 1-based count).
        residue_numbers = self.alignedResidueNumbers(proteinIndex)

        # Track which positions have been connected by lines
        connected_positions = set()
        
        # Plot lines between consecutive non-gap positions
        for i in range(len(y_values) - 1):
            y1 = y_values[i]
            y2 = y_values[i+1]
            
            # Skip if either is a gap
            if y1 is None or y2 is None:
                continue
            
            # Mark these positions as connected
            connected_positions.add(i)
            connected_positions.add(i+1)
            
            # X positions are alignment positions (i+1 for 1-based indexing)
            x1 = i + 1
            x2 = i + 2
            
            avg_flex = (y1 + y2) / 2
            color = _flex_color(avg_flex)
            
            # Get real residue numbers for hover
            res_num_1 = residue_numbers[i]
            res_num_2 = residue_numbers[i+1]
            
            hover_data = [
                f'{self.threeLetterSequence[proteinIndex][i]}{res_num_1}<br>{hov}: {y1:.2f}{unit}',
                f'{self.threeLetterSequence[proteinIndex][i+1]}{res_num_2}<br>{hov}: {y2:.2f}{unit}'
            ]
            
            figure.add_trace(go.Scatter(
                x=[x1, x2],
                y=[y1, y2],
                mode='lines',
                line=dict(color=color, width=2.5, shape='spline'),
                customdata=hover_data,
                hovertemplate='%{customdata}' + self._align_pos_template() + '<extra></extra>',
                hoverlabel=dict(bgcolor=color),
                showlegend=False,
            ), row=row, col=col)
        
        # Add markers for isolated points (not connected by lines)
        # This ensures single values surrounded by gaps are visible and have hover
        isolated_x = []
        isolated_y = []
        isolated_hover = []
        isolated_colors = []
        
        for i in range(len(y_values)):
            y_val = y_values[i]
            
            # Skip gaps and already connected points
            if y_val is None or i in connected_positions:
                continue
            
            # This is an isolated point
            x_pos = i + 1
            res_num = residue_numbers[i]
            
            isolated_x.append(x_pos)
            isolated_y.append(y_val)
            isolated_hover.append(f'{self.threeLetterSequence[proteinIndex][i]}{res_num}<br>{hov}: {y_val:.2f}{unit}')
            
            # Calculate color for this point
            color = _flex_color(y_val)
            isolated_colors.append(color)
        
        # Add isolated points as markers
        if isolated_x:
            figure.add_trace(go.Scatter(
                x=isolated_x,
                y=isolated_y,
                mode='markers',
                marker=dict(
                    size=5,
                    color=isolated_colors,
                    line=dict(width=0.5, color='white')
                ),
                customdata=isolated_hover,
                hovertemplate='%{customdata}' + self._align_pos_template() + '<extra></extra>',
                hoverlabel=dict(bgcolor=isolated_colors),
                showlegend=False,
            ), row=row, col=col)
        
        # Set x-axis to start from 0 and y-axis to be consistent across all proteins
        # Use max length across ALL proteins for consistent x-axis
        maxlen = max([len(seq) for seq in self.oneLetterSequence])
        
        # Calculate y-axis range with padding (same for all proteins, this method's range)
        y_padding = (flex_max - flex_min) * 0.05 if (flex_max is not None and flex_min is not None) else 0.5
        y_min = 0  # Start from 0
        y_max = (flex_max + y_padding) if flex_max else 5
        
        # Update the correct axes based on row position
        if row == 1:
            figure.update_layout(
                xaxis={'range': [1, maxlen + 1], 'dtick': 5, 'tick0': 5},
                yaxis={'range': [y_min, y_max]}
            )
        else:
            axis_x = f'xaxis{row}'
            axis_y = f'yaxis{row}'
            figure.update_layout(**{
                axis_x: {'range': [1, maxlen + 1], 'dtick': 5, 'tick0': 5},
                axis_y: {'range': [y_min, y_max]}
            })




    def makeSASAPlots(self, figure:plotly.graph_objs.Figure, row:int=None, col:int=None) -> None:
        '''Plots relative SASA (surface accessibility) for each protein.'''

        iterator = 0

        if self.sequenceGaps == []:
            figure.add_trace(go.Scatter(
                x=list(map(lambda x: int(x)), self.residueNumList),
                y=self.relSASAList[iterator],
                customdata=self.threeLetterSequence[iterator],
                hovertemplate='%{customdata}%{x}<br>Relative Surface Accessibility: %{y:.2f} ' + self._align_pos_template() + '<extra></extra>',
                showlegend=False,
                line={'color': 'rgba(0,0,0,0.5)', 'shape': 'spline', 'width': 1},
            ), row=row, col=col)

            figure.add_trace(go.Scatter(
                x=list(map(lambda x: int(x)), self.residueNumList),
                y=np.linspace(25, 25, len(self.oneLetterSequence[iterator])),
                customdata=self.threeLetterSequence[iterator],
                showlegend=False,
                line={'width': 0},
                fill='tonexty',
                fillgradient={
                    'colorscale': [(0, 'rgba(95, 140, 200, 0.7)'), (0.24, 'rgba(95, 140, 200, 0.7)'), (0.25, 'rgba(208, 100, 88, 0.7)'),
                                   (1, 'rgba(208, 100, 88, 0.7)')], 'type': 'vertical'},
                hoverinfo='skip',
            ), row=row, col=col)


        else:

            for seq in self.threeLetterSequence:
                residue = 0
                alignment = 0
                rawNoGaps = ''.join(self.threeLetterSequence[iterator])
                rawNoGaps = rawNoGaps.split('---')
                data = "".join(rawNoGaps)
                noGapsThreeLetterSequence = textwrap.wrap(data, 3)
                # Ungapped residues map 1:1 onto the structure's own residue numbers
                protResNums = self.residueNumbers[iterator] if iterator < len(self.residueNumbers) else []
                customdata = []

                for pos, aa in enumerate(noGapsThreeLetterSequence):
                    resNum = protResNums[pos] if pos < len(protResNums) else pos + 1
                    customdata.append(aa + str(resNum))

                for gap in self.sequenceGaps[iterator]:
                    lenght = len(gap)

                    if lenght > 0:
                        values = [v for v in self.relSASAList[iterator][residue:residue+lenght] if v is not None]
                        if values:
                            maxValue = max(values)
                            minValue = min(values)


                            if maxValue != minValue:
                                if maxValue > 25 and minValue < 25:
                                    # Values span across 25 threshold - calculate gradient position
                                    treshold = (25-minValue)/(maxValue-minValue)
                                elif minValue >= 25:
                                    # All values above 25 - should be all red
                                    treshold = 0  # Blue takes 0%, red takes 100%
                                elif maxValue <= 25:
                                    # All values below 25 - should be all blue
                                    treshold = 1  # Blue takes 100%, red takes 0%
                                else:
                                    treshold = 0.25
                            else:
                                # All values are the same
                                if minValue >= 25:
                                    treshold = 0  # All red
                                else:
                                    treshold = 1  # All blue

                            # Prepare y-data and conditional hover colors
                            y_data = self.relSASAList[iterator][residue:residue+lenght]
                            hover_colors = ['rgba(95, 140, 200, 0.7)' if y < 25 else 'rgba(208, 100, 88, 0.7)' for y in y_data]
                            
                            # Check if this is a single isolated value (surrounded by gaps)
                            align_start = residue + alignment
                            align_end = residue + alignment + lenght
                            
                            is_isolated = (lenght == 1 and 
                                         (align_start == 0 or (align_start > 0 and self.oneLetterSequence[iterator][align_start - 1] == '-')) and
                                         (align_end >= len(self.oneLetterSequence[iterator]) or self.oneLetterSequence[iterator][align_end] == '-'))
                            
                            if is_isolated:
                                # Add a vertical line from baseline (25) to the value for isolated single value
                                x_pos = np.linspace(1, len(self.residueNumList), len(self.residueNumList)).tolist()[residue+alignment]
                                marker_color = 'rgba(95, 140, 200, 0.7)' if y_data[0] < 25 else 'rgba(208, 100, 88, 0.7)'
                                
                                # Add vertical line from baseline to value (no hover on line itself)
                                figure.add_trace(go.Scatter(
                                    x=[x_pos, x_pos],
                                    y=[25, y_data[0]],
                                    mode='lines',
                                    line=dict(
                                        width=6,
                                        color=marker_color
                                    ),
                                    hoverinfo='skip',
                                    showlegend=False,
                                ), row=row + iterator, col=col)
                                
                                # Add invisible marker at the top for hover only
                                figure.add_trace(go.Scatter(
                                    x=[x_pos],
                                    y=y_data,
                                    mode='markers',
                                    marker=dict(size=10, opacity=0),  # Invisible but hoverable
                                    customdata=customdata[residue:residue+lenght],
                                    hovertemplate='%{customdata}<br>Relative Surface Accessibility: %{y:.2f} ' + self._align_pos_template() + '<extra></extra>',
                                    hoverlabel=dict(bgcolor=hover_colors[0]),
                                    showlegend=False,
                                ), row=row + iterator, col=col)
                            else:
                                # Normal line plot for connected segments
                                # Check if there's a gap before this segment (or at start)
                                prev_is_gap = (align_start > 0 and 
                                              align_start - 1 < len(self.oneLetterSequence[iterator]) and 
                                              self.oneLetterSequence[iterator][align_start - 1] == '-') or align_start == 0
                                
                                # Check if there's a gap after this segment
                                next_is_gap = (align_end < len(self.oneLetterSequence[iterator]) and 
                                              self.oneLetterSequence[iterator][align_end] == '-') or align_end >= len(self.oneLetterSequence[iterator])
                                
                                # Extend by 0.5 on each side if there are gaps
                                start_offset = -0.5 if prev_is_gap else 0.0
                                end_offset = 0.5 if next_is_gap else 0.0
                                
                                # Calculate x positions with extension for visual representation
                                start_x = residue + alignment + 1 + start_offset
                                end_x = residue + alignment + lenght + end_offset
                                x_positions = np.linspace(start_x, end_x, lenght).tolist()
                                
                                # Original x positions for hover (without extension)
                                x_positions_hover = np.linspace(1, len(self.residueNumList), len(self.residueNumList)).tolist()[residue+alignment:residue+lenght+alignment]
                                
                                # Add line with extended positions (no hover)
                                figure.add_trace(go.Scatter(
                                    x=x_positions,
                                    y=y_data,
                                    hoverinfo='skip',
                                    showlegend=False,
                                    line={'color': 'rgba(0,0,0,0)', 'shape': 'spline', 'width': 1},
                                ), row=row + iterator, col=col)

                                # Add fill gradient with extended positions
                                figure.add_trace(go.Scatter(
                                    x=x_positions,
                                    y=[25] * lenght,
                                    showlegend=False,
                                    line={'width': 0},
                                    fill='tonexty',
                                    fillgradient={
                                        'colorscale': [(0, 'rgba(95, 140, 200, 0.7)'), (treshold, 'rgba(95, 140, 200, 0.7)'), (treshold, 'rgba(208, 100, 88, 0.7)'),
                                                       (1, 'rgba(208, 100, 88, 0.7)')], 'type': 'vertical'},
                                    hoverinfo='skip',
                                    marker = {'opacity':0}
                                ), row=row+iterator, col=col)
                                
                                # Add invisible markers at original positions for hover
                                figure.add_trace(go.Scatter(
                                    x=x_positions_hover,
                                    y=y_data,
                                    mode='markers',
                                    marker=dict(size=10, opacity=0),  # Invisible but hoverable
                                    customdata=customdata[residue:residue+lenght],
                                    hovertemplate='%{customdata}<br>Relative Surface Accessibility: %{y:.2f} ' + self._align_pos_template() + '<extra></extra>',
                                    hoverlabel=dict(bgcolor=hover_colors),
                                    showlegend=False,
                                ), row=row + iterator, col=col)
                            
                            residue += lenght
                            alignment += 1

                        else:
                            residue += lenght
                            alignment += 1
                            continue

                    else:

                        alignment += 1
                iterator += 1

                maxlen = max([len(seq) for seq in self.oneLetterSequence])
                figure.update_layout(xaxis={'range': [1, maxlen + 1]})
        
        # Add horizontal white line at y=25 for all SASA plots
        num_proteins = len(self.threeLetterSequence) if self.sequenceGaps else 1
        for i in range(num_proteins):
            figure.add_trace(go.Scatter(
                x=[1, max([len(seq) for seq in self.oneLetterSequence]) + 1],
                y=[25, 25],
                mode='lines',
                line=dict(color='white', width=2),
                showlegend=False,
                hoverinfo='skip'
            ), row=(row + i if row else 1), col=col)


    def makeBindingSurfacePlots(self, figure:plotly.graph_objs.Figure, row:int=None, col:int=None) -> None:
        iterator = 0
        for seq in self.threeLetterSequence:
            corrected_binding_surface_yn = ['No'] * len(self.oneLetterSequence[iterator])

            for group_name, start_pos, end_pos in self.bindingSurfaceSegments(iterator):
                residue_counter = 0
                for alignment_pos in range(len(self.oneLetterSequence[iterator])):
                    if self.oneLetterSequence[iterator][alignment_pos] != '-':
                        residue_counter += 1
                        if start_pos <= residue_counter <= end_pos:
                            corrected_binding_surface_yn[alignment_pos] = f'Yes [{group_name}]'

            customdata = []
            hover_colors = []
            residue_counter = 0
            authResNums = self.alignedResidueNumbers(iterator)  # original (auth) residue numbers for the hover

            for alignment_pos in range(len(self.oneLetterSequence[iterator])):
                if self.oneLetterSequence[iterator][alignment_pos] != '-':
                    residue_counter += 1
                    auth_num = authResNums[alignment_pos]
                    residue_letter = self.threeLetterSequence[iterator][alignment_pos]
                    binding_status = corrected_binding_surface_yn[alignment_pos]
                    # Extract group name and set hover color
                    if binding_status.startswith('Yes [') and binding_status.endswith(']'):
                        group_name = binding_status[5:-1]
                        customdata.append(f'{residue_letter}{auth_num}<br>Binding surface: {group_name}')
                        if group_name in self.MOLSTAR_COLOR_BY_GROUP:
                            rgb = self.MOLSTAR_COLOR_BY_GROUP[group_name]
                            hover_colors.append(f'rgba({rgb[0]},{rgb[1]},{rgb[2]},0.85)')
                        else:
                            hover_colors.append('rgba(128,128,128,0.7)')
                    else:
                        customdata.append(f'{residue_letter}{auth_num}')
                        hover_colors.append('rgba(176,190,197,0.7)')
                else:
                    customdata.append('')
                    hover_colors.append('rgba(176,190,197,0.7)')

            figure.add_trace(go.Scatter(
                x= [i+1 if self.oneLetterSequence[iterator][i] != '-' else None for i in range(len(self.oneLetterSequence[iterator]))], y=np.linspace(0, 0, len(self.residueNumList)),
                customdata=customdata,
                # text=alignedSeq,
                hovertemplate= '%{customdata}' + self._align_pos_template() + '<extra></extra>',
                hoverlabel=dict(bgcolor=hover_colors),
                showlegend=False,
                mode='lines',
                line=dict(color='#B0BEC5', width=2, dash='dot')  # Grey dashed line like undefined in sec struct
            ), row=row+iterator, col=col)
            
            # Add markers for isolated non-binding surface residues (surrounded by gaps)
            isolated_x = []
            isolated_y = []
            isolated_customdata = []
            isolated_hover_colors = []
            
            for i in range(len(self.oneLetterSequence[iterator])):
                if self.oneLetterSequence[iterator][i] != '-':
                    # Check if this residue is isolated (surrounded by gaps)
                    prev_is_gap = (i == 0) or (i > 0 and self.oneLetterSequence[iterator][i-1] == '-')
                    next_is_gap = (i == len(self.oneLetterSequence[iterator]) - 1) or (i < len(self.oneLetterSequence[iterator]) - 1 and self.oneLetterSequence[iterator][i+1] == '-')
                    
                    # Check if it's NOT part of binding surface
                    is_not_binding = corrected_binding_surface_yn[i] == 'No'
                    
                    if prev_is_gap and next_is_gap and is_not_binding:
                        isolated_x.append(i + 1)
                        isolated_y.append(0)
                        isolated_customdata.append(customdata[i])
                        isolated_hover_colors.append(hover_colors[i])
            
            # Add isolated non-binding surface as dashed line segments
            for i in range(len(isolated_x)):
                x_pos = isolated_x[i]
                # Draw a small dashed line segment from -0.5 to +0.5 around the position
                figure.add_trace(go.Scatter(
                    x=[x_pos - 0.5, x_pos + 0.5],
                    y=[0, 0],
                    mode='lines',
                    line=dict(color='#B0BEC5', width=2, dash='dot'),
                    hoverinfo='skip',
                    showlegend=False,
                ), row=row+iterator, col=col)
                
                # Add invisible marker for hover
                figure.add_trace(go.Scatter(
                    x=[x_pos],
                    y=[0],
                    mode='markers',
                    marker=dict(size=10, opacity=0),
                    customdata=[isolated_customdata[i]],
                    hovertemplate='%{customdata}' + self._align_pos_template() + '<extra></extra>',
                    hoverlabel=dict(bgcolor=[isolated_hover_colors[i]]),
                    showlegend=False,
                ), row=row+iterator, col=col)

            for group_name, start_pos, end_pos in self.bindingSurfaceSegments(iterator):
                binding_x_positions = []
                binding_y_positions = []
                residue_counter = 0
                
                for alignment_pos in range(len(self.oneLetterSequence[iterator])):
                    if self.oneLetterSequence[iterator][alignment_pos] != '-':
                        residue_counter += 1
                        if start_pos <= residue_counter <= end_pos:
                            binding_x_positions.append(alignment_pos + 1)
                            binding_y_positions.append(0)
                    else:
                        if binding_x_positions:
                            temp_counter = residue_counter
                            has_more_binding = False
                            for check_pos in range(alignment_pos + 1, len(self.oneLetterSequence[iterator])):
                                if self.oneLetterSequence[iterator][check_pos] != '-':
                                    temp_counter += 1
                                    if temp_counter <= end_pos:
                                        has_more_binding = True
                                        break
                                    else:
                                        break
                            
                            if has_more_binding:
                                binding_x_positions.append(None)
                                binding_y_positions.append(None)
                
                color = self.MOLSTAR_COLOR_BY_GROUP.get(group_name, (128, 128, 128))
                
                if binding_x_positions:
                    figure.add_trace(go.Scatter(
                        x=binding_x_positions,
                        y=binding_y_positions,
                        line=dict(width=10, color=f'rgb({color[0]}, {color[1]}, {color[2]})'),
                        hoverinfo='skip',
                        showlegend=False,
                        marker={'opacity': 0}
                    ), row=row+iterator, col=col)
                    
                    marker_x = [x for x in binding_x_positions if x is not None]
                    marker_y = [y for y in binding_y_positions if y is not None]
                    
                    if marker_x:
                        figure.add_trace(go.Scatter(
                            x=marker_x,
                            y=marker_y,
                            mode='markers',
                            marker=dict(
                                size=10,
                                color=f'rgb({color[0]}, {color[1]}, {color[2]})',
                                symbol='square',
                                # line=dict(width=2, color='black')
                            ),
                            hoverinfo='skip',
                            showlegend=False
                        ), row=row+iterator, col=col)
            
            iterator += 1

            # this trace allows for charts to align

        figure.add_trace(go.Scatter(x = np.linspace(1, len(self.residueNumList), len(self.residueNumList)),
                                    y= np.linspace(1, len(self.residueNumList), len(self.residueNumList)),
                                    hoverinfo = 'skip',
                                    line_width = 0,
                                    marker_color = 'white'

                         ))
        
        # Set x-axis range to start from 1 for binding surface plots
        maxlen = max([len(seq) for seq in self.oneLetterSequence])
        figure.update_layout(xaxis={'range': [1, maxlen + 1]})
        

    def _addRegionBar(self, figure, row, col, iterator, ranges, y0, y1, region_key, color_of,
                      title, alpha):
        """A bar of named regions above the residue row (between y0 and y1): one faded band per
        region piece with its name. Its invisible hover/click targets are tagged meta.loopBar, so
        the frontend highlights the whole region in 3D on hover and toggles it as sticks on click,
        using molstarColoring[region_key]. ranges: [(name, [start, end])] sequence positions."""
        seq = self.oneLetterSequence[iterator]
        authResNums = self.alignedResidueNumbers(iterator)
        mid_y = (y0 + y1) / 2
        # The frontend finds a band by its index among the plot's filled traces.
        fill_index = sum(1 for t in figure.data if getattr(t, 'fill', None) == 'toself')
        for group, rng in ranges:
            try:
                first, last = int(rng[0]), int(rng[1])
            except (ValueError, TypeError, IndexError):
                continue
            columns = []
            position = 0
            for column, residue in enumerate(seq, start=1):
                if residue != '-':
                    position += 1
                    if first <= position <= last:
                        columns.append(column)
            if not columns:
                continue
            x_start, x_end = columns[0], columns[-1]
            rgb = color_of(group)
            band_color = f'rgba({rgb[0]},{rgb[1]},{rgb[2]},{alpha})'
            x0, x1 = x_start - 0.5, x_end + 0.5
            figure.add_trace(go.Scatter(          # visible band (faded, no frame)
                x=[x0, x1, x1, x0, x0], y=[y0, y0, y1, y1, y0],
                mode='lines', fill='toself', fillcolor=band_color, line=dict(width=0),
                hoverinfo='skip', name='', showlegend=False,
            ), row=row, col=col)
            figure.add_trace(go.Scatter(          # region name
                x=[(x_start + x_end) / 2], y=[mid_y], mode='text', text=[group],
                textfont=dict(size=9, color='#20303a'), hoverinfo='skip', showlegend=False,
            ), row=row, col=col)
            hover_x = [float(c) for c in range(x_start, x_end + 1)]
            tooltip = (f"<b>{title} {group}</b><br>Residues {authResNums[x_start - 1]}–"
                       f"{authResNums[x_end - 1]}<br>click to show as sticks")
            figure.add_trace(go.Scatter(          # invisible hover/click targets
                x=hover_x, y=[mid_y] * len(hover_x), mode='markers',
                marker=dict(size=14, opacity=0),
                customdata=[tooltip] * len(hover_x),
                hovertemplate=' %{customdata}<extra></extra>',
                hoverlabel=dict(bgcolor=f'rgb({rgb[0]},{rgb[1]},{rgb[2]})',
                                font=dict(color='white', size=HOVER_FONT_SIZE)),
                meta={'loopBar': True, 'group': group, 'regionKey': region_key,
                      'segIndex': fill_index, 'baseColor': band_color},
                showlegend=False,
            ), row=row, col=col)
            fill_index += 1

    def disulfideBonds(self, iterator: int) -> list:
        """Disulfide bonds of protein `iterator` whose two Cys are in the shown columns, as
        dicts: 'a', 'b' (Cys residue numbers), 'distance' (S-S, Å), 'group' ('Cys56–Cys178'),
        'columns' (1-based alignment columns) and 'positions' (sequence positions). Bonds of the
        part past the catalytic domain drop out when that part is masked."""
        where = {}
        position = 0
        for column, (residue, number) in enumerate(
                zip(self.oneLetterSequence[iterator], self.alignedResidueNumbers(iterator)), start=1):
            if residue == '-':
                continue
            position += 1
            if residue == 'C' and number is not None:
                where.setdefault(number, (column, position))
        bonds = []
        for bond in self.disulfideList[iterator] if iterator < len(self.disulfideList) else []:
            try:
                a, b, distance = int(bond[0]), int(bond[1]), float(bond[2])
            except (TypeError, ValueError, IndexError):
                continue
            if a in where and b in where:
                bonds.append({'a': a, 'b': b, 'distance': distance, 'group': f'Cys{a}–Cys{b}',
                              'columns': (where[a][0], where[b][0]),
                              'positions': (where[a][1], where[b][1])})
        return bonds

    def _regionAt(self, iterator: int, position: int):
        """The core strand or loop ('β9', 'L2') holding a sequence position, else None."""
        try:
            for name, pieces in self.betaCoreDict[iterator].items():
                if any(int(start) <= position <= int(end) for start, end in pieces):
                    return name
            for name, (start, end) in self.loopRegionsDict[iterator].items():
                if int(start) <= position <= int(end):
                    return name
        except (ValueError, TypeError, IndexError):
            pass
        return None

    def bindingSurfaceSegments(self, iterator: int) -> list:
        """[(region, start, end)] of the binding surface (pipeline step 7c), in sequence
        positions, ordered by start. A region holds several runs, so the list is longer than
        the number of regions; 'other' holds the residues outside every loop and strand."""
        try:
            return sorted(((region, start, end)
                           for region, segments in self.bindingSurfaceDict[iterator].items()
                           for start, end in segments), key=lambda item: item[1])
        except (IndexError, TypeError, ValueError):
            return []

    def _addBindingSurfaceLane(self, figure, row, col, iterator):
        """The binding-surface lane right below the secondary structure row: a tick under every
        residue the substrate plane touches (pipeline step 7c), one trace per region. The
        traces are tagged meta.loopBar (dots=True) with regionKey 'bindingSurface', so hovering
        a tick highlights all the protein's binding-surface ticks and residues, in the plot and
        in 3D, and a click shows them all as sticks."""
        seq = self.oneLetterSequence[iterator]
        three_letter = self.threeLetterSequence[iterator]
        authResNums = self.alignedResidueNumbers(iterator)
        column_of, position = {}, 0
        for column, residue in enumerate(seq, start=1):
            if residue != '-':
                position += 1
                column_of[position] = column

        by_region = {}
        for region, start, end in self.bindingSurfaceSegments(iterator):
            for position in range(start, end + 1):
                if position in column_of:
                    by_region.setdefault(region, []).append(position)
        rgb = BINDING_SURFACE_COLOR
        base_color = f'rgb({rgb[0]},{rgb[1]},{rgb[2]})'
        # the aromatic binding residues (Phe, Trp, Tyr, His; not the brace histidines) get a longer
        # tick (overview.binding_aromatics)
        brace = {int(n) for n in (self.histidineBraceTuple[iterator] or []) if n is not None} \
            if iterator < len(self.histidineBraceTuple) else set()

        def aromatic(column):
            aa = seq[column - 1].upper()
            return aa in 'FWYH' and not (aa == 'H' and authResNums[column - 1] in brace)

        for region, positions in by_region.items():
            columns = [column_of[position] for position in positions]
            tooltips = []
            # 'other' holds the residues outside every loop and core strand -- the N-His
            outside = region == 'other'
            for column in columns:
                name = three_letter[column - 1] if column - 1 < len(three_letter) else ''
                number = authResNums[column - 1]
                tooltips.append(
                    f"<b>Binding surface{' · aromatic' if aromatic(column) else ''}</b><br>{name.capitalize()}{number}"
                    + ("<br>outside the loops and core strands" if outside else f" in {region}")
                    + "<br>click to show all binding-surface residues of this protein as sticks")
            figure.add_trace(go.Scatter(
                x=columns, y=[BINDING_TICK_Y] * len(columns), mode='markers',
                marker=dict(symbol='line-ns', size=[15 if aromatic(c) else 9 for c in columns],
                            color=base_color, line=dict(color=base_color, width=3)),
                customdata=tooltips,
                hovertemplate=' %{customdata}<extra></extra>',
                hoverlabel=dict(bgcolor=base_color,
                                font=dict(color='white', size=HOVER_FONT_SIZE)),
                meta={'loopBar': True, 'dots': True, 'group': region,
                      'regionKey': 'bindingSurface', 'baseColor': base_color},
                showlegend=False,
            ), row=row, col=col)

    def _addDisulfideLane(self, figure, row, col, iterator):
        """The disulfide lane below the secondary structure row: a dot on each Cys of a
        disulfide bond, one trace per bond. The traces are tagged meta.loopBar (dots=True) with
        regionKey 'disulfides': hovering a dot names both Cys and highlights all the protein's
        bonded Cys, in the plot and in 3D; a click shows them all as sticks."""
        rgb = DISULFIDE_COLOR
        base_color = f'rgb({rgb[0]},{rgb[1]},{rgb[2]})'
        for bond in self.disulfideBonds(iterator):
            cys = [(bond['a'], self._regionAt(iterator, bond['positions'][0])),
                   (bond['b'], self._regionAt(iterator, bond['positions'][1]))]
            labels = [f'Cys{number}' + (f' ({region})' if region else '') for number, region in cys]
            tooltips = [f"<b>Disulfide bond {bond['group']}</b><br>{labels[k]} bonded to "
                        f"{labels[1 - k]}<br>S–S {bond['distance']:.2f} Å<br>"
                        f"click to show all disulfide bonds of this protein as sticks" for k in (0, 1)]
            figure.add_trace(go.Scatter(
                x=list(bond['columns']), y=[DISULFIDE_DOT_Y] * 2, mode='markers',
                marker=dict(size=10, color=base_color, line=dict(color='white', width=2)),
                customdata=tooltips,
                hovertemplate=' %{customdata}<extra></extra>',
                hoverlabel=dict(bgcolor='rgb(138,100,8)',
                                font=dict(color='white', size=HOVER_FONT_SIZE)),
                meta={'loopBar': True, 'dots': True, 'group': bond['group'],
                      'regionKey': 'disulfides', 'baseColor': base_color},
                showlegend=False,
            ), row=row, col=col)

    def makeSecStructLoopsPlots(self, figure: plotly.graph_objs.Figure, row: int = None, col: int = 1) -> None:
        """Secondary structure row of each protein with two bars above it -- the core beta
        strands (β1, β2, ...: the same number is the same strand in every LPMO) and the loop
        regions -- and two lanes below it: the binding-surface residues and the disulfide
        bonds."""
        for iterator, struct in enumerate(self.StructList):
            sec_struct_plot(figure, struct, row + iterator, protein_id=self.protId[iterator],
                            three_letter_sequence=self.threeLetterSequence[iterator],
                            auth_numbers=self.alignedResidueNumbers(iterator),
                            glyph_scale=SS_GLYPH_SCALE)
            core_pieces = [(name, piece) for name, pieces in self.betaCoreDict[iterator].items()
                           for piece in pieces]
            self._addRegionBar(figure, row + iterator, col, iterator, core_pieces,
                               CORE_BAR_Y0, CORE_BAR_Y1, 'betaCore',
                               lambda group: CORE_STRAND_COLOR, 'Core strand', 0.40)
            self._addRegionBar(figure, row + iterator, col, iterator,
                               list(self.loopRegionsDict[iterator].items()),
                               LOOP_BAR_Y0, LOOP_BAR_Y1, 'loopRegions',
                               lambda group: self.MOLSTAR_COLOR_BY_GROUP.get(group, (128, 128, 128)),
                               'Loop', 0.20)
            self._addBindingSurfaceLane(figure, row + iterator, col, iterator)
            self._addDisulfideLane(figure, row + iterator, col, iterator)
            figure.update_yaxes(range=[SS_LANE_YMIN, SS_BARS_YMAX], row=row + iterator, col=col)

    def makeHoleHoppingPlots(self, figure: plotly.graph_objs.Figure, row: int = None, col: int = None) -> None:
        '''Plots hole hopping residue roles as a heatmap: donors and bridging residues.
        Each protein uses its own colour pair (see hole_hopping_colors), keyed by its position
        in the comparison; single-protein subsets carry that position in holeHoppingColorPosition.'''

        colorPosition = getattr(self, 'holeHoppingColorPosition', 0)
        iterator = 0
        annotationList = []

        for protein_hh in self.holeHoppingList:
            colors = hole_hopping_colors(colorPosition + iterator)
            DONOR_COLOR, BRIDGE_COLOR = colors['donor'], colors['bridge']
            z_values = []
            customdata = []
            hover_colors = []
            alignedResNums = self.alignedResidueNumbers(iterator)

            for idx, letter in enumerate(self.threeLetterSequence[iterator]):
                role = protein_hh[idx]
                if letter != '---':
                    residue_number = alignedResNums[idx]
                    if role == 'donor':
                        z_values.append(2)
                        customdata.append(f'{letter}{residue_number}<br>Hole hopping: <b>Donor</b> (surface-exposed)')
                        hover_colors.append(f'rgba({DONOR_COLOR[0]},{DONOR_COLOR[1]},{DONOR_COLOR[2]},0.9)')
                    elif role == 'bridge':
                        z_values.append(1)
                        customdata.append(f'{letter}{residue_number}<br>Hole hopping: <b>Bridge</b>')
                        hover_colors.append(f'rgba({BRIDGE_COLOR[0]},{BRIDGE_COLOR[1]},{BRIDGE_COLOR[2]},0.9)')
                    else:
                        z_values.append(None)
                        customdata.append(f'{letter}{residue_number}')
                        hover_colors.append('rgba(176,190,197,0.7)')
                else:
                    z_values.append(None)
                    customdata.append('')
                    hover_colors.append('rgba(176,190,197,0.7)')

            colorscale = [
                [0.0, f'rgb({BRIDGE_COLOR[0]},{BRIDGE_COLOR[1]},{BRIDGE_COLOR[2]})'],
                [0.5, f'rgb({BRIDGE_COLOR[0]},{BRIDGE_COLOR[1]},{BRIDGE_COLOR[2]})'],
                # [0.5, f'rgb({DONOR_COLOR[0]},{DONOR_COLOR[1]},{DONOR_COLOR[2]})'],
                [1.0, f'rgb({DONOR_COLOR[0]},{DONOR_COLOR[1]},{DONOR_COLOR[2]})'],
            ]

            figure.add_trace(go.Heatmap(
                x=[int(x) for x in self.residueNumList],
                y=np.linspace(0, 0, len(self.threeLetterSequence[iterator]), dtype='int').tolist(),
                z=z_values,
                showscale=False,
                colorscale=colorscale,
                zmin=0,
                zmax=2,
                opacity=1,   # soft palette: full opacity so the row matches the 3D colours
                hoverinfo='skip',
                name='hole_hopping',
                showlegend=False,
            ), col=col, row=row + iterator)

            # Add invisible scatter for hover
            x_no_gaps = []
            y_no_gaps = []
            customdata_no_gaps = []
            hover_colors_no_gaps = []
            for idx, cd in enumerate(customdata):
                if cd != '':
                    x_no_gaps.append(int(self.residueNumList[idx]))
                    y_no_gaps.append(0)
                    customdata_no_gaps.append(cd)
                    hover_colors_no_gaps.append(hover_colors[idx])

            figure.add_trace(go.Scatter(
                x=x_no_gaps,
                y=y_no_gaps,
                mode='markers',
                marker=dict(size=10, opacity=0),
                customdata=customdata_no_gaps,
                hovertemplate=' %{customdata}' + self._align_pos_template() + '<extra></extra>',
                hoverlabel=dict(bgcolor=hover_colors_no_gaps),
                showlegend=False,
                hoverinfo='text'
            ), col=col, row=row + iterator)

            for num in self.residueNumList:
                y_axis_idx = (row + iterator) if (row is not None) else (iterator + 1)
                y_axis_name = 'y' if y_axis_idx == 1 else f'y{y_axis_idx}'
                yref_str = f'{y_axis_name} domain'

                annotationList.append({'font': {'color': 'black', 'size': 11, 'weight': 'bold'},
                                       'showarrow': False,
                                       'text': self.oneLetterSequence[iterator][int(num) - 1],
                                       'x': int(num),
                                       'xanchor': 'center',
                                       'xref': 'x',
                                       'y': 0.5,
                                       'yanchor': 'middle',
                                       'yref': yref_str
                                       })

            iterator += 1

        figure.update_layout(annotations=list(figure['layout']['annotations']) + annotationList)

    def _create_subset_instance(self, index: int):
        #instance with data for single protein 

        subset_plotter = makePlots([], None)
        subset_plotter.protId = [self.protId[index]]
        subset_plotter.protName = [self.protName[index]]
        subset_plotter.oneLetterSequence = [self.oneLetterSequence[index]]
        subset_plotter.threeLetterSequence = [self.threeLetterSequence[index]]
        subset_plotter.flexibilityList = [self.flexibilityList[index]]
        subset_plotter.relSASAList = [self.relSASAList[index]]
        subset_plotter.hydrophobicityList = [self.hydrophobicityList[index]]
        subset_plotter.allHydrophobicityLists = {name: [lists[index]] for name, lists in self.allHydrophobicityLists.items()}
        subset_plotter.allHydrophobicityScales = self.allHydrophobicityScales
        subset_plotter.hydrophobicityScaleNames = self.hydrophobicityScaleNames
        subset_plotter.electrostaticsList = [self.electrostaticsList[index]] if hasattr(self, 'electrostaticsList') else [[]]
        subset_plotter.atomChargesList = [self.atomChargesList[index]] if hasattr(self, 'atomChargesList') and index < len(self.atomChargesList) else [[]]
        subset_plotter.maxElectrostaticsVal = self.maxElectrostaticsVal if hasattr(self, 'maxElectrostaticsVal') else None
        subset_plotter.minElectrostaticsVal = self.minElectrostaticsVal if hasattr(self, 'minElectrostaticsVal') else None
        subset_plotter.StructList = [self.StructList[index]]
        subset_plotter.bindingSurfaceYN = [self.bindingSurfaceYN[index]]
        subset_plotter.histidineBraceTuple = [self.histidineBraceTuple[index]]
        subset_plotter.bindingSurfaceDict = [self.bindingSurfaceDict[index]]
        subset_plotter.loopRegionsDict = [self.loopRegionsDict[index]]
        subset_plotter.betaCoreDict = [self.betaCoreDict[index]]
        subset_plotter.disulfideList = [self.disulfideList[index]]
        subset_plotter.loopRegionsYN = [self.loopRegionsYN[index]]
        subset_plotter.holeHoppingList = [self.holeHoppingList[index]] if hasattr(self, 'holeHoppingList') and index < len(self.holeHoppingList) else [[]]
        subset_plotter.holeHoppingPaths = [self.holeHoppingPaths[index]] if hasattr(self, 'holeHoppingPaths') and index < len(self.holeHoppingPaths) else [[]]
        subset_plotter.residueNumList = self.residueNumList
        subset_plotter.residueNumbers = [self.residueNumbers[index]] if index < len(self.residueNumbers) else [[]]
        subset_plotter.domainRanges = [self.domainRanges[index]] if index < len(self.domainRanges) else [[]]
        subset_plotter.totalProteinCount = getattr(self, 'totalProteinCount', 1)
        subset_plotter.sequenceGaps = [self.sequenceGaps[index]] if self.sequenceGaps else []
        subset_plotter.maxFlexVal = self.maxFlexVal
        subset_plotter.minFlexVal = self.minFlexVal
        # Per-method flexibility for this one protein (global per-method min/max stay unchanged).
        subset_plotter.allFlexibilityLists = {m: [lists[index]] for m, lists in self.allFlexibilityLists.items()
                                              if index < len(lists)}
        subset_plotter.flexMinByMethod = self.flexMinByMethod
        subset_plotter.flexMaxByMethod = self.flexMaxByMethod
        subset_plotter.flexMinByMethodCat = self.flexMinByMethodCat
        subset_plotter.flexMaxByMethodCat = self.flexMaxByMethodCat
        subset_plotter.flexibilityMethods = self.flexibilityMethods
        subset_plotter.defaultFlexMethod = self.defaultFlexMethod
        subset_plotter.aaHydrophobicity = self.aaHydrophobicity
        subset_plotter.secondSphereResidues = [self.secondSphereResidues[index]]
        subset_plotter.secondSphereMotifs = [self.secondSphereMotifs[index]]
        subset_plotter.braceMotifs = [self.braceMotifs[index]]
        subset_plotter.braceInteractionResidues = [self.braceInteractionResidues[index]]
        if hasattr(self, 'secondSphereCombinedResidues'):
            subset_plotter.secondSphereCombinedResidues = [self.secondSphereCombinedResidues[index]]
        if hasattr(self, 'axialResiduePositions'):
            subset_plotter.axialResiduePositions = [self.axialResiduePositions[index]]
        if hasattr(self, 'piStackingResiduePositions'):
            subset_plotter.piStackingResiduePositions = [self.piStackingResiduePositions[index]]
        if hasattr(self, 'cationPiResiduePositions'):
            subset_plotter.cationPiResiduePositions = [self.cationPiResiduePositions[index]]
        if hasattr(self, 'hbondResiduePositions'):
            subset_plotter.hbondResiduePositions = [self.hbondResiduePositions[index]]
        
        return subset_plotter

    def make_plot_for_protein(self, plotType: str, protein_index: int, scale_name: str=None) -> plotly.graph_objects.Figure:

        figure = make_subplots(
            rows=1, cols=1,
            specs=[[{"type": "xy"}]],
        )

        temp_plotter = self._create_subset_instance(protein_index)
        temp_plotter.holeHoppingColorPosition = protein_index   # per-protein colour by comparison order

        if plotType == 'Sequence':
            temp_plotter.makeSequencePlots(figure, row=1, col=1)
        elif plotType == 'Flexibility':
            # scale_name doubles as the flexibility method selector (cabsflex/gnm/…).
            temp_plotter.makeFlexibilityPlot(figure, row=1, col=1, showLegend=False, method_name=scale_name)
        elif plotType == 'Hydrophobicity':
            temp_plotter.makeHydrophobicityPlots(figure, row=1, col=1, scale_name=scale_name)
        elif plotType == 'Charge':
            temp_plotter.makeChargePlots(figure, row=1, col=1)
        elif plotType == SS_LOOPS_PLOT_TYPE:
            temp_plotter.makeSecStructLoopsPlots(figure, row=1, col=1)
        elif plotType == 'Relative Surface Accessibility':
            temp_plotter.makeSASAPlots(figure, row=1, col=1)
        elif plotType == 'Binding Surface':
            temp_plotter.makeBindingSurfacePlots(figure, row=1, col=1)
        elif plotType == 'Hole Hopping':
            temp_plotter.makeHoleHoppingPlots(figure, row=1, col=1)

        display_name = temp_plotter.protName[0]

        maxProtLen = max(len(seq) for seq in self.oneLetterSequence)
        PIXELS_PER_RESIDUE = 10
        plot_width = max(maxProtLen * PIXELS_PER_RESIDUE, 1200)

        figure.update_layout(
            # hoverdistance=-1 removes the distance cutoff, so the nearest point always
            # wins no matter how far the cursor is. On the Sequence plot that lets the
            # invisible residue-marker scatter swallow hovers aimed at the domain bar
            # above it, so there we use a finite radius.
            #
            # With the 63px sequence plot, residue markers sit ~43px from the top and the
            # bar ~3-14px, ~30px apart; markers are size 10 so the effective gap is ~25px.
            # The cutoff must stay below that or the residues steal the bar's hover; 15
            # leaves margin both ways. Recheck if DATA_ROW_HEIGHT or DOMAIN_BAR_Y0/Y1 change.
            # The Secondary Structure and Loops plot has two bars ~13px apart above its
            # residue row, so its radius is smaller still.
            hoverdistance=15 if plotType == 'Sequence' else 10 if plotType == SS_LOOPS_PLOT_TYPE else -1,
            hovermode='closest',
            # No drag-zoom / pan / box-select — accidental zooming destroyed the view and
            # we don't use it here. Hover and click still work.
            dragmode=False,
            hoverlabel=dict(font=dict(size=HOVER_FONT_SIZE)),
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            margin=dict(l=10, r=0, t=1, b=2),
            autosize=False,
            width=plot_width,
            showlegend=False,
            # Must match the CSS row height, otherwise the SVG is squashed by
            # `.plotWrapper > div { height: 100% !important }` and the y-axis labels drift
            # out of alignment. Sequence and Secondary Structure and Loops are taller by their
            # bars' share; every other plot is DATA_ROW_HEIGHT so their data rows line up. See
            # compare_protein.css.
            height=(SEQUENCE_PLOT_HEIGHT if plotType == 'Sequence'
                    else SS_LOOPS_PLOT_HEIGHT if plotType == SS_LOOPS_PLOT_TYPE else DATA_ROW_HEIGHT)
        )

        figure.update_xaxes(showgrid=False, zeroline=False, visible = False if protein_index != len(self.oneLetterSequence) - 1 else True,
                            dtick = 5, tick0= 5,  autotickangles = [0], title=None, range = (0.5, maxProtLen + 1))
        figure.update_yaxes(showgrid=True if plotType in ['Flexibility'] else False, 
                            gridcolor='white' if plotType in ['Flexibility'] else None,
                            gridwidth=2 if plotType in ['Flexibility'] else None,
                            zeroline=False if plotType not in ['Flexibility'] else True,
                            title=None,
                            showticklabels=False if plotType not in ['Flexibility'] else True)


        return figure


    def makeChargeColorScale(self, min_charge=None, max_charge=None, fixed_scale=True):
        """
        Create a color scale for charge values.
        Red (negative) -> White (neutral) -> Blue (positive)
        Uses fixed scale from -2 to 2, where 0 corresponds to white.
        Returns a function that maps charge value to RGB tuple.
        """
        def charge_to_rgb(charge_val):
            if charge_val is None:
                return (255, 255, 255)  # White for missing data
            
            # Use fixed scale from -2 to 2, where 0 is white
            charge_min = -2.0
            charge_max = 2.0
            charge_mid = 0.0
            
            # Clamp charge value to the scale range
            charge_val = max(charge_min, min(charge_max, charge_val))
            
            # Normalize to 0-1 range where 0.5 corresponds to charge=0
            if charge_val <= charge_mid:
                # Negative charges: map from -2 (0.0) to 0 (0.5)
                normalized = 0.5 * (charge_val - charge_min) / (charge_mid - charge_min)
            else:
                # Positive charges: map from 0 (0.5) to 2 (1.0)
                normalized = 0.5 + 0.5 * (charge_val - charge_mid) / (charge_max - charge_mid)
            
            # Color mapping:
            # 0.0 (charge = -2) -> Red (255, 0, 0)
            # 0.5 (charge = 0) -> White (255, 255, 255)
            # 1.0 (charge = 2) -> Blue (0, 0, 255)
            
            if normalized <= 0.5:
                # Red to White gradient
                ratio = normalized * 2  # 0 to 1
                r = 255
                g = int(255 * ratio)
                b = int(255 * ratio)
            else:
                # White to Blue gradient
                ratio = (normalized - 0.5) * 2  # 0 to 1
                r = int(255 * (1 - ratio))
                g = int(255 * (1 - ratio))
                b = 255
            
            return (r, g, b)
        
        return charge_to_rgb

    def makeQueryParams(self)-> list:
        # makes query params for mol*, it should be called last

        queryParams = {'sequence': [],'flexibility': [],'hydrophobicity': [],'charge': [],'bindingSurface': [], 'loopRegions': [], 'betaCore': [], 'disulfides': [], 'holeHopping': [], 'secondaryStructure': [], 'rsa': [], 'hydrophobicityScales': {}, 'flexibilityMethods': {}, 'flexibilityMethodsFull': {}}
        # Per-method flexibility colouring for the 3D viewer, so the frontend can swap methods
        # client-side (like hydrophobicity scales). Two normalisations per method: the default
        # `flexibilityMethods` scaled to the CATALYTIC-domain range (matches the default plot
        # view), and `flexibilityMethodsFull` scaled to the FULL-protein range (for the
        # "Full length" toggle). Names in display order for the dropdown.
        queryParams['flexibilityMethodNames'] = list(self.flexibilityMethods)
        for _m in self.flexibilityMethods:
            queryParams['flexibilityMethods'][_m] = []
            queryParams['flexibilityMethodsFull'][_m] = []

        #this dict is initialized in self.makeSequencePlots()
        queryParams['sequence'] = self.sequenceDict
        queryParams['metalResidueNumbers'] = self.metalResidueNumbers
        # Original residue numbers per protein — Mol* selects by auth_seq_id, so the
        # frontend must address residues by these rather than by 1-based position.
        queryParams['residueNumbers'] = self.residueNumbers
        HydrophobColorDict = self.makeHydrophobColorScale(self.aaHydrophobicity, 'aaRGB')

        # Generate molstar coloring for all hydrophobicity scales
        allHydrophobColorDicts = {}
        for scale_name in self.hydrophobicityScaleNames:
            allHydrophobColorDicts[scale_name] = self.makeHydrophobColorScale(self.allHydrophobicityScales[scale_name], 'aaRGB')
            queryParams['hydrophobicityScales'][scale_name] = []
        
        # Create charge color scale function using fixed scale -2 to 2
        chargeColorFunc = self.makeChargeColorScale()

        protLen = [len(''.join(seq)) for seq in self.sequenceGaps]
        queryParams['proteinLenght'] = protLen

        residueIndex = []

        for i, prot in enumerate(self.oneLetterSequence):
            protFlexQuery = []
            protHydrophobQuery = []
            protChargeQuery = []
            protSSQuery = []           # secondary-structure coloring
            protRSAQuery = []          # relative surface accessibility coloring
            protBS = []
            protLoops = []
            residueIndex.append([])

            # Original residue numbers from the structure file — Mol* selects by
            # auth_seq_id, so these must not be renumbered from 1.
            protResNums = self.residueNumbers[i] if i < len(self.residueNumbers) else []
            seqPos = 0
            indexNum = 0
            for residue in prot:
                if residue != '-':
                    residue_number = protResNums[seqPos] if seqPos < len(protResNums) else seqPos + 1
                    seqPos += 1
                    flexNormalized = ((self.flexibilityList[i][indexNum] - self.minFlexVal) / (self.maxFlexVal - self.minFlexVal))
                    # Floor 200 to match the plot gradient (see flexibilityMethods build below).
                    flexRGB = [255, int(200 * (1 - flexNormalized)), int(200 * (1 - flexNormalized))]

                    protFlexQuery.append({'residue_number': residue_number, 'color': {'r': flexRGB[0], 'g': flexRGB[1], 'b': flexRGB[2]}})

                    protHydrophobQuery.append({'residue_number': residue_number, 'color': {'r': HydrophobColorDict[residue][0], 'g': HydrophobColorDict[residue][1], 'b': HydrophobColorDict[residue][2]}})
                    
                    # Add charge coloring (residue-level for plots)
                    if chargeColorFunc and indexNum < len(self.electrostaticsList[i]):
                        charge_val = self.electrostaticsList[i][indexNum]
                        chargeRGB = chargeColorFunc(charge_val)
                        protChargeQuery.append({'residue_number': residue_number, 'color': {'r': chargeRGB[0], 'g': chargeRGB[1], 'b': chargeRGB[2]}})
                    else:
                        protChargeQuery.append({'residue_number': residue_number, 'color': {'r': 255, 'g': 255, 'b': 255}})
                    
                    # (atom-level charge coloring for Mol* was removed — it was never
                    # read by the frontend and its O(atoms²) build dominated makeQueryParams)

                    # Secondary structure — StructList is alignment-length, so index by
                    # the aligned position (indexNum).
                    ss_char = (self.StructList[i][indexNum]
                               if i < len(self.StructList) and indexNum < len(self.StructList[i])
                               else 'U')
                    ss_rgb = SS_MOLSTAR_COLORS.get(ss_char, (255, 255, 255))
                    protSSQuery.append({'residue_number': residue_number,
                                        'color': {'r': ss_rgb[0], 'g': ss_rgb[1], 'b': ss_rgb[2]}})

                    # Relative SASA — relSASAList is ungapped, so index by the ungapped
                    # position (seqPos was already incremented above, so use seqPos - 1).
                    rsa_val = (self.relSASAList[i][seqPos - 1]
                               if i < len(self.relSASAList) and (seqPos - 1) < len(self.relSASAList[i])
                               else None)
                    rsa_rgb = _sasa_to_rgb(rsa_val)
                    protRSAQuery.append({'residue_number': residue_number,
                                         'color': {'r': rsa_rgb[0], 'g': rsa_rgb[1], 'b': rsa_rgb[2]}})

                    residueIndex[i].append(residue_number)

                    indexNum += 1
                else:
                    residueIndex[i].append('NaN')
                    indexNum += 1

            # Binding surface and loop regions are stored as 1-based sequence
            # positions, so map them onto the structure's own residue numbering.
            def toResNum(pos):
                idx = int(pos) - 1
                if 0 <= idx < len(protResNums):
                    return protResNums[idx]
                return int(pos)

            # One entry per run of the binding surface, ordered by start residue number so the
            # JS pointer-based colouring works (a region holds several runs, see step 7c)
            for group_name, first, last in self.bindingSurfaceSegments(i):
                start, end = toResNum(first), toResNum(last)
                color = self.MOLSTAR_COLOR_BY_GROUP.get(group_name, BINDING_SURFACE_COLOR)

                protBS.append({'group': group_name, 'start_residue_number': start, 'end_residue_number': end, 'color': {'r': color[0], 'g': color[1], 'b': color[2]}})

            loop_items_sorted = sorted(
                self.loopRegionsDict[i].items(),
                key=lambda kv: int(kv[1][0])
            )

            for group_name, loop_range in loop_items_sorted:
                try:
                    start = toResNum(loop_range[0])
                    end = toResNum(loop_range[1])
                except (ValueError, TypeError, IndexError):
                    continue
                color = self.MOLSTAR_COLOR_BY_GROUP.get(group_name, (128, 128, 128))
                protLoops.append({'group': group_name, 'start_residue_number': start, 'end_residue_number': end,
                                  'color': {'r': color[0], 'g': color[1], 'b': color[2]}})

            # Core beta strands (same format as the loops, for the bar hover/click)
            protCore = [{'group': name, 'start_residue_number': toResNum(rng[0]),
                         'end_residue_number': toResNum(rng[1]),
                         'color': dict(zip('rgb', CORE_STRAND_COLOR))}
                        for name, pieces in self.betaCoreDict[i].items() for rng in pieces]

            # Disulfide bonds: both Cys of a bond under its name (for the lane hover/click)
            protDisulfides = [{'group': bond['group'], 'start_residue_number': number,
                               'end_residue_number': number,
                               'color': dict(zip('rgb', DISULFIDE_COLOR))}
                              for bond in self.disulfideBonds(i) for number in (bond['a'], bond['b'])]

            # Hole hopping coloring
            protHoleHoppingQuery = []
            hhColors = hole_hopping_colors(i)
            DONOR_RGB = dict(zip('rgb', hhColors['donor']))
            BRIDGE_RGB = dict(zip('rgb', hhColors['bridge']))
            DEFAULT_RGB = {'r': 255, 'g': 255, 'b': 255}
            if hasattr(self, 'holeHoppingList') and i < len(self.holeHoppingList):
                hh_idx = 0
                pos = 0
                for residue in prot:
                    if residue != '-':
                        res_num = protResNums[pos] if pos < len(protResNums) else pos + 1
                        pos += 1
                        role = self.holeHoppingList[i][hh_idx] if hh_idx < len(self.holeHoppingList[i]) else None
                        if role == 'donor':
                            protHoleHoppingQuery.append({'residue_number': res_num, 'color': DONOR_RGB})
                        elif role == 'bridge':
                            protHoleHoppingQuery.append({'residue_number': res_num, 'color': BRIDGE_RGB})
                        else:
                            protHoleHoppingQuery.append({'residue_number': res_num, 'color': DEFAULT_RGB})
                    hh_idx += 1

            # Build per-scale hydrophobicity coloring for molstar
            for scale_name in self.hydrophobicityScaleNames:
                scaleColorDict = allHydrophobColorDicts[scale_name]
                scaleQuery = []
                pos_s = 0
                for residue in prot:
                    if residue != '-':
                        res_num_s = protResNums[pos_s] if pos_s < len(protResNums) else pos_s + 1
                        pos_s += 1
                        if residue in scaleColorDict:
                            rgb = scaleColorDict[residue]
                            scaleQuery.append({'residue_number': res_num_s, 'color': {'r': rgb[0], 'g': rgb[1], 'b': rgb[2]}})
                        else:
                            scaleQuery.append({'residue_number': res_num_s, 'color': {'r': 255, 'g': 255, 'b': 255}})
                queryParams['hydrophobicityScales'][scale_name].append(scaleQuery)

            # Build per-method flexibility colouring for molstar, in two normalisations:
            # catalytic-domain range (default) and full-protein range ("Full length" toggle).
            # Floor 200 so the lightest colour matches the plot gradient (rgb(255,200,200)).
            for method in self.flexibilityMethods:
                fvals = self.allFlexibilityLists[method][i]
                m_meta = FLEX_METHOD_META.get(method, {})
                m_invert = m_meta.get('invert', False)
                m_alphafold = m_meta.get('scheme') == 'alphafold'

                def _flex_query(fmin, fmax):
                    span = (fmax - fmin) if (fmax is not None and fmin is not None and fmax != fmin) else None
                    out = []
                    pos = 0
                    for aln_i, residue in enumerate(prot):
                        if residue == '-':
                            continue
                        rn = protResNums[pos] if pos < len(protResNums) else pos + 1
                        pos += 1
                        v = fvals[aln_i] if aln_i < len(fvals) else None
                        if v is None or (span is None and not m_alphafold):
                            out.append({'residue_number': rn, 'color': {'r': 255, 'g': 255, 'b': 255}})
                        elif m_alphafold:
                            r, g, b = _plddt_rgb(v)
                            out.append({'residue_number': rn, 'color': {'r': r, 'g': g, 'b': b}})
                        else:
                            norm = max(0.0, min(1.0, (v - fmin) / span))   # clamp (full-range residues can exceed the catalytic max)
                            cval = (1 - norm) if m_invert else norm
                            out.append({'residue_number': rn,
                                        'color': {'r': 255, 'g': int(200 * (1 - cval)), 'b': int(200 * (1 - cval))}})
                    return out

                queryParams['flexibilityMethods'][method].append(
                    _flex_query(self.flexMinByMethodCat.get(method), self.flexMaxByMethodCat.get(method)))
                queryParams['flexibilityMethodsFull'][method].append(
                    _flex_query(self.flexMinByMethod.get(method), self.flexMaxByMethod.get(method)))

            queryParams['flexibility'].append(protFlexQuery)
            queryParams['hydrophobicity'].append(protHydrophobQuery)
            queryParams['charge'].append(protChargeQuery)
            queryParams['bindingSurface'].append(protBS)
            queryParams['loopRegions'].append(protLoops)
            queryParams['betaCore'].append(protCore)
            queryParams['disulfides'].append(protDisulfides)
            queryParams['holeHopping'].append(protHoleHoppingQuery)
            queryParams['secondaryStructure'].append(protSSQuery)
            queryParams['rsa'].append(protRSAQuery)
            queryParams['residueIndex'] = residueIndex

        queryParams['holeHoppingPaths'] = self.holeHoppingPaths
        # Brace His residue numbers per protein: the path overlay and PyMOL export place a
        # pseudo-copper at the centre of their nitrogens when a structure has no metal.
        queryParams['histidineBrace'] = [list(brace) if brace else [] for brace in self.histidineBraceTuple]
        # Per-protein hole-hopping colour pairs for the path overlay and PyMOL export.
        queryParams['holeHoppingColors'] = [
            {role: dict(zip('rgb', rgb)) for role, rgb in hole_hopping_colors(i).items()}
            for i in range(len(queryParams['holeHopping']))]
        # Default 3D flexibility colouring = the default method at the CATALYTIC range (matches
        # the default catalytic plot view); the frontend swaps to *Full on the "Full length" toggle.
        if self.defaultFlexMethod in queryParams['flexibilityMethods']:
            queryParams['flexibility'] = list(queryParams['flexibilityMethods'][self.defaultFlexMethod])
        return queryParams

    def makeExportData(self):
        """Raw data for the client-side CSV/FASTA export of each plot section, for both views:
        {'catalytic': {...}, 'full': {...}}; the catalytic one leaves out every residue past
        the protein's catalytic-domain end. residue_number is the structure's own numbering,
        position the sequence position."""
        return {'catalytic': self._exportData(catalytic=True),
                'full': self._exportData(catalytic=False)}

    def _exportData(self, catalytic):
        import re
        strip = lambda code: re.sub(r'^\d+_', '', code)
        export = {}
        n = len(self.oneLetterSequence)

        # Per protein: (alignment index, sequence index, position, residue number, residue)
        # of every residue in the view
        residues = []
        for i in range(n):
            cat_end = self.catalyticEnds[i] if catalytic else None
            numbers = self.alignedResidueNumbers(i)
            rows, seq_idx = [], 0
            for aln_idx, aa in enumerate(self.oneLetterSequence[i]):
                if aa == '-':
                    continue
                if cat_end is None or numbers[aln_idx] is None or numbers[aln_idx] <= cat_end:
                    rows.append((aln_idx, seq_idx, seq_idx + 1, numbers[aln_idx], aa))
                seq_idx += 1
            residues.append(rows)

        # Sequence — FASTA alignment, without the columns no protein has a residue in
        kept = [{r[0] for r in rows} for rows in residues]
        columns = sorted(set().union(*kept)) if kept else []
        fasta_lines = []
        for i in range(n):
            fasta_lines.append('>' + strip(self.protId[i]))
            fasta_lines.append(''.join(self.oneLetterSequence[i][c] if c in kept[i] else '-'
                                       for c in columns))
        export['Sequence'] = {'format': 'fasta', 'content': '\n'.join(fasta_lines)}

        fmt = lambda v: f'{v:.4f}' if isinstance(v, (int, float)) else ''

        def per_residue_csv(column, value):
            """One row per residue; value(i, aln_idx, seq_idx) gives the column's value."""
            rows = [f'protein,residue_number,position,residue,{column}']
            for i in range(n):
                name = strip(self.protId[i])
                for aln_idx, seq_idx, position, number, aa in residues[i]:
                    rows.append(f'{name},{number},{position},{aa},{fmt(value(i, aln_idx, seq_idx))}')
            return {'format': 'csv', 'content': '\n'.join(rows)}

        def at(lst, i, k):
            return lst[i][k] if i < len(lst) and k < len(lst[i]) else None

        # Flexibility — one column per available method (empty where a protein lacks that
        # method). allFlexibilityLists[method][i] is aligned to the sequence (None at gaps),
        # so index by alignment position.
        methods = list(self.flexibilityMethods)
        flex_rows = ['protein,residue_number,position,residue' + ''.join(',' + m for m in methods)]
        for i in range(n):
            name = strip(self.protId[i])
            for aln_idx, seq_idx, position, number, aa in residues[i]:
                cells = [fmt(at(self.allFlexibilityLists.get(m, []), i, aln_idx)) for m in methods]
                flex_rows.append(f'{name},{number},{position},{aa}' + ''.join(',' + c for c in cells))
        export['Flexibility'] = {'format': 'csv', 'content': '\n'.join(flex_rows)}

        # Hydrophobicity and charge are indexed per alignment position, relative SASA per residue
        export['Hydrophobicity'] = per_residue_csv(
            'hydrophobicity', lambda i, a, s: at(self.hydrophobicityList, i, a))
        export['Charge'] = per_residue_csv(
            'charge', lambda i, a, s: at(self.electrostaticsList, i, a))
        export['Relative Surface Accessibility'] = per_residue_csv(
            'relative_sasa', lambda i, a, s: at(self.relSASAList, i, s))

        # Secondary structure, core strands, loop regions and disulfide bonds; the regions are
        # stored by sequence position, residue_number is the structure's own numbering
        ss_map = {'H': 'helix', 'E': 'strand', 'U': 'coil'}
        rows = ['protein,residue_number,position,residue,secondary_structure,core_strand,'
                'loop_region,on_binding_surface,disulfide_partner']
        for i in range(n):
            name = strip(self.protId[i])

            def lookup(regions):
                table = {}
                for group, rng in regions:
                    if len(rng) >= 2:
                        for r in range(int(rng[0]), int(rng[1]) + 1):
                            table[r] = group
                return table
            core_lookup = lookup((name, piece) for name, pieces in self.betaCoreDict[i].items()
                                 for piece in pieces)
            loop_lookup = lookup(self.loopRegionsDict[i].items())
            bs_lookup = lookup((region, (first, last))
                               for region, first, last in self.bindingSurfaceSegments(i))
            partner = {}
            for bond in self.disulfideBonds(i):
                partner[bond['a']], partner[bond['b']] = f"Cys{bond['b']}", f"Cys{bond['a']}"
            for idx, seq_idx, position, number, aa in residues[i]:
                ss_char = self.StructList[i][idx] if idx < len(self.StructList[i]) else ''
                rows.append(f'{name},{number},{position},{aa},{ss_map.get(ss_char, ss_char)},'
                            f'{core_lookup.get(position, "")},{loop_lookup.get(position, "")},'
                            f'{"yes" if position in bs_lookup else "no"},{partner.get(number, "")}')
        export[SS_LOOPS_PLOT_TYPE] = {'format': 'csv', 'content': '\n'.join(rows)}

        return export

