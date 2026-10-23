"""
Step 7 of the LPMOcat pipeline: the beta-sandwich core shared by the LPMOs (7a), found in the
DSSP strands (step 5) of the aligned structures (step 4) within their catalytic domains
(step 6), and the loop regions between its strands (7b).

A core strand is a block of alignment columns where strand residues (DSSP E) are common in
the LPMO families: in each column, the fraction of a family's structures with a strand
residue there, averaged over the families (so that large families do not decide alone), is
at least CORE_MIN_FRACTION; a block bridges up to BLOCK_GAP columns below it. Two consecutive
blocks are one strand when most structures with strand residues in both have them in one
continuous strand. In a structure, each strand (a run of E) takes the block it shares the
most columns with.

One more core strand is not common enough for a block: the hairpin partner of the first
strand (AA9), a strand starting at most HAIRPIN_MAX_GAP residues after it with an antiparallel
bridge to it.

The assignment is then checked against the pairing of the core strands. Its ladders are the
pairs of core strands bridged (in one direction, parallel or antiparallel) in at least
LADDER_MIN_FRACTION of the structures having both. A core strand with none of its ladders
(while a partner is present) is replaced, and a missing core strand found, by the unnamed
strand between the neighbouring core strands with the most bridges of its ladders (5MSZ: b6
is 129-131, paired with b5, not the 120-121 strand in the b6 columns, which pairs with b7;
AA13: b6 is 161-162). Strands after a core strand (up to the next) that form one of its
ladders are further pieces of it, so the next loop starts after them (AaAA15: b6 120-122 and
130-132, LS after the second). Core strands are runs of DSSP E only: isolated bridges (B)
are never core (user decision 2026-09-16: 5MSZ b6 is 129-131 alone, not also the single
bridge 118). The strands are numbered in sequence order, b1 (first block), b2 (the hairpin),
b3... (the other blocks), so a number means the same strand in every protein and a protein
without a strand skips its number. Other strands -- short strands of family-specific insertions
such as the one inside the AA10 L2 -- are not part of the core and stay unnamed.

The loops are the segments between given core strands (LOOPS), which assumes the numbering of
the CORE_STRANDS strands found in the 2026-09 data.
"""

from collections import Counter, defaultdict
from itertools import combinations

# Family-balanced fractions in the 2026-09 data: 0.79-1.00 for the strands of all families,
# 0.54 for the end of b1 (past a 3-column alignment gap; merged into b1), 0.52 for the start of
# the long b8 of AA9, AA11, AA13 and AA15 (merged into b8), then 0.36 and below (the AA9
# hairpin b2 and the strand inside the AA10 L2: 0.12).
CORE_MIN_FRACTION = 0.45
BLOCK_GAP = 2         # alignment columns
HAIRPIN_MAX_GAP = 4   # residues between b1 and its hairpin partner b2 (2-3 in AA9)
LADDER_MIN_FRACTION = 0.5
PIECE_GAP = 3         # residues
CORE_STRANDS = 9

# Loop: (strand before it, strand after it; None = the catalytic domain end). A missing strand
# is replaced by the nearest core strand present before (for the first) or after it (for the
# second); with none before, the loop starts at the first residue. LC, the C-terminal
# segment, needs LC_MIN_LENGTH residues: AA9, AA11 and AA14 have 12-92, the other families
# end 1-6 residues after b9.
LOOPS = {"L2": ("b2", "b3"), "L3": ("b4", "b5"), "LS": ("b6", "b7"), "L8": ("b8", "b9"),
         "LC": ("b9", None)}
LC_MIN_LENGTH = 8


def _segments(ss8, n_residues):
    """[(first, last)] 0-based indices of the E runs among the first n_residues residues."""
    runs, start = [], None
    for i in range(n_residues):
        if ss8[i] == "E" and start is None:
            start = i
        elif ss8[i] != "E" and start is not None:
            runs.append((start, i - 1))
            start = None
    if start is not None:
        runs.append((start, n_residues - 1))
    return runs


def _residue_columns(seq):
    """0-based alignment column of every residue of an aligned sequence."""
    return [col for col, c in enumerate(seq) if c != "-"]


def find_blocks(aligned, ss8, families, cat_end_positions):
    """Core strand blocks as [(first column, last column)] (0-based), in sequence order.

    aligned: {struct_id: aligned sequence}; ss8: {struct_id: DSSP states};
    families: {struct_id: family}; cat_end_positions: {struct_id: number of residues up to
    the catalytic domain end}."""
    segment_of = {}                 # struct_id -> {column: segment index}
    counts = defaultdict(Counter)   # family -> column -> structures with E
    members = Counter()
    for struct_id, seq in aligned.items():
        if struct_id not in ss8:
            continue
        family = families.get(struct_id)
        members[family] += 1
        cols = _residue_columns(seq)
        n = cat_end_positions.get(struct_id, len(cols))
        segment_of[struct_id] = {cols[i]: k
                                 for k, (first, last) in enumerate(_segments(ss8[struct_id], n))
                                 for i in range(first, last + 1)}
        counts[family].update(segment_of[struct_id].keys())

    columns = set().union(*counts.values()) if counts else set()
    good = sorted(c for c in columns
                  if sum(counts[f][c] / members[f] for f in members) / len(members)
                  >= CORE_MIN_FRACTION)
    blocks = []
    for col in good:
        if blocks and col - blocks[-1][1] <= BLOCK_GAP + 1:
            blocks[-1][1] = col
        else:
            blocks.append([col, col])

    merged = blocks[:1]
    for block in blocks[1:]:
        same = other = 0
        for segments in segment_of.values():
            a = {k for c, k in segments.items() if merged[-1][0] <= c <= merged[-1][1]}
            b = {k for c, k in segments.items() if block[0] <= c <= block[1]}
            if a and b:
                if a & b:
                    same += 1
                else:
                    other += 1
        if same > other:
            merged[-1][1] = block[1]
        else:
            merged.append(block)
    return [tuple(b) for b in merged]


def strand_name(block_index):
    """Name of the core strand of a block (0-based): b1, then b3, b4, ... (b2 is the hairpin)."""
    return "b1" if block_index == 0 else f"b{block_index + 2}"


def _span(segments):
    return min(seg[0] for seg in segments), max(seg[1] for seg in segments)


def _bridges_between(segment, other, bridges, parallel):
    """Number of bridges of one direction between two residue ranges."""
    return sum(1 for i, j, label in bridges if label.islower() == parallel
               and ((segment[0] <= i <= segment[1] and other[0] <= j <= other[1])
                    or (segment[0] <= j <= segment[1] and other[0] <= i <= other[1])))


def assign_strands(seq, ss8, sheet, n_residues, blocks, bridges, ladders=None):
    """Core strands of one structure as {"b1": (first, last, sheet label, pieces), ...} with
    0-based residue indices; pieces are the strand's (first, last) runs, pieces at most
    PIECE_GAP residues apart merged (a strand broken by a bulge DSSP does not bridge).
    bridges: [(residue index, partner index, ladder label)] (DSSP; uppercase = antiparallel);
    ladders: {(strand, strand): parallel} from core_ladders, for the pairing check."""
    cols = _residue_columns(seq)
    strands = defaultdict(list)     # strand name -> runs of E
    unnamed = []
    for first, last in _segments(ss8, n_residues):
        overlap = Counter(k for i in range(first, last + 1)
                          for k, (c0, c1) in enumerate(blocks) if c0 <= cols[i] <= c1)
        if overlap:
            strands[strand_name(overlap.most_common(1)[0][0])].append((first, last))
        else:
            unnamed.append((first, last))

    if strands.get("b1"):
        b1 = _span(strands["b1"])
        for segment in unnamed:
            if (1 <= segment[0] - b1[1] - 1 <= HAIRPIN_MAX_GAP
                    and _bridges_between(segment, b1, bridges, parallel=False)):
                strands["b2"].append(segment)
                unnamed.remove(segment)
                break

    for k in range(1, CORE_STRANDS + 1):
        name = f"b{k}"
        partners = {(b if a == name else a): parallel for (a, b), parallel in (ladders or {}).items()
                    if name in (a, b) and (b if a == name else a) in strands}
        if name == "b2" or not partners:
            continue
        if name in strands and any(
                _bridges_between(_span(strands[name]), _span(strands[p]), bridges, parallel)
                for p, parallel in partners.items()):
            continue
        before = max((_span(v)[1] for n, v in strands.items() if int(n[1:]) < k), default=-1)
        after = min((_span(v)[0] for n, v in strands.items() if int(n[1:]) > k),
                    default=n_residues)

        def hits(segment):
            return sum(_bridges_between(segment, _span(strands[p]), bridges, parallel)
                       for p, parallel in partners.items())

        candidates = [seg for seg in unnamed if before < seg[0] and seg[1] < after]
        best = max(((hits(seg), seg) for seg in candidates), default=(0, None))
        if best[0]:
            unnamed.extend(strands.pop(name, []))
            strands[name] = [best[1]]
            unnamed.remove(best[1])

    names = sorted(strands, key=lambda n: int(n[1:]))
    for name, following in zip(names, names[1:] + [None]):
        partners = {(b if a == name else a): parallel for (a, b), parallel in (ladders or {}).items()
                    if name in (a, b) and (b if a == name else a) in strands}
        end = _span(strands[name])[1]
        stop = _span(strands[following])[0] if following else n_residues
        for segment in [seg for seg in unnamed if end < seg[0] and seg[1] < stop]:
            if any(_bridges_between(segment, _span(strands[p]), bridges, parallel)
                   for p, parallel in partners.items()):
                strands[name].append(segment)
                unnamed.remove(segment)

    result = {}
    for name in names:
        pieces = []
        for seg in sorted(strands[name]):
            if pieces and seg[0] - pieces[-1][1] - 1 <= PIECE_GAP:
                pieces[-1] = (pieces[-1][0], max(pieces[-1][1], seg[1]))
            else:
                pieces.append(seg)
        label = Counter(sheet[i] for seg in pieces for i in range(seg[0], seg[1] + 1)
                        if sheet[i] != "-").most_common(1)
        result[name] = (pieces[0][0], pieces[-1][1], label[0][0] if label else "-", pieces)
    return result


def core_ladders(strands, bridges):
    """Ladders of the core: {(strand, strand): parallel} for the pairs bridged in one direction
    in at least LADDER_MIN_FRACTION of the structures having both strands.
    strands: {struct_id: assign_strands result}; bridges: {struct_id: bridges}."""
    both, bridged = Counter(), Counter()
    for struct_id, st in strands.items():
        for a, b in combinations(sorted(st, key=lambda n: int(n[1:])), 2):
            both[(a, b)] += 1
            for parallel in (True, False):
                if _bridges_between(st[a][:2], st[b][:2], bridges[struct_id], parallel):
                    bridged[(a, b, parallel)] += 1
    return {(a, b): parallel for (a, b, parallel), n in bridged.items()
            if n >= LADDER_MIN_FRACTION * both[(a, b)]}


def loop_regions(strands, n_residues):
    """7b: loops of one structure as {loop: (first, last)} (0-based residue indices, up to the
    catalytic domain end n_residues); a loop without residues is left out.
    strands: {name: (first, last, ...)} from assign_strands."""
    numbers = sorted(int(name[1:]) for name in strands)
    loops = {}
    for loop, (before, after) in LOOPS.items():
        earlier = [k for k in numbers if k <= int(before[1:])]
        first = strands[f"b{earlier[-1]}"][1] + 1 if earlier else 0
        later = [k for k in numbers if after and k >= int(after[1:])]
        last = strands[f"b{later[0]}"][0] - 1 if later else n_residues - 1
        if last - first + 1 >= (LC_MIN_LENGTH if after is None else 1):
            loops[loop] = (first, last)
    return loops
