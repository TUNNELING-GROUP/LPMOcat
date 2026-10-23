"""
Descriptive similarity / conservation analysis for a *selected* set of LPMOs, powering the
compare page's sequence logo and similarity heatmap. Same computations as the offline exploration
(lpmocat_pipeline/analysis/binding_surface_similarity.py) but on the user's selection, where
the shared alignment columns are much denser.

Given a REGION (whole catalytic / loops / core / hole-hopping / a named loop) and a
per-residue PROPERTY (charge / hydrophobicity / flexibility / rsa), compute():
  A) a protein x protein similarity matrix + a clustered ordering
  B) per-position tracks: Schneider information content (amino-acid conservation, gap = 21st
     symbol) and property variance (divergence), plus loop/hole-hopping flags
over the US-align structural-alignment columns the selection shares.

This is a descriptive comparison aid, not a predictor of function.
"""
import ast
import json
import math
from collections import Counter

import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram

from .makePlots import CATALYTIC_BOUNDARIES
from .surface import SURFACE_RSA

# Continuous per-residue properties (numeric -> z-scored Euclidean similarity).
CONTINUOUS = ("charge", "hydrophobicity", "flexibility", "rsa")
# Categorical properties (per-residue label -> fraction-matching similarity).
CATEGORICAL = ("sequence", "aa_groups")
PROPERTIES = CONTINUOUS + CATEGORICAL
PROPERTY_LABELS = {
    "charge": "Charge", "hydrophobicity": "Hydrophobicity",
    "flexibility": "Flexibility", "rsa": "Relative SASA",
    "sequence": "Sequence (identity)", "aa_groups": "Amino-acid groups",
}
REGION_LABELS = {
    "whole": "Whole catalytic domain", "loops": "All loops",
    "core": "Core (non-loop)", "hh": "Hole-hopping residues",
    "binding": "Binding surface", "binding_aromatic": "Aromatic binding residues",
    "second": "Second-sphere residues",
}
# Physicochemical classes for the "amino-acid groups" property.
AA_GROUPS = {
    "A": "hydrophobic", "V": "hydrophobic", "L": "hydrophobic", "I": "hydrophobic",
    "M": "hydrophobic", "C": "hydrophobic",
    "F": "aromatic", "W": "aromatic", "Y": "aromatic",
    "S": "polar", "T": "polar", "N": "polar", "Q": "polar",
    "K": "positive", "R": "positive", "H": "positive",
    "D": "negative", "E": "negative",
    "G": "special", "P": "special",
}
# Relative SASA (%) at/above which a residue counts as surface-exposed: the catalogue's one
# threshold (surface.SURFACE_RSA, 25 %); it was 20 % here until 2026-10-06
EXPOSED_MIN = SURFACE_RSA


def available_loops(records):
    """Loop names present across the selected proteins (for the region dropdown)."""
    names = set()
    for r in records:
        if r.loop_regions:
            try:
                names.update(ast.literal_eval(r.loop_regions).keys())
            except Exception:
                pass
    return sorted(names)


def _extract(pm, records, flex_method):
    """Per-protein per-column property arrays + region info (reuses makePlots outputs)."""
    ncol = len(pm.oneLetterSequence[0]) if pm.oneLetterSequence else 0
    method = flex_method if flex_method in pm.allFlexibilityLists else pm.defaultFlexMethod
    data = {}
    for i, r in enumerate(records):
        seq = pm.oneLetterSequence[i]

        def per_col(lst):
            return np.array([np.nan if v is None else float(v) for v in lst], dtype=float)

        props = {
            "hydrophobicity": per_col(pm.hydrophobicityList[i]),
            "charge":         per_col(pm.electrostaticsList[i]),
            "flexibility":    per_col(pm.allFlexibilityLists[method][i]),
        }
        rsa = np.full(ncol, np.nan)
        rl, k = pm.relSASAList[i], 0
        for c, aa in enumerate(seq):
            if aa != "-":
                if k < len(rl):
                    rsa[c] = rl[k]
                k += 1
        props["rsa"] = rsa

        # Two numberings meet here: col2auth holds the structure's own residue numbers (what
        # hole hopping, the brace and the catalytic end use), col2pos the 1-based sequence
        # position (what loop_regions, beta_core and binding_surface use). 49 of 79 AA9
        # structures do not start at residue 1, so mixing them shifts every loop.
        col2auth, col2pos, k, rn = {}, {}, 0, pm.residueNumbers[i]
        for c, aa in enumerate(seq):
            if aa != "-":
                col2auth[c] = rn[k] if k < len(rn) else None
                col2pos[c] = k + 1
                k += 1

        cat_end = (CATALYTIC_BOUNDARIES.get((r.selected_structure or "").strip()) or {}).get("end_residue")
        loop_pos = _ranges_by_name(r.loop_regions)
        core_pos = _ranges_by_name(r.beta_core)
        bind_pos = _ranges_by_name(r.binding_surface)
        hh = set()
        if r.hole_hopping and r.hole_hopping.paths_json:
            try:
                for p in json.loads(r.hole_hopping.paths_json):
                    for res in p.get("path", []):
                        if res.upper().startswith("CU"):
                            continue
                        digits = "".join(ch for ch in res if ch.isdigit())
                        if digits:
                            hh.add(int(digits))
            except Exception:
                pass

        data[r.protein_code] = {"seq": seq, "props": props, "col2auth": col2auth,
                                "col2pos": col2pos, "cat_end": cat_end, "loop_pos": loop_pos,
                                "core_pos": core_pos, "bind_pos": bind_pos,
                                "site_auth": _site_residues(r), "hh_auth": hh,
                                "second_auth": _second_sphere_residues(r),
                                "aromatic_auth": _binding_aromatic_residues(r)}
    return data, ncol, method


def _ranges_by_name(field):
    """{name: {sequence position, ...}} from a Protein field holding "start-end" ranges
    (loop_regions, beta_core, binding_surface). These are 1-based positions, not residue
    numbers; see _extract."""
    out = {}
    if not field:
        return out
    try:
        groups = ast.literal_eval(field) if str(field).lstrip().startswith("{'") else json.loads(field)
    except Exception:
        return out
    for name, ranges in (groups or {}).items():
        s = set()
        for rng in ranges or []:
            try:
                a, b = str(rng).split("-")
                s.update(range(int(a), int(b) + 1))
            except ValueError:
                continue
        if s:
            out[name] = s
    return out


# Active-site categories shown as their own lane, by the structure's residue numbers. They
# split the same way the compare page's Sequence plot colours them, so a category is the same
# colour in both: the two brace histidines, the aromatic on the copper's axial position, the
# rest of the second-sphere motif, and the shielding residue on the solvent side.
SITE_CATEGORIES = ("N-His", "Internal-His", "Axial aromatic", "Second sphere", "Shielding")


def _binding_aromatic_residues(protein):
    """Residue numbers of the aromatic binding residues (overview.binding_aromatics: Phe, Trp, Tyr,
    His on the binding surface, the brace histidines left out)."""
    from .overview import protein_binding_aromatics
    try:
        return {a["residue"] for a in (protein_binding_aromatics(protein) or [])}
    except Exception:
        return set()


def _second_sphere_residues(protein):
    """Residue numbers of the second sphere, the set the 3D viewer shows as "Second-sphere
    residues": the hydrogen-bonding motif with the axial aromatic, the shielding residue and the
    side chains in contact with a brace histidine (H-bond, pi-stacking, cation-pi). Not the brace
    histidines themselves, nor the Gly/Pro/Arg of the brace motifs."""
    from .second_sphere import second_sphere
    out = set()
    try:
        site = second_sphere(protein)
        for res in (site.get("residues") or []) + [site.get("distal_residue") or {}]:
            if res.get("residue"):
                out.add(int(res["residue"]))
    except Exception:
        pass
    try:
        contacts = json.loads(protein.histidine_brace_interactions or "{}")
        for partners in (contacts.values() if isinstance(contacts, dict) else []):
            for c in partners or []:
                if c.get("other_res_id") is not None:
                    out.add(int(c["other_res_id"]))
    except (TypeError, ValueError):
        pass
    try:
        brace = ast.literal_eval(protein.histidine_brace) if protein.histidine_brace else ()
        out -= {int(n) for n in brace if n}
    except (ValueError, SyntaxError, TypeError):
        pass
    return out


def _site_residues(protein):
    """{residue number: category} for the brace histidines, the axial aromatic, the other
    second-sphere motif residues and the shielding residue. These fields store residue
    numbers, not sequence positions."""
    from .brace_motifs import brace_motifs
    from .second_sphere import second_sphere
    out = {}
    try:
        brace = ast.literal_eval(protein.histidine_brace) if protein.histidine_brace else None
        if brace:
            for num, label in zip(brace, ("N-His", "Internal-His")):
                if num:
                    out[int(num)] = label
    except (ValueError, SyntaxError):
        pass
    try:
        site = second_sphere(protein)
        for res in site.get("residues") or []:
            if res.get("residue"):
                out.setdefault(int(res["residue"]),
                               "Axial aromatic" if res.get("role") == "aromatic" else "Second sphere")
        shield = site.get("distal_residue")
        if shield and shield.get("residue"):
            out[int(shield["residue"])] = "Shielding"
    except Exception:
        pass
    try:
        for motif in brace_motifs(protein).values():
            for res in motif.get("residues") or []:
                if res.get("residue"):
                    out.setdefault(int(res["residue"]), "Second sphere")
    except Exception:
        pass
    return out


def _in_region(d, col, region, exposed=False):
    auth = d["col2auth"].get(col)
    if auth is None:
        return False
    if d["cat_end"] is not None and auth > d["cat_end"]:
        return False
    if exposed:                                   # surface-exposed residues only
        rsa = d["props"]["rsa"][col]
        if not (rsa == rsa) or rsa < EXPOSED_MIN:  # NaN or below threshold
            return False
    if region == "whole":
        in_r = True
    else:
        pos = d["col2pos"].get(col)             # loops are sequence positions, not residue numbers
        is_loop = pos is not None and any(pos in s for s in d["loop_pos"].values())
        if region == "loops":
            in_r = is_loop
        elif region == "core":
            in_r = not is_loop
        elif region == "hh":
            in_r = auth in d["hh_auth"]
        elif region == "second":
            in_r = auth in d["second_auth"]
        elif region == "binding_aromatic":
            in_r = auth in d["aromatic_auth"]
        elif region == "binding":
            in_r = pos is not None and any(pos in s for s in d["bind_pos"].values())
        else:
            in_r = pos is not None and pos in d["loop_pos"].get(region, set())
    return in_r


def _schneider_ic(chars):
    """Schneider information content (bits); gap counted as a 21st symbol."""
    n = len(chars)
    if not n:
        return 0.0
    H = 0.0
    for cnt in Counter(chars).values():
        p = cnt / n
        H -= p * math.log2(p)
    s = 21
    en = (s - 1) / (2 * n * math.log(2))
    return round(max(0.0, math.log2(s) - (H + en)), 4)


LOG2_20 = math.log2(20)   # max information content per column (bits), gaps excluded


def _logo_column(chars):
    """Frequency-logo stack for a column: [[aa, fraction], ...] descending, where fraction =
    (count / TOTAL selected proteins), i.e. the letter height IS the residue's percentage, so
    it matches the composition shown on hover. Gaps leave empty space at the top of the stack
    (the total = occupancy), and a conserved column is one full-height letter while a variable
    one is several short letters."""
    n_total = len(chars)
    aas = [c for c in chars if c not in ('-', 'X', 'x', 'B', 'Z') and c.isalpha()]
    if n_total == 0 or not aas:
        return []
    stack = [[aa, round(cnt / n_total, 3)]
             for aa, cnt in Counter(aas).items() if cnt / n_total > 0.01]
    stack.sort(key=lambda x: -x[1])
    return stack


def _composition(chars):
    """Full residue composition of a column as [[symbol, pct], ...] descending, incl. 'gap',
    over ALL selected proteins (so percentages sum to 100)."""
    total = len(chars)
    if not total:
        return []
    out = [['gap' if sym == '-' else sym, round(100 * cnt / total, 1)]
           for sym, cnt in Counter(chars).items()]
    out.sort(key=lambda x: -x[1])
    return out


def _column_members(data, codes, keep):
    """Per kept column, which selected proteins carry each residue and at which of their own
    residue numbers: [{symbol: [[code, residue number | None], ...]}], 'gap' for absent. The
    logo's hover reads this, so a letter names its proteins with the number each one uses (they
    are numbered differently) and can highlight them in the selection tree."""
    out = []
    for c in keep:
        by_symbol = {}
        for cc in codes:
            aa = data[cc]["seq"][c]
            by_symbol.setdefault("gap" if aa == "-" else aa, []).append(
                [cc, data[cc]["col2auth"].get(c)])
        out.append(dict(sorted(by_symbol.items(), key=lambda kv: (-len(kv[1]), kv[0]))))
    return out


def _lane_names(data, codes):
    """The lanes to draw under the logo: the loops present in the selection (drawn as one
    lane, coloured per loop) and the active-site categories. Core strands are counted only
    as a whole ("core_any"), not per strand."""
    loops = sorted({n for cc in codes for n in data[cc]["loop_pos"]})
    return {"loops": loops, "site": list(SITE_CATEGORIES)}


def _column_lanes(data, codes, keep):
    """Per kept column, how many selected proteins have it in each region — the conservation
    of the regions themselves, alongside the residue conservation the logo shows. Counts and
    the member codes, so a lane hovers like a logo letter."""
    names = _lane_names(data, codes)
    lanes = {"binding": [], "core_any": [], "loop_any": []}
    for key in names["loops"]:
        lanes["loop:" + key] = []
    for key in names["site"]:
        lanes["site:" + key] = []

    for c in keep:
        hit = {k: [] for k in lanes}
        for cc in codes:
            d = data[cc]
            pos, auth = d["col2pos"].get(c), d["col2auth"].get(c)
            if pos is None:
                continue
            if any(pos in s for s in d["bind_pos"].values()):
                hit["binding"].append(cc)
            in_loop = False
            for name, s in d["loop_pos"].items():
                if pos in s:
                    in_loop = True
                    if "loop:" + name in hit:
                        hit["loop:" + name].append(cc)
            if in_loop:
                hit["loop_any"].append(cc)
            if any(pos in s for s in d["core_pos"].values()):
                hit["core_any"].append(cc)
            category = d["site_auth"].get(auth)
            if category and "site:" + category in hit:
                hit["site:" + category].append(cc)
        for key, members in hit.items():
            lanes[key].append(members)
    return lanes


def _cluster_order(D):
    """Leaf order (indices into codes) from average-linkage on the distance matrix."""
    n = D.shape[0]
    kept = [i for i in range(n) if np.isfinite(D[i]).sum() > 1]
    if len(kept) < 3:
        return list(range(n))
    sub = D[np.ix_(kept, kept)].copy()
    fm = np.nanmax(sub[np.isfinite(sub)]) if np.isfinite(sub).any() else 1.0
    sub[~np.isfinite(sub)] = fm * 1.25
    np.fill_diagonal(sub, 0.0)
    link = linkage(sub[np.triu_indices(len(kept), 1)], method="average")
    leaves = dendrogram(link, no_plot=True)["leaves"]
    ordered = [kept[l] for l in leaves]
    ordered += [i for i in range(n) if i not in kept]      # any dropped proteins at the end
    return ordered


def _round_or_none(v):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return None
    return round(float(v), 4)


def _similarity_from_distance(D, absolute=False):
    """Clustered order + a 0..1 similarity matrix from a distance matrix.

    absolute=True (categorical distances, already 0..1) -> similarity = 1 - D, i.e. the true
    fraction of matching positions. absolute=False (continuous, unbounded RMS z-diff) ->
    1 - D/max, a similarity relative to the most-distant pair in the set."""
    order = _cluster_order(D)
    n = D.shape[0]
    if absolute:
        S = 1.0 - D
    else:
        dmax = np.nanmax(D[np.isfinite(D)]) if np.isfinite(D).any() else 1.0
        S = 1.0 - D / dmax if dmax > 0 else np.ones_like(D)
    return order, [[_round_or_none(S[a, b]) for b in range(n)] for a in range(n)]


def compute(records, msa, region, prop, flex_method="cabsflex", min_occ=0.3, exposed=False):
    """Return the similarity matrix + clustered order (A) and per-position tracks (B).

    Continuous properties (charge/hydrophobicity/flexibility/rsa) use z-scored Euclidean
    similarity; categorical ones (sequence/aa_groups) use fraction-matching (identity /
    same-class). `exposed` restricts the region to surface-exposed residues (relSASA>=EXPOSED_MIN).
    """
    from .makePlots import makePlots
    if prop not in PROPERTIES:
        prop = "charge"
    categorical = prop in CATEGORICAL
    pm = makePlots(records, msa, no_trim=True)
    data, ncol, method = _extract(pm, records, flex_method)
    codes = [r.protein_code for r in records]
    n = len(codes)

    # A column enters the analysis if it is in the region for AT LEAST ONE selected protein
    # (so sparse networks like hole-hopping keep their single-protein positions) AND enough
    # proteins have a residue there (drops rare-insertion noise). Every protein's residue at
    # those columns is then compared — so presence/absence counts, not only positions where
    # both proteins are in the region.
    # min_occ <= 0 keeps every column of the region, so the sequence logo spans exactly the
    # columns the compare page's plots show (identical sets) and the two line up position for
    # position. The similarity matrix keeps the default floor, where a column only one or two
    # proteins occupy is noise rather than signal.
    need_res = 1 if min_occ <= 0 else max(2, math.ceil(min_occ * n))
    keep = [c for c in range(ncol)
            if sum(data[cc]["seq"][c] != '-' for cc in codes) >= need_res
            and any(_in_region(data[cc], c, region, exposed) for cc in codes)]

    D = np.full((n, n), np.nan)
    np.fill_diagonal(D, 0.0)
    # per pair: the positions compared (the heatmap's tooltip) and, for categorical properties,
    # how many of them match
    positions = [[None] * n for _ in range(n)]
    matching = [[None] * n for _ in range(n)]
    if categorical:
        # Per-protein token per kept column: amino acid (or its class), or "GAP" where absent.
        T = [[None] * len(keep) for _ in range(n)]
        for i, cc in enumerate(codes):
            seq = data[cc]["seq"]
            for j, c in enumerate(keep):
                aa = seq[c]
                T[i][j] = "GAP" if aa == "-" else (AA_GROUPS.get(aa, "other") if prop == "aa_groups" else aa)
        # Distance = fraction of columns (where at least one has a residue) that DON'T match.
        # gap-vs-residue counts as a mismatch; gap-vs-gap is skipped as uninformative.
        for a in range(n):
            for b in range(a + 1, n):
                same = tot = 0
                for j in range(len(keep)):
                    ta, tb = T[a][j], T[b][j]
                    if ta == "GAP" and tb == "GAP":
                        continue
                    tot += 1
                    if ta == tb:
                        same += 1
                positions[a][b] = positions[b][a] = tot
                matching[a][b] = matching[b][a] = same
                if tot >= 5:
                    D[a, b] = D[b, a] = 1.0 - same / tot
    else:
        # Property matrix, z-scored per column within the selection -> RMS z-difference distance.
        M = np.full((n, len(keep)), np.nan)
        for i, cc in enumerate(codes):
            arr = data[cc]["props"][prop]
            for j, c in enumerate(keep):
                M[i, j] = arr[c]
        Z = M.copy()
        for j in range(len(keep)):
            col = M[:, j]
            sd = np.nanstd(col)
            Z[:, j] = (col - np.nanmean(col)) / sd if sd > 1e-9 else 0.0
        for a in range(n):
            for b in range(a + 1, n):
                m = ~np.isnan(Z[a]) & ~np.isnan(Z[b])
                positions[a][b] = positions[b][a] = int(m.sum())
                if m.sum() >= 5:
                    D[a, b] = D[b, a] = float(np.sqrt(np.mean((Z[a, m] - Z[b, m]) ** 2)))

    order, similarity = _similarity_from_distance(D, absolute=categorical)

    # B) per-position tracks over the region columns.
    ic, pstd, occ, is_loop, is_hh, logo, comp = [], [], [], [], [], [], []
    for c in keep:
        col_chars = [data[cc]["seq"][c] for cc in codes]
        ic.append(_schneider_ic(col_chars))
        logo.append(_logo_column(col_chars))
        comp.append(_composition(col_chars))
        if categorical:
            pstd.append(None)                        # spread is not defined for a label
        else:
            vals = np.array([data[cc]["props"][prop][c] for cc in codes], dtype=float)
            vals = vals[~np.isnan(vals)]
            pstd.append(_round_or_none(np.std(vals)) if len(vals) >= 2 else None)
        occ.append(round(float(np.mean([data[cc]["seq"][c] != "-" for cc in codes])), 3))
        is_loop.append(bool(np.mean([_in_region(data[cc], c, "loops") for cc in codes]) >= 0.5))
        is_hh.append(bool(np.mean([_in_region(data[cc], c, "hh") for cc in codes]) >= 0.15))

    # A column has no single residue number: the proteins are numbered differently (some start
    # at 1, others at 16, 18, 20 ...), so the modal number this used to report repeated itself
    # and ran backwards. Positions are the column itself (1..n_cols, in `members` order) and
    # each protein's own residue number travels with it in `members`.
    col_labels = [c + 1 for c in range(len(keep))]

    members, lanes = _column_members(data, codes, keep), _column_lanes(data, codes, keep)

    return {
        "ok": True,
        "codes": codes,
        "order": order,
        "similarity": similarity,
        # per pair: positions compared, matching positions (categorical), the distance itself
        # (continuous: RMS z-difference; similarity = 1 - distance / the largest one)
        "positions": positions, "matching": matching,
        "distance": [[_round_or_none(D[a, b]) for b in range(n)] for a in range(n)],
        "n_cols": len(keep),
        "col_labels": col_labels,
        "members": members, "lanes": lanes, "lane_names": _lane_names(data, codes),
        "ic": ic, "pstd": pstd, "occ": occ, "is_loop": is_loop, "is_hh": is_hh,
        "logo": logo, "logo_max_bits": round(LOG2_20, 3), "comp": comp,
        "region": region, "property": prop, "flex_method": method,
        "categorical": categorical, "exposed": bool(exposed),
        "property_label": PROPERTY_LABELS.get(prop, prop),
        "region_label": REGION_LABELS.get(region, region)
                        + (" · exposed only" if exposed else ""),
    }
