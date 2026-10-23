#!/usr/bin/env python
"""
Offline exploration (NOT the web app): property-based SIMILARITY + per-position
CONSERVATION over the US-align structural alignment, restricted to a chosen REGION of
the catalytic domain, across the whole catalogue.

Question being probed (see the design brainstorm): do the biophysical properties of the
LPMO binding-surface loops / hole-hopping network cluster proteins by substrate or
regioselectivity? If yes, it's worth a compare-page feature; if not, we stop here.

REGION  (which structural-alignment columns enter the analysis; all are restricted to the
         catalytic domain first):
    whole  - whole catalytic domain
    loops  - all loop-region residues (binding surface)
    core   - catalytic domain minus loops (scaffold)
    hh     - hole-hopping path residues
    L2/L3/LS/...  - one named loop

PROPERTY (per residue, using the same definitions the web plots use, pulled via makePlots):
    charge, hydrophobicity, flexibility, rsa

Two outputs per (region, property):
    A) protein x protein similarity + hierarchical clustering  -> clustermap PNG, rows
       colour-coded by substrate & regioselectivity, with a cluster/substrate purity print.
    B) per-position tracks: Schneider information content (amino-acid conservation, gap =
       21st symbol) and property variance -> both extremes are interesting (conserved =
       essential, divergent = specificity-determining).

Files land in lpmocat_pipeline/analysis/similarity_out/ (gitignored).

Usage:
    python lpmocat_pipeline/analysis/binding_surface_similarity.py               # default sweep
    python lpmocat_pipeline/analysis/binding_surface_similarity.py --region loops --property charge
    python lpmocat_pipeline/analysis/binding_surface_similarity.py --flex-method cabsflex --min-occ 0.3
"""
import os
import sys
import ast
import json
import math
import argparse
from pathlib import Path
from collections import Counter

# --- Django setup (read-only) ---------------------------------------------------------
HERE = Path(__file__).resolve()
DJANGO_DIR = HERE.parents[2] / "lpmocat_web"
sys.path.insert(0, str(DJANGO_DIR))
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "config.settings")
import django  # noqa: E402
django.setup()

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from scipy.cluster.hierarchy import linkage, dendrogram  # noqa: E402

from catalogue.models import Protein  # noqa: E402
from catalogue.makePlots import makePlots, CATALYTIC_BOUNDARIES  # noqa: E402

OUT = HERE.parent / "similarity_out"
OUT.mkdir(parents=True, exist_ok=True)
OUT_TAG = [""]   # filename prefix (set to "<family>_" for subset runs so they don't overwrite)

PROPERTIES = ["charge", "hydrophobicity", "flexibility", "rsa"]
SUBSTRATE_COLORS = {
    "Cellulose": "#2c7fb8", "Chitin": "#d95f0e", "Starch": "#31a354",
    "Xylan": "#756bb1", "Pectin": "#c51b8a", "Mannan": "#636363",
}
REGIO_COLORS = {
    "C1-hydroxylating": "#1b9e77", "C4-dehydrogenating": "#d95f02",
    "C1-hydroxylating/C4-dehydrogenating": "#7570b3", "Unknown": "#bdbdbd",
}


# ============================ data extraction =========================================
def load_dataset(flex_method, only=None, family=None):
    """makePlots over all proteins -> per-protein per-column property arrays + region info."""
    recs = [r for r in Protein.objects.all().order_by("protein_code") if r.aligned_sequence]
    if only:
        want = set(only)
        recs = [r for r in recs if r.protein_code in want]
    if family:
        recs = [r for r in recs if family.lower() in r.protein_code.lower()]
    msa = {"name": [r.protein_code for r in recs], "seq": [r.aligned_sequence for r in recs]}
    pm = makePlots(recs, msa, no_trim=True)
    ncol = len(pm.oneLetterSequence[0])
    method = flex_method if flex_method in pm.allFlexibilityLists else pm.defaultFlexMethod

    data = {}
    for i, r in enumerate(recs):
        seq = pm.oneLetterSequence[i]

        def per_col(lst):
            return np.array([np.nan if v is None else float(v) for v in lst], dtype=float)

        props = {
            "hydrophobicity": per_col(pm.hydrophobicityList[i]),
            "charge":         per_col(pm.electrostaticsList[i]),
            "flexibility":    per_col(pm.allFlexibilityLists[method][i]),
        }
        # relSASAList is per non-gap residue -> expand to per column
        rsa = np.full(ncol, np.nan)
        rlist, k = pm.relSASAList[i], 0
        for c, aa in enumerate(seq):
            if aa != "-":
                if k < len(rlist):
                    rsa[c] = rlist[k]
                k += 1
        props["rsa"] = rsa

        # column -> auth residue number
        col2auth, k, resnums = {}, 0, pm.residueNumbers[i]
        for c, aa in enumerate(seq):
            if aa != "-":
                col2auth[c] = resnums[k] if k < len(resnums) else None
                k += 1

        cat_end = (CATALYTIC_BOUNDARIES.get((r.selected_structure or "").strip()) or {}).get("end_residue")
        loop_auth = {}
        if r.loop_regions:
            try:
                for lname, ranges in ast.literal_eval(r.loop_regions).items():
                    s = set()
                    for rng in ranges:
                        a, b = str(rng).split("-")
                        s.update(range(int(a), int(b) + 1))
                    loop_auth[lname] = s
            except Exception:
                pass
        hh_auth = set()
        if r.hole_hopping and r.hole_hopping.paths_json:
            try:
                for path in json.loads(r.hole_hopping.paths_json):
                    for res in path.get("path", []):
                        if res.upper().startswith("CU"):
                            continue
                        digits = "".join(ch for ch in res if ch.isdigit())
                        if digits:
                            hh_auth.add(int(digits))
            except Exception:
                pass

        data[r.protein_code] = {
            "seq": seq, "props": props, "col2auth": col2auth, "cat_end": cat_end,
            "loop_auth": loop_auth, "hh_auth": hh_auth,
            "substrate": [str(s) for s in r.substrate.all()],
            "regio": [str(x) for x in r.regioselectivity.all()],
        }
    return [r.protein_code for r in recs], ncol, data, method


def in_region(d, col, region):
    """Is this protein's residue at column `col` part of `region` (within catalytic domain)?"""
    auth = d["col2auth"].get(col)
    if auth is None:
        return False
    if d["cat_end"] is not None and auth > d["cat_end"]:
        return False                                   # everything is catalytic-domain only
    if region == "whole":
        return True
    is_loop = any(auth in s for s in d["loop_auth"].values())
    if region == "loops":
        return is_loop
    if region == "core":
        return not is_loop
    if region == "hh":
        return auth in d["hh_auth"]
    return auth in d["loop_auth"].get(region, set())    # a named loop


# ============================ similarity + clustering =================================
def build_matrix(codes, ncol, data, region, prop, min_occ):
    """proteins x kept-columns property matrix (NaN where absent), z-scored per column."""
    M = np.full((len(codes), ncol), np.nan)
    for i, c in enumerate(codes):
        d, arr = data[c], data[c]["props"][prop]
        for col in range(ncol):
            if in_region(d, col, region):
                M[i, col] = arr[col]
    occ = np.mean(~np.isnan(M), axis=0)
    keep = np.where(occ >= min_occ)[0]
    M = M[:, keep]
    # z-score each column across the whole catalogue (fixed frame)
    Z = M.copy()
    for j in range(M.shape[1]):
        col = M[:, j]
        sd = np.nanstd(col)
        Z[:, j] = (col - np.nanmean(col)) / sd if sd > 1e-9 else 0.0
    return M, Z, keep


def pairwise_distance(Z, min_shared=5):
    """RMS z-difference over columns both proteins share (scale-free in #columns)."""
    n = Z.shape[0]
    D = np.full((n, n), np.nan)
    for a in range(n):
        D[a, a] = 0.0
        for b in range(a + 1, n):
            mask = ~np.isnan(Z[a]) & ~np.isnan(Z[b])
            if mask.sum() >= min_shared:
                diff = Z[a, mask] - Z[b, mask]
                D[a, b] = D[b, a] = float(np.sqrt(np.mean(diff ** 2)))
    return D


def cluster_order(D):
    """Drop proteins with no usable distances, fill remaining gaps, return (order, kept_idx, linkage)."""
    keep = [i for i in range(D.shape[0]) if np.isfinite(D[i]).sum() > 1]
    if len(keep) < 3:
        return keep, keep, None
    sub = D[np.ix_(keep, keep)].copy()
    finite_max = np.nanmax(sub[np.isfinite(sub)]) if np.isfinite(sub).any() else 1.0
    sub[~np.isfinite(sub)] = finite_max * 1.25         # unknown pairs = far apart
    np.fill_diagonal(sub, 0.0)
    condensed = sub[np.triu_indices(len(keep), 1)]
    Z = linkage(condensed, method="average")
    order = dendrogram(Z, no_plot=True)["leaves"]
    return [keep[o] for o in order], keep, Z


def primary(labels):
    return labels[0] if labels else "Unknown"


def clustermap(codes, data, Zmat, D, region, prop, method):
    row_order, kept_idx, link = cluster_order(D)
    if link is None:
        print(f"  [{region}/{prop}] too few proteins with data — skipped")
        return
    row_codes = [codes[i] for i in row_order]
    Zplot = Zmat[row_order]

    fig = plt.figure(figsize=(min(22, 5 + Zplot.shape[1] * 0.03), 2 + len(row_order) * 0.14))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.5, 12], wspace=0.02)

    # substrate + regio colour bars (+ protein labels)
    axc = fig.add_subplot(gs[0, 0])
    bar = np.ones((len(row_order), 2, 3))
    for r, i in enumerate(row_order):
        bar[r, 0] = matplotlib.colors.to_rgb(SUBSTRATE_COLORS.get(primary(data[codes[i]]["substrate"]), "#ffffff"))
        bar[r, 1] = matplotlib.colors.to_rgb(REGIO_COLORS.get(primary(data[codes[i]]["regio"]), "#ffffff"))
    axc.imshow(bar, aspect="auto")
    axc.set_xticks([0, 1]); axc.set_xticklabels(["sub", "regio"], fontsize=6, rotation=90)
    axc.set_yticks(range(len(row_order))); axc.set_yticklabels(row_codes, fontsize=4)

    # property heatmap (z-scored), rows in clustered order
    axh = fig.add_subplot(gs[0, 1])
    vlim = np.nanpercentile(np.abs(Zplot), 98) if np.isfinite(Zplot).any() else 2
    im = axh.imshow(np.ma.masked_invalid(Zplot), aspect="auto", cmap="RdBu_r",
                    vmin=-vlim, vmax=vlim, interpolation="nearest")
    axh.set_yticks([]); axh.set_xlabel(f"{region} columns  ({Zplot.shape[1]})", fontsize=8)
    fig.colorbar(im, ax=axh, fraction=0.015, pad=0.02, label=f"{prop} (z)")

    legend = ([Patch(fc=c, label=s) for s, c in SUBSTRATE_COLORS.items()]
              + [Patch(fc=c, label=r) for r, c in REGIO_COLORS.items()])
    axh.legend(handles=legend, fontsize=5, ncol=3, loc="upper left",
               bbox_to_anchor=(0, -0.04), frameon=False)

    kept_codes = [codes[i] for i in kept_idx]
    sub_p, sub_b = cluster_purity(link, kept_codes, data, "substrate")
    reg_p, reg_b = cluster_purity(link, kept_codes, data, "regio")
    fig.suptitle(f"{region}  ·  {prop}"
                 + (f"  (flex={method})" if prop == "flexibility" else "")
                 + f"    substrate {sub_p:.2f} (base {sub_b:.2f})"
                 + f"  ·  regio {reg_p:.2f} (base {reg_b:.2f})", fontsize=10)
    f = OUT / f"{OUT_TAG[0]}{region}_{prop}_clustermap.png"
    fig.savefig(f, dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"  [{region}/{prop}] {f.name}  |  n={len(row_order)} cols={Zplot.shape[1]}  "
          f"substrate={sub_p:.2f}/base{sub_b:.2f}  regio={reg_p:.2f}/base{reg_b:.2f}")


def cluster_purity(link, kept_codes, data, key, k=6):
    """Cut the tree into k clusters; fraction of proteins in their cluster's majority `key`
    label ('substrate' or 'regio'). Also returns the baseline (largest single class)."""
    from scipy.cluster.hierarchy import fcluster
    if len(kept_codes) <= k:
        k = max(2, len(kept_codes) // 2)
    labels = fcluster(link, k, criterion="maxclust")
    correct = total = 0
    for cl in set(labels):
        members = [kept_codes[i] for i in range(len(kept_codes)) if labels[i] == cl]
        cnt = Counter(primary(data[c][key]) for c in members)
        if cnt:
            correct += cnt.most_common(1)[0][1]
            total += len(members)
    overall = Counter(primary(data[c][key]) for c in kept_codes)
    baseline = overall.most_common(1)[0][1] / len(kept_codes) if kept_codes else 0.0
    return (correct / total if total else 0.0), baseline


# ============================ per-position conservation ================================
def schneider_ic(chars):
    """Schneider information content (bits) for a column; gap counted as a 21st symbol."""
    n = len(chars)
    if n == 0:
        return 0.0
    H = 0.0
    for cnt in Counter(chars).values():
        p = cnt / n
        H -= p * math.log2(p)
    s = 21
    en = (s - 1) / (2 * n * math.log(2))               # small-sample correction
    return max(0.0, math.log2(s) - (H + en))


def conservation_tracks(codes, ncol, data, prop, method, min_occ):
    """Per catalytic-domain column: Schneider IC (aa) + property std, shaded by region."""
    cat_cols = [c for c in range(ncol)
                if np.mean([in_region(data[cd], c, "whole") for cd in codes]) >= min_occ]
    if not cat_cols:
        print(f"  [tracks/{prop}] no catalytic columns above occupancy — skipped")
        return

    ic, pstd, occ, is_loop, is_hh = [], [], [], [], []
    for c in cat_cols:
        chars = [data[cd]["seq"][c] for cd in codes]
        ic.append(schneider_ic(chars))
        vals = np.array([data[cd]["props"][prop][c] for cd in codes
                         if in_region(data[cd], c, "whole")], dtype=float)
        vals = vals[~np.isnan(vals)]
        pstd.append(float(np.std(vals)) if len(vals) >= 2 else np.nan)
        occ.append(np.mean([data[cd]["seq"][c] != "-" for cd in codes]))
        is_loop.append(np.mean([in_region(data[cd], c, "loops") for cd in codes]) >= 0.5)
        is_hh.append(np.mean([in_region(data[cd], c, "hh") for cd in codes]) >= 0.15)
    x = np.arange(len(cat_cols))
    ic, pstd = np.array(ic), np.array(pstd)

    fig, axes = plt.subplots(2, 1, figsize=(min(26, 6 + len(cat_cols) * 0.03), 6), sharex=True)
    for ax, y, ylab, hi_hi in [(axes[0], ic, "Schneider IC (bits)\naa conservation", True),
                               (axes[1], pstd, f"{prop} std\n(divergence)", True)]:
        ax.bar(x, np.nan_to_num(y), width=1.0, color="#4a6b5a")
        for xi, lp in zip(x, is_loop):
            if lp:
                ax.axvspan(xi - 0.5, xi + 0.5, color="#f0c419", alpha=0.18, lw=0)
        for xi, hh in zip(x, is_hh):
            if hh:
                ax.axvspan(xi - 0.5, xi + 0.5, color="#e0559b", alpha=0.22, lw=0)
        ax.set_ylabel(ylab, fontsize=8)
        finite = y[np.isfinite(y)]
        if finite.size:
            thr = np.nanpercentile(finite, 90)
            for xi, yi in zip(x, y):
                if np.isfinite(yi) and yi >= thr:
                    ax.plot(xi, yi, ".", color="#c0392b", ms=3)
    axes[1].set_xlabel("catalytic-domain alignment columns (occupied)  ·  "
                       "yellow = loop, pink = hole-hopping  ·  red dots = top 10%", fontsize=8)
    axes[0].set_title(f"Per-position conservation / divergence — {prop}"
                      + (f" (flex={method})" if prop == "flexibility" else ""), fontsize=11)
    f = OUT / f"{OUT_TAG[0]}tracks_{prop}.png"
    fig.tight_layout(); fig.savefig(f, dpi=140); plt.close(fig)
    print(f"  [tracks/{prop}] {f.name}  ({len(cat_cols)} columns)")


# ============================ main ====================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--region", default=None, help="whole|loops|core|hh|<loop name>")
    ap.add_argument("--property", default=None, choices=PROPERTIES)
    ap.add_argument("--flex-method", default="cabsflex")
    ap.add_argument("--min-occ", type=float, default=0.3,
                    help="min fraction of proteins present in a column to keep it")
    ap.add_argument("--proteins", nargs="+", default=None,
                    help="restrict the whole analysis to these protein codes (denser alignment)")
    ap.add_argument("--family", default=None,
                    help="restrict to protein codes containing this substring, e.g. AA9")
    args = ap.parse_args()

    print("Loading dataset (makePlots over all proteins)...")
    codes, ncol, data, method = load_dataset(args.flex_method, only=args.proteins, family=args.family)
    if args.family:
        OUT_TAG[0] = args.family + "_"
    print(f"  {len(codes)} proteins, {ncol} alignment columns, flex method = {method}")

    if args.region and args.property:
        combos = [(args.region, args.property)]
        tracks = [args.property]
    else:
        combos = [("loops", p) for p in PROPERTIES] + [("hh", "charge"), ("whole", "charge")]
        tracks = PROPERTIES
        print("No region/property given -> default sweep:", [f"{r}/{p}" for r, p in combos])

    print("\nSimilarity + clustering:")
    for region, prop in combos:
        _, Z, _ = build_matrix(codes, ncol, data, region, prop, args.min_occ)
        D = pairwise_distance(Z)
        clustermap(codes, data, Z, D, region, prop, method)

    print("\nConservation tracks:")
    for prop in tracks:
        conservation_tracks(codes, ncol, data, prop, method, args.min_occ)

    print(f"\nDone. Figures in {OUT}")


if __name__ == "__main__":
    main()
