#!/usr/bin/env python
"""
Standalone analysis: across ALL proteins, how well do the per-residue flexibility methods
agree? Computes Spearman rank correlation between every pair of methods per protein, then
aggregates (mean +/- std, n proteins). Answers "which method correlates most with MD /
B-factor / pLDDT?".

NOT part of the web app. Run directly:
    python lpmocat_pipeline/analysis/flexibility_correlation_stats.py

pLDDT is INVERTED (100 - pLDDT) so it points the same direction as the flexibility methods
(higher = more flexible / less confident); correlations are then positive when methods agree.

Outputs into lpmocat_pipeline/analysis/flexibility_comparison/:
  correlation_stats.png   mean Spearman-correlation matrix across all proteins (annotated with n)
"""
import os, sys, json
from pathlib import Path
from collections import defaultdict
from itertools import combinations

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2] / "lpmocat_web"))
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "config.settings")
import django
django.setup()

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from catalogue.models import Protein
from catalogue.makePlots import FLEX_METHOD_ORDER, FLEX_METHOD_META

OUT = HERE.parent / "flexibility_comparison"
OUT.mkdir(parents=True, exist_ok=True)


def load_methods(protein):
    out = {}
    for fx in protein.flexibility.all():
        try:
            vals = np.array([float(v) for v in json.loads(fx.rmsf)["values"]], dtype=float)
        except (ValueError, KeyError, TypeError):
            continue
        if fx.method == "plddt":
            vals = 100.0 - vals          # invert -> higher = more flexible/disordered
        out[fx.method] = vals
    return out


def main():
    pair_rhos = defaultdict(list)     # (A,B) sorted -> [rho per protein]
    n_proteins = 0
    for p in Protein.objects.all():
        methods = load_methods(p)
        if len(methods) < 2:
            continue
        n_proteins += 1
        for a, b in combinations(methods, 2):
            va, vb = methods[a], methods[b]
            n = min(len(va), len(vb))
            if n < 5:
                continue
            rho = spearmanr(va[:n], vb[:n]).statistic
            if rho is not None and not np.isnan(rho):
                pair_rhos[tuple(sorted((a, b)))].append(rho)

    methods = [m for m in FLEX_METHOD_ORDER
               if any(m in pair for pair in pair_rhos)]
    print(f"Proteins with >=2 methods: {n_proteins}")
    print("\nPairwise Spearman correlation (pLDDT inverted), across proteins:")
    print("  %-9s %-9s %7s %7s %6s" % ("method A", "method B", "mean", "std", "n"))
    mean_mat = np.full((len(methods), len(methods)), np.nan)
    for i, a in enumerate(methods):
        mean_mat[i, i] = 1.0
    for (a, b), rhos in sorted(pair_rhos.items(), key=lambda kv: -np.mean(kv[1])):
        arr = np.array(rhos)
        print("  %-9s %-9s %7.2f %7.2f %6d" % (a, b, arr.mean(), arr.std(), len(arr)))
        i, j = methods.index(a), methods.index(b)
        mean_mat[i, j] = mean_mat[j, i] = arr.mean()

    # Rank candidates by correlation with each reference method
    for ref in ["md_rmsf", "bfactor", "plddt"]:
        if ref not in methods:
            continue
        ranked = []
        for other in methods:
            if other == ref:
                continue
            key = tuple(sorted((ref, other)))
            if key in pair_rhos:
                arr = np.array(pair_rhos[key])
                ranked.append((arr.mean(), arr.std(), len(arr), other))
        ranked.sort(reverse=True)
        label = FLEX_METHOD_META.get(ref, {}).get("label", ref)
        print(f"\nBest agreement with {label}{' (inverted)' if ref=='plddt' else ''}:")
        for mean, std, n, other in ranked:
            print(f"    {other:9s} rho={mean:+.2f} +/- {std:.2f}  (n={n})")

    # heatmap
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(mean_mat, vmin=-1, vmax=1, cmap="RdBu_r")
    ax.set_xticks(range(len(methods))); ax.set_xticklabels(methods, rotation=45, ha="right")
    ax.set_yticks(range(len(methods))); ax.set_yticklabels(methods)
    for i in range(len(methods)):
        for j in range(len(methods)):
            if not np.isnan(mean_mat[i, j]):
                ax.text(j, i, f"{mean_mat[i,j]:.2f}", ha="center", va="center", fontsize=9,
                        color="white" if abs(mean_mat[i, j]) > 0.6 else "black")
    ax.set_title(f"Mean Spearman correlation across {n_proteins} proteins\n(pLDDT inverted)", fontsize=10)
    fig.colorbar(im, fraction=0.046)
    plt.tight_layout()
    f = OUT / "correlation_stats.png"; plt.savefig(f, dpi=140); plt.close()
    print(f"\nWrote {f}")


if __name__ == "__main__":
    main()
