#!/usr/bin/env python
"""
Standalone, throwaway analysis: compare the per-residue flexibility methods against each
other for a few proteins, to sanity-check the data (especially GNM, whose values look odd).

NOT part of the web app. Run it directly; it sets up Django to read the Flexibility rows.

Usage:
    python lpmocat_pipeline/analysis/compare_flexibility_methods.py [CODE ...]

With no args it uses a default set of crystal structures that have all four methods
(bfactor / cabsflex / gnm / md_rmsf) so GNM can be compared to experiment + MD.

Outputs (PNG) into lpmocat_pipeline/analysis/flexibility_comparison/:
  <CODE>_curves.png   overlaid per-residue curves, each method min-max normalised 0..1
  correlation.png     Spearman rank-correlation matrix between methods, per protein
"""
import os, sys, json
from pathlib import Path

# --- Django setup (read-only DB access) -----------------------------------------------
HERE = Path(__file__).resolve()
DJANGO_DIR = HERE.parents[2] / "lpmocat_web"      # dir containing manage.py / settings pkg
sys.path.insert(0, str(DJANGO_DIR))
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "config.settings")
import django
django.setup()

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from catalogue.models import Protein
from catalogue.makePlots import parse_pdb_residue_numbers, FLEX_METHOD_ORDER, FLEX_METHOD_META

OUT = HERE.parent / "flexibility_comparison"
OUT.mkdir(parents=True, exist_ok=True)

DEFAULTS = ["TrAA9", "TmAA10", "TdAA15A", "TtAA9E"]   # all have bfactor+cabsflex+gnm+md_rmsf


def load_methods(protein):
    """method -> np.array of per-residue values (in structure/residue order)."""
    out = {}
    for fx in protein.flexibility.all():
        try:
            vals = json.loads(fx.rmsf)["values"]
        except (ValueError, KeyError, TypeError):
            continue
        out[fx.method] = np.array([float(v) for v in vals], dtype=float)
    # keep display order
    return {m: out[m] for m in FLEX_METHOD_ORDER if m in out}


def minmax(a):
    lo, hi = np.nanmin(a), np.nanmax(a)
    return (a - lo) / (hi - lo) if hi > lo else np.zeros_like(a)


def main(codes):
    proteins = []
    for code in codes:
        p = Protein.objects.filter(protein_code=code).first()
        if not p:
            print(f"  ! {code} not found"); continue
        methods = load_methods(p)
        if len(methods) < 2:
            print(f"  ! {code} has <2 methods ({list(methods)}) — skipping"); continue
        # auth residue numbers for the x-axis (fall back to 1..N)
        res = parse_pdb_residue_numbers(p.pdb_structure)
        n = min([len(v) for v in methods.values()])
        if len(res) != n:
            res = list(range(1, n + 1))
        methods = {m: v[:n] for m, v in methods.items()}
        proteins.append((code, np.array(res[:n]), methods))
        print(f"  {code:9s} methods={list(methods)}  n_res={n}")

    # ---- per-protein overlaid normalised curves --------------------------------------
    for code, res, methods in proteins:
        plt.figure(figsize=(11, 3.2))
        for m, v in methods.items():
            plt.plot(res, minmax(v), label=FLEX_METHOD_META.get(m, {}).get("label", m), lw=1.3, alpha=0.85)
        plt.title(f"{code} — per-residue flexibility (each method min-max normalised)")
        plt.xlabel("residue number"); plt.ylabel("normalised 0–1")
        plt.legend(fontsize=8, ncol=len(methods)); plt.tight_layout()
        f = OUT / f"{code}_curves.png"; plt.savefig(f, dpi=130); plt.close()
        print(f"    wrote {f.name}")

    # ---- Spearman correlation matrix per protein -------------------------------------
    ncol = len(proteins)
    fig, axes = plt.subplots(1, ncol, figsize=(3.5 * ncol, 3.6), squeeze=False)
    for ax, (code, res, methods) in zip(axes[0], proteins):
        ms = list(methods)
        M = np.array([methods[m] for m in ms])
        rho = np.ones((len(ms), len(ms)))
        for i in range(len(ms)):
            for j in range(len(ms)):
                rho[i, j] = spearmanr(M[i], M[j]).statistic
        im = ax.imshow(rho, vmin=-1, vmax=1, cmap="RdBu_r")
        ax.set_xticks(range(len(ms))); ax.set_xticklabels(ms, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(len(ms))); ax.set_yticklabels(ms, fontsize=8)
        for i in range(len(ms)):
            for j in range(len(ms)):
                ax.text(j, i, f"{rho[i,j]:.2f}", ha="center", va="center", fontsize=7,
                        color="white" if abs(rho[i, j]) > 0.6 else "black")
        ax.set_title(code, fontsize=9)
    fig.suptitle("Spearman rank correlation between flexibility methods (per residue)", fontsize=10)
    fig.colorbar(im, ax=axes[0].tolist(), fraction=0.02)
    f = OUT / "correlation.png"; plt.savefig(f, dpi=130, bbox_inches="tight"); plt.close()
    print(f"    wrote {f.name}")

    # ---- text summary: how well GNM correlates with the others -----------------------
    print("\nGNM vs other methods (Spearman rho):")
    for code, res, methods in proteins:
        if "gnm" not in methods:
            continue
        line = "  %-9s " % code
        for m in methods:
            if m == "gnm":
                continue
            r = spearmanr(methods["gnm"], methods[m]).statistic
            line += f"{m}={r:+.2f}  "
        print(line)
    print(f"\nDone. Figures in {OUT}")


if __name__ == "__main__":
    codes = sys.argv[1:] or DEFAULTS
    print("Comparing flexibility methods for:", codes)
    main(codes)
