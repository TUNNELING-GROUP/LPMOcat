#!/usr/bin/env python
"""
Standalone analysis: find proteins where GNM agrees POORLY with the other flexibility
methods, and test whether restricting to the catalytic domain (residues <= end_residue)
improves it — GNM assigns huge fluctuations to floppy termini/linkers, which may drag the
full-protein correlation down even when GNM is fine inside the domain.

    python lpmocat_pipeline/analysis/gnm_diagnostics.py

Prints, per protein: GNM's Spearman rho vs each other method, over the FULL protein and over
the CATALYTIC domain, sorted worst-first. Flags proteins where GNM rho < 0.3 with an
experimental/MD reference (B-factor / MD RMSF), or vs pLDDT for AlphaFold models.
"""
import os, sys, json
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2] / "lpmocat_web"))
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "config.settings")
import django
django.setup()

import numpy as np
from scipy.stats import spearmanr
from catalogue.models import Protein
from catalogue.makePlots import parse_pdb_residue_numbers, FLEX_METHOD_ORDER

bounds = json.load(open(HERE.parents[2] / "data" / "06_catalytic_domain" / "catalytic_domain_boundaries.json"))
POOR = 0.3   # rho below this vs a reference = poor agreement


def load(protein):
    out = {}
    for fx in protein.flexibility.all():
        try:
            vals = np.array([float(v) for v in json.loads(fx.rmsf)["values"]], dtype=float)
        except (ValueError, KeyError, TypeError):
            continue
        if fx.method == "plddt":
            vals = 100.0 - vals
        out[fx.method] = vals
    return {m: out[m] for m in FLEX_METHOD_ORDER if m in out}


def rho(a, b, mask=None):
    n = min(len(a), len(b))
    a, b = a[:n], b[:n]
    if mask is not None:
        m = mask[:n]
        a, b = a[m], b[m]
    if len(a) < 5 or np.ptp(a) == 0 or np.ptp(b) == 0:
        return None
    r = spearmanr(a, b).statistic
    return None if (r is None or np.isnan(r)) else r


rows = []
for p in Protein.objects.all():
    methods = load(p)
    if "gnm" not in methods or len(methods) < 2:
        continue
    end = bounds.get(p.selected_structure, {}).get("end_residue")
    res = parse_pdb_residue_numbers(p.pdb_structure)
    n = len(methods["gnm"])
    if len(res) != n:
        res = list(range(1, n + 1))
    cat_mask = np.array([r <= end for r in res]) if end is not None else np.ones(n, bool)

    others = [m for m in methods if m != "gnm"]
    full = {m: rho(methods["gnm"], methods[m]) for m in others}
    cat = {m: rho(methods["gnm"], methods[m], cat_mask) for m in others}
    # reference = experimental/MD if present, else pLDDT
    refs = [m for m in ("md_rmsf", "bfactor", "plddt") if m in others]
    ref_full = np.mean([full[m] for m in refs if full[m] is not None]) if refs else None
    rows.append((ref_full if ref_full is not None else 1.0, p.protein_code, p.selected_structure,
                 end, others, full, cat))

rows.sort(key=lambda r: r[0])   # worst reference-agreement first

print("GNM per-residue agreement, worst-first (rho vs each method; FULL protein | CATALYTIC only):")
print("  %-10s %-13s %5s  %s" % ("code", "struct", "end", "gnm vs ... (full | catalytic)"))
n_poor = 0
for ref_full, code, sid, end, others, full, cat in rows:
    parts = []
    for m in others:
        f = f"{full[m]:+.2f}" if full[m] is not None else "  na"
        c = f"{cat[m]:+.2f}" if cat[m] is not None else "  na"
        parts.append(f"{m}={f}|{c}")
    flag = ""
    exp = [m for m in ("md_rmsf", "bfactor", "plddt") if m in others and full[m] is not None]
    if exp and min(full[m] for m in exp) < POOR:
        flag = "  <-- POOR"; n_poor += 1
    print(f"  {code:10s} {sid:13s} {str(end):>5}  " + "  ".join(parts) + flag)

print(f"\n{n_poor} proteins with GNM rho < {POOR} vs an experimental/MD/pLDDT reference (full protein).")
print("Look at the FULL vs CATALYTIC columns: if catalytic is much better, the floppy tail is the culprit.")
