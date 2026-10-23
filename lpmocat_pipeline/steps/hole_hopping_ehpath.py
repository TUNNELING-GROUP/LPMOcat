"""
Step 12 of the LPMOcat pipeline: hole hopping pathways from the protein surface to the copper
with EHPath, using the scripts in tools/hole_hopping.

For every aligned structure with a histidine brace (step 6f):
  12a copper: the active-site metal; without a metal, the pseudo-Cu of step 6f (the centroid
      of the brace nitrogens: His-1 N, ND1, NE2 and His-2 ND1, NE2)
  12b classify the residues of the catalytic domain (up to the step-5 end residue):
      surface-exposed residues (relative SASA of the whole residue >= 10 %, step 5's values)
      are donors, Trp/Tyr/Cys/Met bridges and the Cu the
      acceptor (01_classify_resids.py)
  12c EHPath for every donor (EHPath/EHPath.py, as in 02_get_EH_results.py)
  12d paths ranked by mean residence time (03_analyse_paths.py) -> <id>_paths.txt
"""

import contextlib
import glob
import hashlib
import importlib.util
import io
import json
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd

from paths import TOOLS_DIR

HH_SCRIPTS_DIR = os.path.join(TOOLS_DIR, "hole_hopping")

# Structures excluded after manual inspection
EXCLUDED_STRUCT_IDS = {
    "3UAM": "BpAA10A: the first His of the brace is misplaced, so the Cu cannot be positioned",
}

# Per-structure donor thresholds (relative SASA, %): none. Until 2026-10-01 three structures had a
# lower one (9HDG 9, 6Z5Y 7, 8S6S 7), set by hand because the former side-chain measure left them
# without a donor at 10 %; with step 5's whole-residue values they have 1, 5 and 6 (user: 10 %
# everywhere, no exceptions).
SASA_THRESHOLD_OVERRIDES = {}

_METALS = ("CU", "ZN", "NI", "CO", "MN", "FE", "MO", "MG")   # = step 3g's BRACE_METALS
_SCRIPT_FILES = ("01_classify_resids.py", "02_get_EH_results.py", "03_analyse_paths.py",
                 os.path.join("EHPath", "EHPath.py"))
_modules = {}


def _script(name):
    """One of the numbered hole_hopping scripts, imported as a module."""
    if name not in _modules:
        spec = importlib.util.spec_from_file_location(
            name.lstrip("0123456789_"), os.path.join(HH_SCRIPTS_DIR, f"{name}.py"))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        _modules[name] = module
    return _modules[name]


def scripts_fingerprint():
    """Hash of the hole hopping scripts, so that changing them invalidates earlier results."""
    h = hashlib.sha256()
    for name in _SCRIPT_FILES:
        with open(os.path.join(HH_SCRIPTS_DIR, name), "rb") as f:
            h.update(f.read())
    return h.hexdigest()


def input_fingerprint(pdb_path, his1, his2, end_residue, sasa_threshold, scripts_sha, rsa):
    """Hash of everything a structure's hole hopping result depends on."""
    h = hashlib.sha256()
    with open(pdb_path, "rb") as f:
        h.update(f.read())
    h.update(json.dumps([his1, his2, end_residue, sasa_threshold, scripts_sha, rsa]).encode())
    return h.hexdigest()[:16]


def inject_metal_as_cu(pdb_path, cu_coord, out_dir, drop_metals=False):
    """Path of a copy of the structure with a CU HETATM (residue 999) before the first
    CONECT/MASTER/END record -- Biopython stops parsing at MASTER. drop_metals: remove the
    structure's own metal HETATMs first (a metal step 6f did not take as the copper)."""
    x, y, z = (float(v) for v in cu_coord)
    cu_line = f"HETATM 9999  CU   CU A 999    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          CU\n"
    with open(pdb_path) as f:
        lines = f.readlines()
    if drop_metals:
        lines = [line for line in lines if not (line.startswith("HETATM")
                 and line[17:20].strip().upper().rstrip("0123456789") in _METALS)]
    insert_at = next((i for i, line in enumerate(lines)
                      if line.startswith(("CONECT", "MASTER", "END"))), len(lines))
    lines.insert(insert_at, cu_line)
    out_path = os.path.join(out_dir, os.path.splitext(os.path.basename(pdb_path))[0]
                            + "_pseudocu.pdb")
    with open(out_path, "w") as f:
        f.writelines(lines)
    return out_path


def parse_paths(paths_file):
    """The paths ranked by exact mean residence time in a <id>_paths.txt file:
    [{"rank", "path", "mean_residence_time", "donor_residue", "donor_distance"}]."""
    paths = []
    in_exact = False
    with open(paths_file) as f:
        for line in f:
            if "Paths ranked according to exact mean residence time" in line:
                in_exact = True
                continue
            if not in_exact or line.startswith("---") or not line.strip():
                continue
            m = re.match(
                r"\s*(\d+)\.\s+\[([^\]]+)\]\s+([\d.eE+\-]+)\s+\(donor:\s+(\w+),\s+dist:\s+([\d.]+)",
                line)
            if m:
                paths.append({
                    "rank": int(m.group(1)),
                    "path": [r.strip() for r in m.group(2).split(",")],
                    "mean_residence_time": float(m.group(3)),
                    "donor_residue": m.group(4),
                    "donor_distance": float(m.group(5)),
                })
    return paths


def clear_outputs(struct_id, hh_dir):
    """Remove a structure's files of an earlier run."""
    stale = [os.path.join(hh_dir, f"{struct_id}_paths.txt"),
             os.path.join(hh_dir, f"{struct_id}_pseudocu.pdb")]
    for sub in ("EH_files", "EH_results"):
        stale += glob.glob(os.path.join(hh_dir, sub, glob.escape(struct_id) + "_*"))
    for path in stale:
        if os.path.exists(path):
            os.remove(path)


def _run(struct_id, pdb_path, end_residue, hh_dir, sasa_threshold, rsa, copper):
    """rsa: step 5's relative SASA (%) of the structure, as (residue numbers, values);
    copper: step 6f's copper position [x, y, z] (the structure's metal or the pseudo-Cu), or
    None. The copper is the structure's metal (any of _METALS, as find_copper reads it) when it
    sits at that position, else a CU injected there (the pseudo-Cu)."""
    classify = _script("01_classify_resids")
    get_results = _script("02_get_EH_results")
    analyse = _script("03_analyse_paths")

    clear_outputs(struct_id, hh_dir)
    eh_files_dir = os.path.join(hh_dir, "EH_files")
    eh_results_dir = os.path.join(hh_dir, "EH_results")
    os.makedirs(eh_files_dir, exist_ok=True)
    os.makedirs(eh_results_dir, exist_ok=True)
    aligned_dir = os.path.dirname(pdb_path)

    # 12a-12b: copper and residue classification
    try:
        structure = classify.load_structure(pdb_path)
        sasa_df = pd.DataFrame({"resid": rsa[0], "Relative SASA": rsa[1]})
        model = list(structure.get_models())[0]
        if copper is None:
            print("no copper position (step 6f)")
            return None
        cu_res, cu_coord = classify.find_copper(model)
        if cu_coord is None or np.linalg.norm(np.asarray(cu_coord) - np.asarray(copper)) > 0.01:
            if cu_coord is not None:
                print("its metal is not held by the brace: the pseudo-Cu instead")
            else:
                print("no metal: pseudo-Cu at the brace nitrogens")
            pdb_path = inject_metal_as_cu(pdb_path, copper, hh_dir, drop_metals=cu_coord is not None)
            structure = classify.load_structure(pdb_path)
            model = list(structure.get_models())[0]
            cu_res, cu_coord = classify.find_copper(model)
            if cu_coord is None:
                print("Cu injection failed")
                return None
        his_brace = classify.find_histidine_brace(model, cu_coord)
        classified = classify.classify_residues(model, sasa_df, cu_coord, his_brace, cu_res,
                                                end_residue=end_residue,
                                                sasa_threshold=sasa_threshold)
        classify.print_EH_input(classified, f"{struct_id}.pdb", eh_files_dir)
    except Exception as e:
        print(f"classification failed: {e}")
        return None

    # 12c: EHPath per donor
    donor_files = get_results.get_donor_files(eh_files_dir, f"{struct_id}.pdb")
    if not donor_files:
        print("no donor residues")
        return None
    ehpath = os.path.join(HH_SCRIPTS_DIR, "EHPath", "EHPath.py")
    bridge_csv = os.path.join(eh_files_dir, f"{struct_id}_bridge.csv")
    acceptor_csv = os.path.join(eh_files_dir, f"{struct_id}_acceptor.csv")
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    for donor_csv, resname, resid in donor_files:
        try:
            result = subprocess.run([sys.executable, ehpath, donor_csv, bridge_csv, acceptor_csv, resid],
                                    capture_output=True, text=True, timeout=300,
                                    cwd=HH_SCRIPTS_DIR, env=env)
        except subprocess.TimeoutExpired:
            print(f"EHPath timed out for donor {resname}{resid}")
            continue
        if result.returncode == 0:
            with open(os.path.join(eh_results_dir, f"{struct_id}_donor_{resname}{resid}.log"), "w") as f:
                f.write(result.stdout)
        else:
            print(f"EHPath failed for donor {resname}{resid}: {result.stderr.strip()[-200:]}")

    # 12d: rank the paths
    try:
        analyse.analyze_eh_results(eh_results_dir=eh_results_dir, output_dir=hh_dir,
                                   eh_files_dir=eh_files_dir, pdb_dir=aligned_dir,
                                   protein_filter=struct_id)
    except Exception as e:
        print(f"path analysis failed: {e}")
    paths_file = os.path.join(hh_dir, f"{struct_id}_paths.txt")
    return parse_paths(paths_file) if os.path.exists(paths_file) else None


def run_structure(args):
    """12a-12d for one structure (a multiprocessing worker).

    args: (struct_id, aligned pdb path, His-1, His-2, catalytic domain end residue, output
    directory, SASA threshold or None, step 5's relative SASA as (residue numbers, values),
    step 6f's pseudo-Cu or None). Returns (struct_id, paths or None, log text)."""
    struct_id, pdb_path, his1, his2, end_residue, hh_dir, sasa_threshold, rsa, copper = args
    log = io.StringIO()
    with contextlib.redirect_stdout(log):
        try:
            paths = _run(struct_id, pdb_path, end_residue, hh_dir, sasa_threshold, rsa, copper)
        except Exception as e:
            print(f"failed: {e}")
            paths = None
    return struct_id, paths, log.getvalue()
