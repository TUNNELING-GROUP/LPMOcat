"""
Step 5 of the LPMOcat pipeline: per-residue solvent accessibility and secondary structure of
the aligned structures (04_alignment/aligned/<id>.pdb); step 6 uses the secondary structure too.

  5a  SASA: solvent-accessible surface area of every residue (Shrake-Rupley, Biopython),
      computed on the protein atoms only -- the active-site metal is left out -- and the
      relative SASA: SASA / the residue type's maximum (Tien et al. 2013, PLoS One 8:e80635,
      theoretical values in a Gly-X-Gly tripeptide; Biopython's 'Wilke' scale), in %
  5b  secondary structure: DSSP (mkdssp 4.4) with beta bulges joining strand segments, sheet
      labels and bridge partners; 3-state H (H, G, I) / E (E, B) / - as Mol* draws them

mkdssp (BSD-2, https://github.com/PDB-REDO/dssp) lives in the conda env 'dssp'
(conda create -n dssp -c conda-forge dssp=4.4.11); set MKDSSP to use another binary. It is run
in a temporary directory holding its mmCIF dictionary and tools/dssp/components.cif
(the 20 amino acids of the wwPDB Chemical Component Dictionary): the conda-forge build does
not find these data files elsewhere.
"""

import functools
import json
import os
import shutil
import subprocess
import tempfile

from paths import TOOLS_DIR
DSSP_ENV = "dssp"
COMPONENTS_CIF = os.path.join(TOOLS_DIR, "dssp", "components.cif")
STANDARD_RESIDUES = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
                     "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
THREE_STATE = {"H": "H", "G": "H", "I": "H", "E": "E", "B": "E"}
# 5a: maximum SASA per residue type (A^2) for the relative SASA
from Bio.Data.PDBData import residue_sasa_scales

MAX_ASA = residue_sasa_scales["Wilke"]


def relative_sasa(name, sasa):
    """Relative SASA (%, 2 decimals) of a residue type's SASA; None for a type without a maximum."""
    return round(100 * sasa / MAX_ASA[name], 2) if name in MAX_ASA else None


def residue_sasa(pdb_path):
    """5a: SASA (A^2) and relative SASA (%) of every polymer residue of the first model, as
    [(residue number, chain, residue name, SASA, relative SASA)]. Hetero residues (the metal)
    are removed before the calculation so that they do not bury part of the surface."""
    from Bio.PDB import PDBParser
    from Bio.PDB.SASA import ShrakeRupley

    structure = PDBParser(QUIET=True).get_structure("X", pdb_path)
    model = list(structure.get_models())[0]
    for other in list(structure)[1:]:
        structure.detach_child(other.id)
    for chain in model:
        for res in [r for r in chain if r.id[0].strip()]:
            chain.detach_child(res.id)
    ShrakeRupley().compute(model, level="R")
    return [(res.id[1], chain.id, res.get_resname(), res.sasa,
             relative_sasa(res.get_resname(), res.sasa))
            for chain in model for res in chain]


@functools.lru_cache(maxsize=None)
def mkdssp_path():
    """The mkdssp binary: $MKDSSP, else mkdssp on PATH, else the one of the conda env
    DSSP_ENV; None if there is none."""
    path = os.environ.get("MKDSSP") or shutil.which("mkdssp")
    if not path and shutil.which("conda"):
        try:
            r = subprocess.run(["conda", "env", "list", "--json"], capture_output=True,
                               text=True, timeout=60)
            envs = json.loads(r.stdout)["envs"]
        except (OSError, ValueError, KeyError, subprocess.TimeoutExpired):
            envs = []
        prefix = next((e for e in envs if os.path.basename(e) == DSSP_ENV), None)
        if prefix:
            path = os.path.join(prefix, "bin", "mkdssp")
    return path if path and os.path.exists(path) else None


def mkdssp_version():
    """'mkdssp version 4.4.11', or None when mkdssp is not installed."""
    exe = mkdssp_path()
    if not exe:
        return None
    r = subprocess.run([exe, "--version"], capture_output=True, text=True, timeout=60)
    lines = r.stdout.strip().splitlines()
    return lines[0] if lines else None


def _parse_dssp(text):
    """Residue lines of a classic DSSP file as
    [(DSSP number, chain, residue number, insertion code, state, sheet label,
      [(partner DSSP number, bridge label), ...])]; chain break lines ('!') are skipped."""
    rows, started = [], False
    for line in text.splitlines():
        if line.startswith("  #  RESIDUE"):
            started = True
            continue
        if not started or len(line) < 38 or line[13] == "!":
            continue
        partners = [(int(line[25:29]), line[23]), (int(line[29:33]), line[24])]
        rows.append((int(line[0:5]), line[11], int(line[5:10]), line[10].strip(),
                     line[16].strip() or "-", line[33].strip() or "-",
                     [(p, label) for p, label in partners if p]))
    return rows


def run_dssp(pdb_path):
    """5b: DSSP of the standard amino acids of a structure, as
    {"residues": [(chain, number, insertion code)], "ss8": "..." (DSSP states, '-' for none),
     "sheet": "..." (sheet labels, '-' for none),
     "bridges": [[residue, partner residue, ladder label], ...] (each pair once; lowercase
     label = parallel, uppercase = antiparallel)}.
    Residues DSSP skips (incomplete backbone) are '-'. None when mkdssp is not installed;
    RuntimeError when it fails."""
    from Bio.PDB import PDBParser

    exe = mkdssp_path()
    if exe is None:
        return None
    model = next(PDBParser(QUIET=True).get_structure("x", pdb_path).get_models())
    residues = [(chain.id, res.id[1], res.id[2].strip()) for chain in model for res in chain
                if res.id[0] == " " and res.get_resname() in STANDARD_RESIDUES]
    with open(pdb_path) as f:
        atoms = [line for line in f
                 if line.startswith("ATOM  ") and line[17:20] in STANDARD_RESIDUES]

    dictionary = os.path.join(os.path.dirname(os.path.realpath(exe)), os.pardir, "share",
                              "libcifpp", "mmcif_pdbx.dic")
    with tempfile.TemporaryDirectory() as tmp:
        if os.path.exists(dictionary):
            os.symlink(os.path.realpath(dictionary), os.path.join(tmp, "mmcif_pdbx.dic"))
        os.symlink(COMPONENTS_CIF, os.path.join(tmp, "components.cif"))
        with open(os.path.join(tmp, "in.pdb"), "w") as f:
            f.writelines(atoms)
            f.write("END\n")
        r = subprocess.run([exe, "--output-format", "dssp", "in.pdb", "out.dssp"], cwd=tmp,
                           capture_output=True, text=True, timeout=300)
        out = os.path.join(tmp, "out.dssp")
        if r.returncode != 0 or not os.path.exists(out):
            message = (r.stderr or r.stdout).strip().splitlines()
            raise RuntimeError(f"mkdssp: {message[0] if message else f'exit code {r.returncode}'}")
        with open(out) as f:
            rows = _parse_dssp(f.read())

    by_number = {row[0]: (row[1], row[2], row[3]) for row in rows}
    state = {(row[1], row[2], row[3]): (row[4], row[5]) for row in rows}
    ss8 = "".join(state.get(res, ("-", "-"))[0] for res in residues)
    sheet = "".join(state.get(res, ("-", "-"))[1] for res in residues)
    bridges = set()
    for row in rows:
        for partner, label in row[6]:
            pair = sorted([row[2], by_number[partner][1]])
            bridges.add((pair[0], pair[1], label))
    return {"residues": residues, "ss8": ss8, "sheet": sheet,
            "bridges": [list(b) for b in sorted(bridges)]}


def three_state(ss8):
    """DSSP states reduced to H (H, G, I), E (E, B) and - (the rest), as Mol* draws them."""
    return [THREE_STATE.get(s, "-") for s in ss8]


def secondary_structure(pdb_path):
    """Per-residue 3-state secondary structure (see three_state) of the standard amino acids,
    or None when mkdssp is not installed."""
    dssp = run_dssp(pdb_path)
    return None if dssp is None else three_state(dssp["ss8"])

