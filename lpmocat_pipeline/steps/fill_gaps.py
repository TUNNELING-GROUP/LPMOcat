"""
Fill unresolved internal loops of experimental structures with MODELLER.

Step 3a of the pipeline (`prepare`), before 3b-3i. For each selected PDB structure, the
polymer residues that the depositors list as unobserved in the selected chain (PDB REMARK 465 /
mmCIF _pdbx_unobs_or_zero_occ_residues) are rebuilt when they form an internal gap of at most
MAX_GAP_LENGTH residues. Missing termini and longer gaps are left as they are.

MODELLER builds only the missing residues: every observed atom is a fixed template atom
(AutoModel with select_atoms restricted to the gap residues). N_MODELS models are built and the
lowest-DOPE one is kept. Only the rebuilt residues are then spliced into the downloaded file, so
the deposited coordinates, B-factors, ligands (Cu, glycans), modified residues (HIC) and other
chains stay exactly as they were. Rebuilt atoms get occupancy 0.00 and B-factor 0.00 and are
listed in filled/modelled_residues.json.

MODELLER (academic licence) lives in its own conda env and runs in a subprocess through
modeller_worker.py: set MODELLER_PYTHON to that env's python, otherwise the python of the conda
env 'modeller' is used.
"""
import functools
import json
import os
import shutil
import subprocess

import numpy as np
from Bio.PDB import PDBParser, parse_pdb_header
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from preprocess_structures import cif_to_pdb_lines, keep_first_altloc, lookup_parent_residue

MAX_GAP_LENGTH = 30   # longer internal gaps (e.g. whole unresolved domains) are left unfilled
N_MODELS = 10
RAND_SEED = -8123     # MODELLER's default seed, fixed so that reruns rebuild the same loops

WORKER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "modeller_worker.py")
MODELLER_ENV = "modeller"

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
_ATOM_RECORDS = ("ATOM  ", "HETATM")
_ATOM_EXTRAS = ("ANISOU", "SIGATM", "SIGUIJ")   # per-atom records following an ATOM line


# ---------------------------------------------------------------------------
# MODELLER subprocess
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=None)
def _worker_cmd():
    # The env's own interpreter: `conda run -n modeller python` would pick up an active venv's
    # python first on PATH.
    python = os.environ.get("MODELLER_PYTHON")
    if not python and shutil.which("conda"):
        try:
            r = subprocess.run(["conda", "env", "list", "--json"], capture_output=True,
                               text=True, timeout=60)
            envs = json.loads(r.stdout)["envs"]
        except (OSError, ValueError, KeyError, subprocess.TimeoutExpired):
            envs = []
        prefix = next((e for e in envs if os.path.basename(e) == MODELLER_ENV), None)
        if prefix:
            python = os.path.join(prefix, "bin", "python")
    return (python, WORKER) if python and os.path.exists(python) else None


def modeller_version():
    """(version, None) if MODELLER runs with a valid licence, else (None, reason)."""
    cmd = _worker_cmd()
    if not cmd:
        return None, f"no MODELLER_PYTHON set and no conda env '{MODELLER_ENV}'"
    try:
        r = subprocess.run([*cmd, "--version"], capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired) as e:
        return None, str(e)
    if r.returncode != 0:
        lines = (r.stderr or r.stdout).strip().splitlines()
        return None, lines[-1] if lines else f"exit code {r.returncode}"
    return r.stdout.strip().splitlines()[-1], None


# ---------------------------------------------------------------------------
# Gap detection
# ---------------------------------------------------------------------------

def _res_key(line):
    return int(line[22:26]), line[26].strip()


def _fmt_key(key):
    return f"{key[0]}{key[1]}"


def _element(line):
    # Element column, else the first letter of the atom name (no standard heavy atom starts with H)
    return line[76:78].strip() or line[12:16].strip()[:1]


def _unobserved_residues(path, chain):
    """Polymer residues listed as unobserved in `chain`: [((resseq, icode), resname)]."""
    if path.endswith(".cif"):
        d = MMCIF2Dict(path)
        p = "_pdbx_unobs_or_zero_occ_residues."
        if p + "auth_seq_id" not in d:
            return []

        def col(name, default="?"):
            v = d.get(p + name, default)
            return v if isinstance(v, list) else [v] * len(d[p + "auth_seq_id"])

        out = []
        for asym, seq, ins, comp, poly, occ, model in zip(
                col("auth_asym_id"), col("auth_seq_id"), col("PDB_ins_code"),
                col("auth_comp_id"), col("polymer_flag"), col("occupancy_flag"),
                col("PDB_model_num", "1")):
            if asym == chain and poly == "Y" and occ == "1" and model in ("1", "?", "."):
                out.append(((int(seq), "" if ins in ("?", ".") else ins), comp))
        return out
    header = parse_pdb_header(path)
    return [((m["ssseq"], (m["insertion"] or "").strip()), m["res_name"])
            for m in header.get("missing_residues", [])
            if m["chain"] == chain and m["model"] in (None, 1)]


def _read_lines(path):
    """Structure as PDB-format lines (mmCIF converted with the pipeline's converter)."""
    if path.endswith(".cif"):
        return cif_to_pdb_lines(path)
    with open(path) as f:
        return f.readlines()


def _observed_polymer(lines, pdb_chain):
    """Polymer residues of the chain in file order (first model, first conformer per residue --
    see keep_first_altloc): [(key, resname, atom_lines)]. HETATM residues count as polymer when
    they have N, CA and C (modified amino acids such as HIC or MSE)."""
    residues = {}
    for line in keep_first_altloc(lines):
        if line.startswith("ENDMDL"):
            break
        if line[:6] not in _ATOM_RECORDS or line[21] != pdb_chain:
            continue
        resname = line[17:20].strip()
        if resname in ("HOH", "WAT", "DOD"):
            continue
        key = _res_key(line)
        entry = residues.setdefault(key, [resname, []])
        if entry[0] == resname:
            entry[1].append(line)
    out = []
    for key, (resname, atoms) in residues.items():
        names = {l[12:16].strip() for l in atoms}
        if any(l.startswith("ATOM") for l in atoms) or {"N", "CA", "C"} <= names:
            out.append((key, resname, atoms))
    return out


def _parent(resname):
    """Standard parent of a residue name, or None."""
    if resname in THREE_TO_ONE:
        return resname
    parent = lookup_parent_residue(resname)
    return parent if parent in THREE_TO_ONE else None


def plan_gaps(path, chain):
    """Which unobserved residues of the selected chain can be rebuilt.

    Returns dict with 'items' (observed and to-be-built residues in chain order:
    (key, resname, 'obs' | 'fill' | 'break')), 'fill' / 'skip' gap lists, the counts of missing
    terminal residues, and 'error' when the structure cannot be processed."""
    pdb_chain = chain[0]
    lines = _read_lines(path)
    observed = _observed_polymer(lines, pdb_chain)
    plan = {"lines": lines, "pdb_chain": pdb_chain, "observed": observed,
            "fill": [], "skip": [], "n_term_missing": 0, "c_term_missing": 0, "items": []}
    keys = [k for k, _, _ in observed]
    if keys != sorted(keys):
        plan["error"] = "residue numbering is not monotonic - not filled"
        return plan
    obs_keys = set(keys)
    missing = [(k, n) for k, n in _unobserved_residues(path, chain) if k not in obs_keys]
    if not observed or not missing:
        return plan

    merged = sorted([(k, n, "obs") for k, n, _ in observed] + [(k, n, "miss") for k, n in missing],
                    key=lambda r: r[0])
    first = next(i for i, r in enumerate(merged) if r[2] == "obs")
    last = max(i for i, r in enumerate(merged) if r[2] == "obs")
    plan["n_term_missing"] = first
    plan["c_term_missing"] = len(merged) - 1 - last
    merged = merged[first:last + 1]

    items, i = [], 0
    while i < len(merged):
        if merged[i][2] == "obs":
            items.append(merged[i])
            i += 1
            continue
        j = i
        while merged[j][2] == "miss":
            j += 1
        run = merged[i:j]
        gap = {"start": _fmt_key(run[0][0]), "end": _fmt_key(run[-1][0]), "length": len(run),
               "sequence": "".join(THREE_TO_ONE.get(_parent(n) or "", "X") for _, n, _ in run)}
        if len(run) > MAX_GAP_LENGTH:
            gap["reason"] = f"longer than {MAX_GAP_LENGTH} residues"
        elif "X" in gap["sequence"]:
            gap["reason"] = "unknown residue type"
        if "reason" in gap:
            plan["skip"].append(gap)
            items.append((None, None, "break"))
        else:
            plan["fill"].append(gap)
            items.extend((k, _parent(n), "fill") for k, n, _ in run)
        i = j
    plan["items"] = items
    return plan


# ---------------------------------------------------------------------------
# MODELLER input
# ---------------------------------------------------------------------------

def _write_modeller_inputs(plan, workdir):
    """Template PDB (observed residues, renumbered, one chain per segment between unfilled gaps)
    and the PIR alignment. Returns the 0-based target indices of the residues to build."""
    atoms_by_key = {k: (n, atoms) for k, n, atoms in plan["observed"]}
    tmpl_lines, tmpl_seq, target_seq, select = [], [], [], []
    chain_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    segment, resnum, serial, target_idx = 0, 0, 0, 0
    for key, resname, kind in plan["items"]:
        if kind == "break":
            segment += 1
            tmpl_seq.append("/")
            target_seq.append("/")
            continue
        if kind == "fill":
            tmpl_seq.append("-")
            target_seq.append(THREE_TO_ONE[resname])
            select.append(target_idx)
            target_idx += 1
            continue
        name, atoms = atoms_by_key[key]
        parent = _parent(name) or "ALA"
        resnum += 1
        for line in atoms:
            if _element(line) in ("H", "D"):
                continue
            serial += 1
            tmpl_lines.append(f"ATOM  {serial:>5d} {line[12:16]} {parent:>3s} "
                              f"{chain_ids[segment]}{resnum:>4d}    {line[30:54]}"
                              f"  1.00  0.00          {line[76:78]}\n")
        tmpl_seq.append(THREE_TO_ONE[parent])
        target_seq.append(THREE_TO_ONE[parent])
        target_idx += 1
    tmpl_lines.append("END\n")
    with open(os.path.join(workdir, "template.pdb"), "w") as f:
        f.writelines(tmpl_lines)

    def wrap(seq):
        s = "".join(seq) + "*"
        return "\n".join(s[i:i + 75] for i in range(0, len(s), 75))

    with open(os.path.join(workdir, "gaps.ali"), "w") as f:
        f.write(f">P1;template\nstructureX:template:1:A:{resnum}:{chain_ids[segment]}::::\n"
                f"{wrap(tmpl_seq)}\n>P1;target\nsequence:target::::::::\n{wrap(target_seq)}\n")
    return select


# ---------------------------------------------------------------------------
# Splicing the rebuilt residues into the downloaded file
# ---------------------------------------------------------------------------

def _atom_line(atom, resname, pdb_chain, key):
    x, y, z = atom.coord
    element = (atom.element or atom.get_name()[0]).strip()
    return (f"ATOM  {0:>5d} {atom.get_fullname()} {resname:>3s} {pdb_chain}{key[0]:>4d}"
            f"{key[1]:1s}   {x:8.3f}{y:8.3f}{z:8.3f}{0.0:6.2f}{0.0:6.2f}          {element:>2s}\n")


def _renumber(lines):
    """Sequential atom serials after the insertion; CONECT records remapped, MASTER dropped."""
    out, serial_map, n, last = [], {}, 0, 0
    for line in lines:
        rec = line[:6]
        if rec in _ATOM_RECORDS:
            n += 1
            old = line[6:11].strip()
            if old and old != "0":
                serial_map[old] = n
            line = f"{rec}{n % 100000:>5d}{line[11:]}"
            last = n
        elif rec.startswith("TER"):
            n += 1
            line = f"TER   {n % 100000:>5d}{line[11:]}" if len(line.rstrip("\n")) > 6 else "TER\n"
        elif rec in _ATOM_EXTRAS:
            line = f"{rec}{last % 100000:>5d}{line[11:]}"
        elif rec == "CONECT":
            body = line.rstrip("\n")
            fields = [body[i:i + 5].strip() for i in range(6, len(body), 5)]
            fields = [f for f in fields if f]
            if n >= 100000 or any(f not in serial_map for f in fields):
                continue
            line = "CONECT" + "".join(f"{serial_map[f]:>5d}" for f in fields) + "\n"
        elif rec == "MASTER":
            continue
        out.append(line)
    return out


def _splice(plan, model_residues, remarks):
    """Downloaded lines with the rebuilt residues inserted after their preceding residue."""
    lines, pdb_chain = plan["lines"], plan["pdb_chain"]
    polymer = {k: n for k, n, _ in plan["observed"]}

    # Rebuilt residues grouped by the observed residue they follow.
    inserts, prev_obs = {}, None
    for (key, resname, kind), res in zip([it for it in plan["items"] if it[2] != "break"],
                                         model_residues):
        if kind == "obs":
            prev_obs = key
        else:
            inserts.setdefault(prev_obs, []).extend(
                _atom_line(a, resname, pdb_chain, key) for a in res if a.element not in ("H", "D"))

    last_line, cur = {}, None
    for i, line in enumerate(lines):
        if line.startswith("ENDMDL"):
            break
        rec = line[:6]
        if rec in _ATOM_RECORDS and line[21] == pdb_chain:
            key = _res_key(line)
            cur = key if polymer.get(key) == line[17:20].strip() else None
            if cur is not None:
                last_line[cur] = i
        elif rec in _ATOM_EXTRAS and cur is not None:
            last_line[cur] = i
        else:
            cur = None
    insert_after = {last_line[k]: v for k, v in inserts.items()}

    out, remarks_done = [], False
    last_remark = max((i for i, l in enumerate(lines) if l.startswith("REMARK")), default=-1)
    if last_remark < 0:
        out.extend(remarks)
        remarks_done = True
    for i, line in enumerate(lines):
        out.append(line)
        if i == last_remark and not remarks_done:
            out.extend(remarks)
        if i in insert_after:
            out.extend(insert_after[i])
    return _renumber(out)


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

def _check(plan, model_residues):
    """Geometry sanity checks of the best model.

    fixed_max_dev: largest CA shift of an observed residue (should be ~0 -- they are fixed);
    junctions: C(i)-N(i+1) bond lengths where rebuilt residues meet observed ones (~1.33 A);
    min_contact: closest heavy-atom distance between rebuilt residues and non-neighbouring
    residues or ligands of the chain."""
    items = [it for it in plan["items"] if it[2] != "break"]
    obs_atoms = {k: {l[12:16].strip(): np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
                     for l in atoms} for k, _, atoms in plan["observed"]}
    devs, junctions = [], []
    for idx, ((key, _, kind), res) in enumerate(zip(items, model_residues)):
        if kind == "obs" and "CA" in res and "CA" in obs_atoms[key]:
            devs.append(float(np.linalg.norm(res["CA"].coord - obs_atoms[key]["CA"])))
        if idx and (kind == "fill") != (items[idx - 1][2] == "fill"):
            prev = model_residues[idx - 1]
            if "C" in prev and "N" in res:
                junctions.append(round(float(np.linalg.norm(res["N"].coord - prev["C"].coord)), 2))

    fill_idx = [i for i, it in enumerate(items) if it[2] == "fill"]
    near = set(fill_idx) | {i - 1 for i in fill_idx} | {i + 1 for i in fill_idx}
    built = np.array([a.coord for i in fill_idx for a in model_residues[i] if a.element != "H"])
    others = [a.coord for i, res in enumerate(model_residues) if i not in near
              for a in res if a.element != "H"]
    for line in plan["lines"]:   # ligands / metal of the chain (not water)
        if line.startswith("ENDMDL"):
            break
        if (line.startswith("HETATM") and line[21] == plan["pdb_chain"]
                and line[17:20].strip() not in ("HOH", "WAT", "DOD")
                and _res_key(line) not in obs_atoms):
            others.append(np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]))
    others = np.array(others)
    min_contact = None
    if len(built) and len(others):
        d = np.sqrt(((built[:, None, :] - others[None, :, :]) ** 2).sum(-1))
        min_contact = round(float(d.min()), 2)
    return {"fixed_max_dev": round(max(devs), 3) if devs else None,
            "junctions": junctions, "min_contact": min_contact}


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def fill_structure_gaps(path, chain, out_dir, version, plan=None):
    """Rebuild the fillable gaps of one structure. Writes out_dir/<id>.pdb and returns the
    record for modelled_residues.json, or None if there is nothing to fill."""
    plan = plan or plan_gaps(path, chain)
    if plan.get("error") or not plan["fill"]:
        return None
    struct_id = os.path.splitext(os.path.basename(path))[0]
    workdir = os.path.join(out_dir, "modeller", struct_id)
    shutil.rmtree(workdir, ignore_errors=True)
    os.makedirs(workdir)

    select = _write_modeller_inputs(plan, workdir)
    spec = {"workdir": os.path.abspath(workdir), "alignment": "gaps.ali",
            "select": select, "n_models": N_MODELS, "seed": RAND_SEED}
    with open(os.path.join(workdir, "spec.json"), "w") as f:
        json.dump(spec, f, indent=1)
    with open(os.path.join(workdir, "modeller.log"), "w") as log:
        r = subprocess.run([*_worker_cmd(), os.path.join(workdir, "spec.json")],
                           stdout=log, stderr=subprocess.STDOUT, timeout=3600)
    result_path = os.path.join(workdir, "result.json")
    if r.returncode != 0 or not os.path.exists(result_path):
        raise RuntimeError(f"MODELLER failed for {struct_id} - see {workdir}/modeller.log")
    with open(result_path) as f:
        models = [m for m in json.load(f)["models"] if not m["failure"]]
    if not models:
        raise RuntimeError(f"MODELLER built no model for {struct_id} - see {workdir}/modeller.log")
    best = min(models, key=lambda m: m["dope"])

    structure = PDBParser(QUIET=True).get_structure("m", os.path.join(workdir, best["name"]))
    model_residues = [r for r in structure[0].get_residues() if r.id[0] == " "]
    n_target = sum(1 for it in plan["items"] if it[2] != "break")
    if len(model_residues) != n_target:
        raise RuntimeError(f"{struct_id}: model has {len(model_residues)} residues, "
                           f"expected {n_target}")

    remarks = [f"REMARK 999 LPMOcat: unresolved residues rebuilt with MODELLER {version}\n",
               f"REMARK 999 LPMOcat: (best DOPE of {len(models)} models), occupancy 0.00:\n"]
    remarks += [f"REMARK 999 LPMOcat:   {plan['pdb_chain']} {g['start']}-{g['end']} {g['sequence']}\n"
                for g in plan["fill"]]
    out_path = os.path.join(out_dir, f"{struct_id}.pdb")
    with open(out_path, "w") as f:
        f.writelines(_splice(plan, model_residues, remarks))

    return {"source": os.path.basename(path), "chain": chain, "output": os.path.basename(out_path),
            "method": f"MODELLER {version} AutoModel, observed atoms fixed, "
                      f"best DOPE of {len(models)} models (seed {RAND_SEED})",
            "filled": plan["fill"], "skipped": plan["skip"],
            "n_term_missing": plan["n_term_missing"], "c_term_missing": plan["c_term_missing"],
            "dope": round(best["dope"], 1), **_check(plan, model_residues)}
