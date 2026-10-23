"""
Step 3 of the LPMOcat pipeline: preprocess the selected structures.

Every selected structure becomes one clean single-chain PDB file,
03_preparation/preprocessed/<id>.pdb, which US-align (step 4) and every per-structure
calculation after it use:

  3a  rebuild unresolved internal loops of crystal structures (MODELLER, fill_gaps.py);
      the rebuilt file replaces the download as input for 3b-3h
  3b  convert mmCIF to PDB format
  3c  keep the first model (NMR ensembles)
  3d  keep the selected chain
  3e  keep the first alternate conformer of each residue
  3f  rename modified residues to their standard parent (HIC -> HIS)
  3g  remove waters, ligands and all ions except the active-site metal
  3h  trim the N-terminus to the brace histidine His-1: a curated first residue, else the
      UniProt annotation, else the end of the SignalP 6 signal peptide (models), else the
      first histidine
  3i  remove a C-terminal expression tag: a curated last residue, else the residues the PDB
      entry annotates as an expression tag, cloning artifact or linker, else a C-terminal
      match of KNOWN_TAGS (His-tag, Strep-tag II, myc, ...)

3b-3i work on lists of PDB-format lines; preprocess_structure() runs them in order.
"""

import math
import os
import re
import shutil
import subprocess

import requests

# Standard amino acid residue names (3-letter codes)
STANDARD_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

# Same set populate_db recognises as the brace metal (3EJA has Mg in the brace)
BRACE_METALS = {"CU", "ZN", "NI", "CO", "MN", "FE", "MO", "MG"}

WATERS = ("HOH", "WAT", "DOD")

# His-1 is expected among the first residues of the structure
HIS_SEARCH_LIMIT = 35


# ---------------------------------------------------------------------------
# 3b: mmCIF -> PDB
# ---------------------------------------------------------------------------

def cif_to_pdb_lines(cif_path):
    """Convert mmCIF atom records to PDB-format lines.

    Extracts ATOM/HETATM records from mmCIF _atom_site loop and formats
    them as PDB fixed-width lines. Multi-character chain IDs are truncated
    to first character.
    """
    with open(cif_path) as f:
        lines = f.readlines()

    # Find _atom_site column order
    col_names = []
    in_atom_site = False
    data_start = None
    for i, line in enumerate(lines):
        if line.startswith("_atom_site."):
            in_atom_site = True
            col_names.append(line.strip().split(".")[1])
        elif in_atom_site and not line.startswith("_atom_site."):
            if not line.startswith("#") and line.strip():
                data_start = i
            break

    if not col_names or data_start is None:
        return []

    col_idx = {name: i for i, name in enumerate(col_names)}

    def get(fields, name, default=""):
        idx = col_idx.get(name)
        if idx is not None and idx < len(fields):
            v = fields[idx]
            return default if v == "?" or v == "." else v
        return default

    pdb_lines = []
    for i in range(data_start, len(lines)):
        line = lines[i].strip()
        if not line or line.startswith("#") or line.startswith("loop_"):
            break
        fields = line.split()
        record = get(fields, "group_PDB", "ATOM")
        if record not in ("ATOM", "HETATM"):
            continue

        serial = get(fields, "id", "1")
        atom_name = get(fields, "label_atom_id", "CA")
        alt_loc = get(fields, "label_alt_id", "")
        if alt_loc == ".":
            alt_loc = ""
        resname = get(fields, "label_comp_id", "UNK")
        chain = get(fields, "auth_asym_id",
                     get(fields, "label_asym_id", "A"))
        chain = chain[0] if chain else "A"  # truncate multi-char chain
        resseq = get(fields, "auth_seq_id",
                      get(fields, "label_seq_id", "1"))
        icode = get(fields, "pdbx_PDB_ins_code", "")
        if icode == "?":
            icode = ""
        x = get(fields, "Cartn_x", "0.0")
        y = get(fields, "Cartn_y", "0.0")
        z = get(fields, "Cartn_z", "0.0")
        occ = get(fields, "occupancy", "1.00")
        bfac = get(fields, "B_iso_or_equiv", "0.00")
        element = get(fields, "type_symbol", "")

        # Format atom name: 4 chars, left-padded if 1-3 chars
        if len(atom_name) < 4:
            atom_field = f" {atom_name:<3s}"
        else:
            atom_field = f"{atom_name:<4s}"

        pdb_line = (
            f"{record:<6s}"
            f"{int(serial):>5d} "
            f"{atom_field}"
            f"{alt_loc:>1s}"
            f"{resname:>3s} "
            f"{chain:>1s}"
            f"{int(resseq):>4d}"
            f"{icode:>1s}   "
            f"{float(x):>8.3f}"
            f"{float(y):>8.3f}"
            f"{float(z):>8.3f}"
            f"{float(occ):>6.2f}"
            f"{float(bfac):>6.2f}"
            f"          "
            f"{element:>2s}\n"
        )
        pdb_lines.append(pdb_line)

    pdb_lines.append("END\n")
    return pdb_lines


def read_structure(path):
    """3b: the structure as PDB-format lines, converting mmCIF files."""
    if path.endswith(".cif"):
        return cif_to_pdb_lines(path)
    with open(path) as f:
        return f.readlines()


# ---------------------------------------------------------------------------
# 3c-3e: model, chain, alternate conformers
# ---------------------------------------------------------------------------

def first_model(lines):
    """3c: the records before the first MODEL and those of model 1; everything after the
    first ENDMDL is dropped. Returns (lines, number of models)."""
    if not any(line.startswith("MODEL ") for line in lines):
        return lines, 0
    out, n_models, in_first = [], 0, False
    for line in lines:
        if line.startswith("MODEL "):
            n_models += 1
            if n_models == 1:
                in_first = True
            continue
        if line.startswith("ENDMDL"):
            in_first = False
            continue
        if in_first or n_models == 0:
            out.append(line)
    out.append("END\n")
    return out, n_models


def select_chain(lines, chain_id):
    """3d: drop the ATOM/HETATM records of all chains but chain_id, so that the alignment
    sees a single polypeptide. PDB format holds a one-character chain ID, so mmCIF chain IDs
    such as "AAA" are compared by their first character. Returns (lines, dropped chain IDs)."""
    if not chain_id:
        return lines, set()
    pdb_chain = chain_id[0]
    out, others = [], set()
    for line in lines:
        if line[:6].strip() in ("ATOM", "HETATM"):
            line_chain = line[21] if len(line) > 21 else " "
            if line_chain != pdb_chain:
                if line_chain.strip():
                    others.add(line_chain)
                continue
        out.append(line)
    return out, others


def keep_first_altloc(lines):
    """3e: resolve alternate conformations in PDB-format lines: per residue (chain, number,
    insertion code) keep the first conformer listed and clear the altloc indicator.

    The first conformer is normally A, but not always -- 6HAQ Arg58 was deposited as B/C/D
    only. Choosing per residue rather than per atom keeps every residue one consistent
    conformer, also under microheterogeneity (A and B being different amino acids).
    ANISOU/SIGATM/SIGUIJ records follow the atom they belong to.
    """
    kept, out, dropped = {}, [], False
    for line in lines:
        record = line[:6]
        if record in ("ATOM  ", "HETATM"):
            dropped = False
            alt_loc = line[16:17]
            if alt_loc.strip():
                if kept.setdefault(line[21:27], alt_loc) != alt_loc:
                    dropped = True
                    continue
                line = line[:16] + " " + line[17:]
        elif record in ("ANISOU", "SIGATM", "SIGUIJ"):
            if dropped:
                continue
            if line[16:17].strip():
                line = line[:16] + " " + line[17:]
        out.append(line)
    return out


# ---------------------------------------------------------------------------
# 3f-3g: modified residues, waters, ligands, ions
# ---------------------------------------------------------------------------

# Cache for parent residue lookups from RCSB Chemical Component Dictionary
_parent_residue_cache = {}


def lookup_parent_residue(resname):
    """Look up the standard parent amino acid for a modified residue
    using the RCSB Chemical Component Dictionary.

    Returns the 3-letter parent code (e.g. 'HIS' for 'HIC') or None
    if not found or not a standard amino acid modification.
    """
    if resname in _parent_residue_cache:
        return _parent_residue_cache[resname]

    try:
        r = requests.get(
            f"https://data.rcsb.org/rest/v1/core/chemcomp/{resname}",
            timeout=10)
        if r.status_code == 200:
            data = r.json()
            parents = (data.get("chem_comp", {})
                       .get("mon_nstd_parent_comp_id", []))
            if parents and parents[0] in STANDARD_RESIDUES:
                _parent_residue_cache[resname] = parents[0]
                return parents[0]
    except Exception:
        pass

    _parent_residue_cache[resname] = None
    return None


def _atom_residue_range(lines):
    """(first, last) residue number of the ATOM records."""
    numbers = set()
    for line in lines:
        if line[:4] == "ATOM":
            try:
                numbers.add(int(line[22:26]))
            except ValueError:
                pass
    return (min(numbers), max(numbers)) if numbers else (0, 0)


def _in_chain(line, residue_range):
    """A record numbered within, or right next to, the ATOM residues of the chain. HETATM
    residues like HIC at residue 1 are part of the chain even though they come before the
    first ATOM record."""
    try:
        return residue_range[0] - 1 <= int(line[22:26]) <= residue_range[1] + 1
    except ValueError:
        return False


def rename_modified_residues(lines):
    """3f: modified amino acids of the chain (e.g. HIC, the N-methylated His-1 of fungal
    LPMOs) become ATOM records of their standard parent from the RCSB Chemical Component
    Dictionary. Returns (lines, {modified name: (parent, number of atom lines)})."""
    residue_range = _atom_residue_range(lines)
    parents = {}
    for line in lines:
        if line[:6] == "HETATM" and _in_chain(line, residue_range):
            resname = line[17:20].strip()
            if resname not in STANDARD_RESIDUES and resname not in parents:
                parents[resname] = lookup_parent_residue(resname)

    out, renamed = [], {}
    for line in lines:
        if line[:6] == "HETATM" and _in_chain(line, residue_range):
            resname = line[17:20].strip()
            parent = parents.get(resname)
            if parent:
                line = "ATOM  " + line[6:17] + f"{parent:>3}" + line[20:]
                renamed[resname] = (parent, renamed.get(resname, (parent, 0))[1] + 1)
        out.append(line)
    return out, renamed


def active_site_metal(lines, max_dist=3.0):
    """The metal ion held by the histidine brace: a Cu (or other BRACE_METALS) HETATM with
    at least two histidine nitrogens within max_dist A -- surface ions bound by a single His
    (e.g. the extra Zn in 4OY6 and 5LSV) do not qualify. Returns its line, or None.
    If several qualify, the one with the most (then the closest) His nitrogens wins."""
    his_n, metals = [], []
    for line in lines:
        record = line[:6]
        if record not in ("ATOM  ", "HETATM"):
            continue
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        resname = line[17:20].strip()
        if resname in ("HIS", "HIC") and line[12:16].strip().startswith("N"):
            his_n.append(xyz)
        elif record == "HETATM":
            element = line[76:78].strip().upper() or resname.rstrip("0123456789")
            if element in BRACE_METALS:
                metals.append((line, xyz))
    best = None
    for line, xyz in metals:
        dists = sorted(math.dist(xyz, n) for n in his_n)
        n_bound = sum(d <= max_dist for d in dists)
        if n_bound >= 2 and (best is None or (-n_bound, dists[0]) < best[0]):
            best = ((-n_bound, dists[0]), line)
    return best[1] if best else None


def remove_non_protein(lines):
    """3g: keep the protein chain and the active-site metal only.

    Waters, ligands (glycans, buffer molecules) and other ions take no part in the alignment
    or in any later calculation; the brace metal stays for brace detection, protonation
    states and hole hopping. ANISOU records go too, and CONECT records survive only between
    remaining atoms (e.g. metal-His). Returns (lines, metal line or None,
    {removed residue name: number of atom lines})."""
    metal_line = active_site_metal(lines)
    residue_range = _atom_residue_range(lines)
    out, removed, kept_serials = [], {}, set()
    for line in lines:
        record = line[:6].strip()
        if record == "ANISOU":
            continue
        if record == "HETATM":
            resname = line[17:20].strip()
            if not ((_in_chain(line, residue_range) and resname in STANDARD_RESIDUES)
                    or line == metal_line):
                removed[resname] = removed.get(resname, 0) + 1
                continue
        if record in ("ATOM", "HETATM"):
            kept_serials.add(line[6:11].strip())
        elif record == "CONECT":
            body = line.rstrip("\n")
            serials = {body[k:k + 5].strip() for k in range(6, len(body), 5)} - {""}
            if not serials <= kept_serials:
                continue
        out.append(line)
    return out, metal_line, removed


# ---------------------------------------------------------------------------
# 3h: N-terminus
# ---------------------------------------------------------------------------

def _atom_residues(lines):
    """(residue number, residue name) of the ATOM residues in file order."""
    residues, seen = [], set()
    for line in lines:
        if line[:4] != "ATOM":
            continue
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        if resseq not in seen:
            seen.add(resseq)
            residues.append((resseq, line[17:20].strip()))
    return residues


def annotated_starts(annotation):
    """Candidates for the first residue of the mature protein from a domain annotation (UniProt
    numbering), best first: the UniProt chain start after a signal peptide, then the Pfam
    catalytic domain start. The chain start comes first because the mature LPMO begins at the
    brace histidine, while the Pfam match may begin a few residues later, on another histidine
    (G2RB72: chain 20 is the brace His, the PF03443 match starts at His 23, 10 A away from the
    internal His). The caller takes the first candidate that is a histidine.
    Returns [(residue, description)]."""
    out = []
    if not annotation:
        return out
    if annotation.get("signal_peptide") and annotation.get("chain"):
        out.append((annotation["chain"]["start"], "signal peptide (UniProt chain start)"))
    if annotation.get("catalytic_domain"):
        domain = annotation["catalytic_domain"]
        out.append((domain["start"], f"catalytic domain {domain['pfam_id']}"))
    return out


def annotated_start(annotation):
    """The first annotated start (annotated_starts), or (None, None)."""
    starts = annotated_starts(annotation)
    return starts[0] if starts else (None, None)


def _drop_before(lines, first_residue, keep_line=None):
    """Remove the ATOM/HETATM records numbered below first_residue (keep_line is kept)."""
    out = []
    for line in lines:
        if line[:6].strip() in ("ATOM", "HETATM") and line != keep_line:
            try:
                if int(line[22:26]) < first_residue:
                    continue
            except ValueError:
                pass
        out.append(line)
    return out


def trim_n_terminus(lines, annotation=None, metal_line=None, predicted_start=None, manual=None):
    """3h: the mature LPMO starts with the brace histidine His-1, so the residues before it
    (signal peptide, tags) are removed. The cut is, in order of precedence, at
      1. manual: a curated first residue (curation/n_terminus.csv), or "keep" for a chain
         that is right as deposited (6RS9 has no His-1)
      2. an annotated start (annotated_starts: the UniProt chain start, then the Pfam
         catalytic domain start) when that residue is a histidine
      3. predicted_start, the residue after the SignalP 6 signal peptide, when there is no
         annotated start and that residue is a histidine -- for models without annotation
      4. the first histidine among the first HIS_SEARCH_LIMIT residues -- a wrong annotation,
         a crystal structure numbered differently from UniProt, no annotation and no usable
         prediction; a structure without one is left untrimmed
    The active-site metal is never removed. Returns (lines, message or None)."""
    residues = _atom_residues(lines)
    names = dict(residues)

    if manual is not None:
        if str(manual).strip().lower() == "keep":
            return lines, "N-terminus kept as deposited (curation/n_terminus.csv)"
        try:
            first = int(manual)
        except ValueError:
            first = None
        if first not in names:
            return lines, (f"WARNING: curated first residue {manual!r} is not in the structure "
                           f"- not trimmed (curation/n_terminus.csv)")
        return (_drop_before(lines, first, metal_line),
                f"trimmed at residue {first} {names[first]} (curation/n_terminus.csv)")

    starts = annotated_starts(annotation)
    for start, source in starts:
        if names.get(start) in ("HIS", "HIC"):
            return (_drop_before(lines, start, metal_line),
                    f"trimmed signal peptide at residue {start} ({source})")
    start, source = starts[0] if starts else (None, None)

    if start is None and predicted_start is not None and names.get(predicted_start) in ("HIS", "HIC"):
        return (_drop_before(lines, predicted_start, metal_line),
                f"trimmed signal peptide at residue {predicted_start} (SignalP 6)")

    first_his = next((resseq for resseq, name in residues[:HIS_SEARCH_LIMIT]
                      if name in ("HIS", "HIC")), None)
    if start is not None:
        if first_his is None:
            return lines, (f"skipping trim — residue {start} is {names.get(start)}, not HIS/HIC, "
                           f"no HIS in first {HIS_SEARCH_LIMIT} residues ({source})")
        return (_drop_before(lines, first_his, metal_line),
                f"annotation start ({start}) is {names.get(start)} — trimmed to first HIS "
                f"at residue {first_his} ({source})")
    if first_his is not None and residues[0][1] not in ("HIS", "HIC"):
        return (_drop_before(lines, first_his, metal_line),
                f"trimmed to first HIS at residue {first_his} (no annotation)")
    return lines, None


THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# ---------------------------------------------------------------------------
# 3i: C-terminal expression tags
# ---------------------------------------------------------------------------

# Details a PDB entry gives for residues that are not part of the protein: SEQADV's conflict
# field (PDB format) and _struct_ref_seq_dif.details (mmCIF), both upper/lower case.
TAG_DETAILS = ("expression tag", "cloning artifact", "linker")

# Purification and epitope tags, and protease sites, as they are cloned onto the C-terminus.
# A tag is often only partly resolved, so a match is any prefix of the motif (the tag is read
# from the protein outwards) of at least MIN_TAG_MATCH residues -- two for the repeat of a
# poly-His tag, which is unambiguous, four for the rest.
KNOWN_TAGS = {
    "His-tag": ("HHHHHHHHHH",),
    "Strep-tag II": ("WSHPQFEK",),
    "myc": ("LEQKLISEEDL", "EQKLISEEDL"),      # LE: the XhoI site of the pET/pcDNA vectors
    "FLAG": ("DYKDDDDK",),
    "HA": ("YPYDVPDYA",),
    "V5": ("GKPIPNPLLGLDST",),
    "S-tag": ("KETAAAKFERQHMDS",),
    "AviTag": ("GLNDIFEAQKIEWHE",),
    "thrombin site": ("LVPRGS",),
    "TEV site": ("ENLYFQG",),
    "3C/PreScission site": ("LEVLFQGP",),
    "enterokinase site": ("DDDDK",),
    "Factor Xa site": ("IEGR",),
}
MIN_TAG_MATCH = {"His-tag": 3}
MIN_TAG_MATCH_DEFAULT = 4
# Cloning scars vectors leave between the protein and the tag; removed with the tag.
TAG_LINKERS = ("LE", "GS", "GG", "SA", "AAA", "GSG", "GGS", "SG", "GSA")
# Guard rails: a tag is short and the catalytic domain must survive it.
MAX_TAG_RESIDUES = 30
MIN_KEPT_RESIDUES = 100


def tagged_residues(path, chain_id=None):
    """3i: the author residue numbers a PDB entry annotates as an expression tag, cloning
    artifact or linker -- SEQADV records in PDB files, _struct_ref_seq_dif in mmCIF ones.

    These say which residues of the deposited chain are not part of the protein; they are the
    evidence for the C-terminal trim (an N-terminal one is 3h's business: 5O2X annotates its
    brace histidine as an expression tag). Returns a set of residue numbers, empty when the
    entry annotates none (or is a model, which has no such records)."""
    chain = chain_id[0] if chain_id else None
    out = set()
    if path.endswith(".cif"):
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
        cif = MMCIF2Dict(path)
        details = cif.get("_struct_ref_seq_dif.details") or []
        chains = cif.get("_struct_ref_seq_dif.pdbx_pdb_strand_id") or []
        numbers = cif.get("_struct_ref_seq_dif.pdbx_auth_seq_num") or []
        for detail, strand, number in zip(details, chains, numbers):
            if not any(d in str(detail).lower() for d in TAG_DETAILS):
                continue
            if chain and str(strand)[:1] != chain:
                continue
            try:
                out.add(int(number))
            except ValueError:
                pass
        return out
    with open(path) as f:
        for line in f:
            if not line.startswith("SEQADV"):
                continue
            if not any(d in line[49:].lower() for d in TAG_DETAILS):
                continue
            if chain and line[16:17] != chain:
                continue
            try:
                out.add(int(line[18:22]))
            except ValueError:
                pass
    return out


def _drop_after(lines, last_residue, keep_line=None):
    """Remove the ATOM/HETATM records numbered above last_residue (keep_line is kept)."""
    out = []
    for line in lines:
        if line[:6].strip() in ("ATOM", "HETATM") and line != keep_line:
            try:
                if int(line[22:26]) > last_residue:
                    continue
            except ValueError:
                pass
        out.append(line)
    return out


def _annotated_tag_start(residues, tagged):
    """The first residue of the trailing run of annotated tag residues, or None."""
    if not tagged:
        return None
    i = len(residues)
    while i > 0 and residues[i - 1][0] in tagged:
        i -= 1
    return i if i < len(residues) else None


def _known_tag_start(sequence):
    """(index of the first tag residue, tag name) for a C-terminal tag of KNOWN_TAGS read off
    the sequence, or (None, None). The resolved part of the tag is a prefix of the motif (the
    end of a tag is often disordered); a cloning scar right before it (TAG_LINKERS) is taken
    as part of the tag, and of two matches the longer one wins."""
    best, best_name = None, None
    for name, motifs in KNOWN_TAGS.items():
        least = MIN_TAG_MATCH.get(name, MIN_TAG_MATCH_DEFAULT)
        for motif in motifs:
            for length in range(min(len(motif), MAX_TAG_RESIDUES), least - 1, -1):
                start = len(sequence) - length
                if start <= 0 or sequence[start:] != motif[:length]:
                    continue
                for linker in sorted(TAG_LINKERS, key=len, reverse=True):
                    if sequence[start - len(linker):start] == linker:
                        start -= len(linker)
                        break
                if best is None or start < best:
                    best, best_name = start, name
                break
    return best, best_name


def trim_c_terminus(lines, tagged=(), manual=None, metal_line=None):
    """3i: remove a purification or epitope tag cloned onto the C-terminus, so that the
    catalogue holds the protein and not the construct. The tag is, in order of precedence,
      1. manual: a curated last residue (curation/c_terminus.csv), or "keep" for a chain that
         is right as deposited
      2. the trailing run of residues the PDB entry annotates as an expression tag, cloning
         artifact or linker (tagged_residues) -- the entry's own statement
      3. a C-terminal match of KNOWN_TAGS, for entries that annotate none (5NO7's myc tag)
    Only a trailing block is removed, at most MAX_TAG_RESIDUES residues and never so many
    that fewer than MIN_KEPT_RESIDUES are left. The active-site metal is never removed.

    Every structure goes through all three rules (user decision 2026-09-23). A model is
    folded from a sequence, and that sequence can carry a tag as easily as a construct can:
    TrEMBL entries are translations of EMBL/GenBank CDS records and CAZyme3D models a GenBank
    protein, and tagged expression constructs are deposited there. Only rule 2 is in practice
    PDB-only, since a model has no entry annotation to read.
    Returns (lines, message or None)."""
    residues = _atom_residues(lines)
    names = dict(residues)

    if manual is not None:
        if str(manual).strip().lower() == "keep":
            return lines, "C-terminus kept as deposited (curation/c_terminus.csv)"
        try:
            last = int(manual)
        except ValueError:
            last = None
        if last not in names:
            return lines, (f"WARNING: curated last residue {manual!r} is not in the structure "
                           f"- not trimmed (curation/c_terminus.csv)")
        return (_drop_after(lines, last, metal_line),
                f"trimmed after residue {last} {names[last]} (curation/c_terminus.csv)")

    start = _annotated_tag_start(residues, tagged)
    source, tag_name = "PDB expression-tag annotation", None
    if start is None:
        sequence = "".join(THREE_TO_ONE.get(name, "X") for _, name in residues)
        start, tag_name = _known_tag_start(sequence)
        source = f"{tag_name} sequence" if tag_name else None
    if start is None:
        return lines, None

    removed = residues[start:]
    if len(removed) > MAX_TAG_RESIDUES:
        return lines, (f"WARNING: {len(removed)} residues look like a C-terminal tag "
                       f"({source}) - more than {MAX_TAG_RESIDUES}, not trimmed")
    if start < MIN_KEPT_RESIDUES:
        return lines, (f"WARNING: trimming the C-terminal tag ({source}) would leave only "
                       f"{start} residues - not trimmed")
    sequence = "".join(THREE_TO_ONE.get(name, "X") for _, name in removed)
    return (_drop_after(lines, residues[start - 1][0], metal_line),
            f"removed a {len(removed)}-residue C-terminal tag, residues "
            f"{removed[0][0]}-{removed[-1][0]} ({sequence}; {source})")


SIGNALP6 = os.environ.get("SIGNALP6") or shutil.which("signalp6")


def structure_sequence(path, chain_id=None):
    """(residue numbers, one-letter sequence) of a structure's chain after 3b-3f -- the input
    for SignalP 6."""
    lines, _ = first_model(read_structure(path))
    lines, _ = select_chain(lines, chain_id)
    lines, _ = rename_modified_residues(keep_first_altloc(lines))
    residues = _atom_residues(lines)
    return [n for n, _ in residues], "".join(THREE_TO_ONE.get(name, "X") for _, name in residues)


def predict_signal_peptides(sequences, organism, out_dir):
    """3h: SignalP 6.0 (fast mode) predictions for {id: sequence}. organism: "eukarya" or
    "other" (bacteria, archaea). The SignalP output stays in out_dir.

    Returns ({id: {"prediction": e.g. "SP" / "LIPO" / "OTHER", "sp_end": last signal peptide
    position (1-based, None without one), "cleavage_probability"}}, error or None)."""
    if not SIGNALP6:
        return {}, "signalp6 not found (install SignalP 6, or set SIGNALP6 to its executable)"
    os.makedirs(out_dir, exist_ok=True)
    fasta = os.path.join(out_dir, "input.fasta")
    with open(fasta, "w") as f:
        f.writelines(f">{sid}\n{seq}\n" for sid, seq in sequences.items())
    try:
        result = subprocess.run([SIGNALP6, "--fastafile", fasta, "--output_dir", out_dir,
                                 "--format", "none", "--organism", organism, "--mode", "fast"],
                                capture_output=True, text=True, timeout=7200)
    except (OSError, subprocess.TimeoutExpired) as e:
        return {}, str(e)
    results_path = os.path.join(out_dir, "prediction_results.txt")
    if result.returncode != 0 or not os.path.exists(results_path):
        return {}, (result.stderr or result.stdout).strip()[-300:]

    predictions, header = {}, None
    with open(results_path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if line.startswith("# ID"):
                header = [h.lstrip("# ") for h in fields]
                continue
            if line.startswith("#") or not header:
                continue
            row = dict(zip(header, fields))
            site = re.search(r"CS pos: (\d+)-\d+\. Pr: ([\d.]+)", row.get("CS Position", ""))
            predictions[row["ID"]] = {"prediction": row["Prediction"],
                                      "sp_end": int(site.group(1)) if site else None,
                                      "cleavage_probability": float(site.group(2)) if site else None}
    return predictions, None


# ---------------------------------------------------------------------------
# 3b-3i for one structure
# ---------------------------------------------------------------------------

def preprocess_structure(src_path, out_path, chain_id=None, annotation=None,
                         predicted_start=None, manual=None, tagged=(), manual_c=None):
    """Run 3b-3i on one structure file and write the result to out_path (predicted_start and
    manual: see trim_n_terminus; tagged and manual_c: see trim_c_terminus).

    Returns (messages, metal) -- the log lines describing what changed and the residue name
    of the active-site metal kept (None if there is none)."""
    messages = []
    lines = read_structure(src_path)
    if src_path.endswith(".cif"):
        messages.append("converted mmCIF to PDB format")

    lines, n_models = first_model(lines)
    if n_models > 1:
        messages.append(f"extracted Model 1 from {n_models} NMR models")

    lines, other_chains = select_chain(lines, chain_id)
    if other_chains:
        messages.append(f"keeping chain {chain_id[0]}, dropped chains "
                        f"{', '.join(sorted(other_chains))}")

    lines = keep_first_altloc(lines)

    lines, renamed = rename_modified_residues(lines)
    for old, (parent, count) in renamed.items():
        messages.append(f"renamed {old} -> {parent} ({count} atom lines)")

    lines, metal_line, removed = remove_non_protein(lines)
    if removed:
        waters = sum(c for r, c in removed.items() if r in WATERS)
        others = ", ".join(f"{r}({c})" for r, c in sorted(removed.items()) if r not in WATERS)
        messages.append(f"dropped {waters} water atoms"
                        + (f"; ligands/ions: {others}" if others else ""))
    metal = metal_line[17:20].strip() if metal_line else None
    messages.append(f"active-site metal: {metal or 'none'}")

    lines, trim_message = trim_n_terminus(lines, annotation, metal_line,
                                          predicted_start=predicted_start, manual=manual)
    if trim_message:
        messages.append(trim_message)

    lines, tag_message = trim_c_terminus(lines, tagged, manual=manual_c, metal_line=metal_line)
    if tag_message:
        messages.append(tag_message)

    with open(out_path, "w") as f:
        f.writelines(lines)
    return messages, metal
