"""
LPMOcat data pipeline: from the CAZy lists of characterized LPMOs to aligned structures
with per-residue data.

Steps (--steps accepts names and numbers):

  1  scrape          1a  the CAZy characterized and structure pages
                     1b-1d  characterized LPMOs CAZy's lists miss (other_sources.py):
                         1b  Swiss-Prot entries with an LPMO EC number shown experimentally
                         1c  BRENDA proteins with a substrate row citing a paper
                         1d  curation/extra_proteins.csv (unreviewed entries with published
                             experimental results)
                         1e  curation/local_models.csv (proteins in no public database, with a
                             model file made for the catalogue)
                       -> 01_proteins/cazy_lpmos.csv, other_lpmos.csv
  2  fetch           One structure per entry (CAZy's first, then the other lists; a strain
                     variant of a CAZy entry only if that entry gets none; a local model of 1e
                     is taken as it is), the first found of
                       2a  the best PDB entry (RCSB)
                       2b  AlphaFold DB, by the UniProt IDs the entry lists
                       2c  AlphaFold DB, by UniProt entries with the sequence of CAZy's
                           GenBank protein (UniParc); also detects duplicate CAZy entries
                       2d  CAZyme3D, by GenBank ID
                     2e  drop structures with the same sequence as a better one
                     2f  UniProt accessions: CAZy's, the model's and the active entries with
                         the same sequence (UniParc), with status, EC numbers and papers
                     2g  InterPro/UniProt domain annotations of the accessions
                     2h  EC numbers: CAZy's, else reviewed UniProt's, else the family's most
                         common EC
                     2i  NCBI Taxonomy lineage of the source organisms
                     2j  substrates and oxidised positions with the papers showing them
                         (CAZy EC numbers, Swiss-Prot, BRENDA, curation/activities.csv)
                       -> 02_data/<id>.pdb|cif, pdb_selection.csv, uniprot_ids.json,
                          uniprot_references.json, domain_annotations.json, taxonomy.json,
                          activities.json
  3  prepare         3a  rebuild unresolved internal loops (<= 30 aa) of crystal structures
                         with MODELLER
                     3b  mmCIF -> PDB format          3c  first model (NMR)
                     3d  selected chain               3e  first alternate conformer
                     3f  modified residues -> standard parent (HIC -> HIS)
                     3g  remove waters, ligands and ions except the active-site metal
                     3h  trim the N-terminus to His-1: curated start
                         (curation/n_terminus.csv), else the annotation, else the SignalP 6
                         signal peptide end (models), else the first histidine
                     3i  remove a C-terminal expression tag: curated last residue
                         (curation/c_terminus.csv), else the residues the PDB entry annotates
                         as a tag, cloning artifact or linker, else a known tag motif
                       -> 03_preparation/filled/<id>.pdb + modelled_residues.json,
                          preprocessed/<id>.pdb, signal_peptides.json
  4  align           US-align -mm 4 multiple structural alignment
                       -> 04_alignment/alignment.fasta, matrices.json, aligned/<id>.pdb
  5  dssp_sasa       Secondary structure and solvent accessibility (sasa_dssp.py)
                     5a  SASA and relative SASA per residue (Shrake-Rupley, protein atoms only)
                     5b  secondary structure, sheets and bridge partners (DSSP, mkdssp)
                       -> 05_dssp_sasa/sasa.json, secondary_structure.json, dssp.json
  6  catalytic_domain
                     Catalytic domain end, disulfide bonds and histidine brace, through the
                     alignment
                     6a  disulfide bonds (Cys SG-SG <= 3.0 A)
                     6b  catalytic domain end: family consensus of annotated AlphaFold models
                     6c  cap C-terminal tails reaching far past the family consensus
                     6d  move an end lying in a long coil (linker) back to the last
                         helix/strand residue (the secondary structure of 5b)
                     6e  move an end lying inside a disulfide bond to its second Cys
                     6f  histidine brace: from the active-site metal, else from the alignment
                         columns of the brace in metal-bound structures; pseudo-Cu: the
                         centroid of the brace nitrogens, the copper position steps 7c, 8
                         and 12 use where there is no metal
                       -> 06_catalytic_domain/disulfides.json, catalytic_domain_boundaries.json,
                          histidine_brace.json
  7  core_fold       7a  beta-sandwich core: strands common to the families in the alignment
                         and the hairpin partner of the first, numbered b1-b9 (the same number
                         = the same strand), checked against the strand pairing (ladders)
                     7b  loop regions: the segments between given core strands
                     7c  binding surface: the residues a plane tangent to the copper face
                         touches (binding_surface.py)
                       -> 07_core_fold/beta_core.json, loop_regions.json, binding_surface.json
  8  second_sphere   8a  second-sphere aromatic: the Phe/Tyr closest to the metal, else at
                         the alignment columns of the metal-bound structures' aromatic
                     8b  hydrogen bonds, pi-pi stacking and cation-pi contacts of the brace
                     8c  second-sphere hydrogen-bonding motif (H-X8-Q-X-Y, E…N-X-F, ...)
                     8d  residue on the copper's solvent-facing axial position (AA10 Ala)
                     8e  sequence motifs after the brace histidines (N-His: H-G-X4-P-X2-R,
                         Internal-His: H-X-G-P, ...)
                       -> 08_second_sphere/brace_interactions.json
  9  distance_matrix pairwise distances between the catalytic domains (BLOSUM62 over the
                     aligned residues up to each catalytic domain end), for the dendrogram
                       -> 09_distance_matrix/distance_matrix.csv
  10 flexibility     10a B-factors (crystal structures) / pLDDT (models)
                     10b GNM (ProDy)
                     10c CABSflex RMSF (missing or outdated runs only with --cabsflex)
                       -> 10_flexibility/flexibility_raw.json, cabsflex_results/<id>/
  11 electrostatics  pdb2pqr (AMBER) with brace His-1 as HIE and His-2 as HID
                       -> 11_electrostatics/<id>.pqr
  12 hole_hopping    EHPath hole hopping paths to the Cu (6f's pseudo-Cu without a metal);
                     structures with unchanged inputs are reused
                       -> 12_hole_hopping/hole_hopping.json, <id>_paths.txt

Each step writes into the data folder carrying its number (logs/ holds the run logs) and
reads only the files of earlier steps, so steps can be rerun separately; a failing step
stops the run. Steps 5-8 annotate the aligned structures (the catalytic domain, the fold and
the active site), steps 9-12 compute their properties. Loading the results into the Django
catalogue (manage.py populate_db) is not part of this script.

Usage:
    python lpmocat_pipeline/pipeline.py                            # all steps
    python lpmocat_pipeline/pipeline.py --steps scrape,fetch       # by name
    python lpmocat_pipeline/pipeline.py --steps 5-12               # by number, ranges too
    python lpmocat_pipeline/pipeline.py --steps 10 --cabsflex      # run missing CABSflex
    python lpmocat_pipeline/pipeline.py --steps 1-2 --families AA15   # only these families
                                                                    # (overwrites the lists)
    python lpmocat_pipeline/pipeline.py --data-dir /path/to/run    # data elsewhere
"""

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import sys
from collections import Counter, defaultdict
from datetime import datetime
from itertools import combinations

import pandas as pd

# the step modules live in steps/ (imported by name inside the steps below)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "steps"))
from paths import CURATION_DIR, DEFAULT_DATA_DIR, PIPELINE_DIR  # noqa: E402

# 3h: curated N-termini for exceptions (struct_id, first_residue or "keep", reason); curated
# input, kept with the code rather than in a data directory
N_TERMINUS_CSV = os.path.join(CURATION_DIR, "n_terminus.csv")
C_TERMINUS_CSV = os.path.join(CURATION_DIR, "c_terminus.csv")

# 6a: A; S-S bonds are 2.03-2.05 A, radiation-damaged ones in crystal structures longer (4YN2
# Cys14-Cys33 and 7ZJB Cys64-Cys175 at 2.72 A, both annotated in the PDB); no other Cys pair of
# the catalogue comes closer than 3.5 A
DISULFIDE_MAX_DISTANCE = 3.0
# 6c: alignment columns an end may lie past its family consensus before it is capped
CAP_MARGIN = 100
# 6d: a coil run at least this long (residues) is taken as an interdomain linker
LINKER_MIN_COIL = 60
# 6f: His-2 of a metal-free structure is searched at the His-2 columns of metal-bound
# structures and up to this many columns next to them
BRACE_COLUMN_WINDOW = 3
# 6f: braces whose His-1 and His-2 nitrogens are further apart (A) are reported
BRACE_WARN_N_DIST = 7.0
# 6f: the nitrogens that hold the copper, (brace histidine, atom); their centroid is the pseudo-Cu
BRACE_NITROGENS = (("his1", "N"), ("his1", "ND1"), ("his1", "NE2"),
                   ("his2", "ND1"), ("his2", "NE2"))
# 6f: the structure's metal is the copper position when a ring nitrogen (ND1 or NE2) of each brace
# histidine lies within this distance (A); otherwise the pseudo-Cu. The amine is not required: in
# 7VFC the N-His ring is flipped and the amine 5.9 A away, yet both rings still hold the Cu
COPPER_MAX_RING_N_DISTANCE = 3.0

# Structure sources whose B-factor column holds pLDDT
MODEL_SOURCES = ("AlphaFold", "CAZyme3D", "Local model")
# 5a: relative SASA (%) from which a residue counts as surface exposed (only reported here;
# the app uses the same threshold)
SURFACE_RSA = 25
# 10c: fingerprint of the structure a CABSflex run simulated, kept in its run directory
CABSFLEX_INPUT_KEY = "input.sha1"


def set_data_dir(data_dir):
    """Point every input and output path at data_dir (default: the script's directory).
    Each step writes into the folder carrying its number."""
    global DATA_DIR, CAZY_CSV, OTHER_CSV, STRUCTURES_DIR, PDB_SELECTION_CSV, UNIPROT_IDS_JSON, \
        UNIPROT_REFERENCES_JSON, DOMAIN_ANNOTATIONS_JSON, TAXONOMY_JSON, ACTIVITIES_JSON, \
        PREPROCESSING_DIR, FILLED_DIR, \
        SIGNAL_PEPTIDES_JSON, PREPROCESSED_DIR, ALIGNMENT_DIR, ALIGNED_DIR, ALIGNMENT_FASTA, \
        ANNOTATION_DIR, BOUNDARIES_JSON, BRACE_JSON, DISTANCE_MATRIX_CSV, FLEXIBILITY_JSON, \
        CABSFLEX_DIR, PQR_DIR, SASA_JSON, SECONDARY_STRUCTURE_JSON, DSSP_JSON, DISULFIDES_JSON, \
        BRACE_INTERACTIONS_JSON, BETA_CORE_JSON, LOOPS_JSON, BINDING_SURFACE_JSON, \
        HOLE_HOPPING_DIR, HOLE_HOPPING_JSON, LOG_DIR
    DATA_DIR = data_dir
    CAZY_CSV = os.path.join(data_dir, "01_proteins", "cazy_lpmos.csv")
    OTHER_CSV = os.path.join(data_dir, "01_proteins", "other_lpmos.csv")
    STRUCTURES_DIR = os.path.join(data_dir, "02_data")
    PDB_SELECTION_CSV = os.path.join(STRUCTURES_DIR, "pdb_selection.csv")
    UNIPROT_IDS_JSON = os.path.join(STRUCTURES_DIR, "uniprot_ids.json")
    UNIPROT_REFERENCES_JSON = os.path.join(STRUCTURES_DIR, "uniprot_references.json")
    DOMAIN_ANNOTATIONS_JSON = os.path.join(STRUCTURES_DIR, "domain_annotations.json")
    TAXONOMY_JSON = os.path.join(STRUCTURES_DIR, "taxonomy.json")
    ACTIVITIES_JSON = os.path.join(STRUCTURES_DIR, "activities.json")
    PREPROCESSING_DIR = os.path.join(data_dir, "03_preparation")
    FILLED_DIR = os.path.join(PREPROCESSING_DIR, "filled")
    SIGNAL_PEPTIDES_JSON = os.path.join(PREPROCESSING_DIR, "signal_peptides.json")
    PREPROCESSED_DIR = os.path.join(PREPROCESSING_DIR, "preprocessed")
    ALIGNMENT_DIR = os.path.join(data_dir, "04_alignment")
    ALIGNED_DIR = os.path.join(ALIGNMENT_DIR, "aligned")
    ALIGNMENT_FASTA = os.path.join(ALIGNMENT_DIR, "alignment.fasta")
    ANNOTATION_DIR = os.path.join(data_dir, "06_catalytic_domain")
    BOUNDARIES_JSON = os.path.join(ANNOTATION_DIR, "catalytic_domain_boundaries.json")
    BRACE_JSON = os.path.join(ANNOTATION_DIR, "histidine_brace.json")
    DISULFIDES_JSON = os.path.join(ANNOTATION_DIR, "disulfides.json")
    DISTANCE_MATRIX_CSV = os.path.join(data_dir, "09_distance_matrix", "distance_matrix.csv")
    FLEXIBILITY_JSON = os.path.join(data_dir, "10_flexibility", "flexibility_raw.json")
    CABSFLEX_DIR = os.path.join(data_dir, "10_flexibility", "cabsflex_results")
    PQR_DIR = os.path.join(data_dir, "11_electrostatics")
    SASA_JSON = os.path.join(data_dir, "05_dssp_sasa", "sasa.json")
    SECONDARY_STRUCTURE_JSON = os.path.join(data_dir, "05_dssp_sasa", "secondary_structure.json")
    DSSP_JSON = os.path.join(data_dir, "05_dssp_sasa", "dssp.json")
    BRACE_INTERACTIONS_JSON = os.path.join(data_dir, "08_second_sphere",
                                           "brace_interactions.json")
    BETA_CORE_JSON = os.path.join(data_dir, "07_core_fold", "beta_core.json")
    LOOPS_JSON = os.path.join(data_dir, "07_core_fold", "loop_regions.json")
    BINDING_SURFACE_JSON = os.path.join(data_dir, "07_core_fold", "binding_surface.json")
    HOLE_HOPPING_DIR = os.path.join(data_dir, "12_hole_hopping")
    HOLE_HOPPING_JSON = os.path.join(HOLE_HOPPING_DIR, "hole_hopping.json")
    LOG_DIR = os.path.join(data_dir, "logs")


set_data_dir(DEFAULT_DATA_DIR)


class StepError(Exception):
    """A step cannot run or failed; the pipeline stops."""


class TeeStream:
    """Write to a terminal stream and a log file."""

    def __init__(self, log_file, original_stream):
        self.log_file = log_file
        self.original = original_stream

    def write(self, text):
        self.original.write(text)
        self.log_file.write(text)
        self.log_file.flush()

    def flush(self):
        self.original.flush()
        self.log_file.flush()


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _listed_entries():
    """The entries of step 1: cazy_lpmos.csv, then other_lpmos.csv (when there), with
    list_source ("CAZy", "Swiss-Prot", "BRENDA" or "curated") and variant_of."""
    frames = [pd.read_csv(CAZY_CSV).assign(list_source="CAZy", variant_of="")]
    if os.path.exists(OTHER_CSV):
        frames.append(pd.read_csv(OTHER_CSV))
    return pd.concat(frames, ignore_index=True)


def _selected_structures():
    """Rows of pdb_selection.csv that have a structure."""
    if not os.path.exists(PDB_SELECTION_CSV):
        raise StepError(f"{PDB_SELECTION_CSV} not found. Run the 'fetch' step first.")
    sel_df = pd.read_csv(PDB_SELECTION_CSV)
    return [row for _, row in sel_df.iterrows() if row.source != "NONE"]


def _label(row):
    return f"{row.family}_{row.organism}_{row.protein_name}"[:60]


def _chain(row):
    return str(row.selected_chain).strip() if pd.notna(row.selected_chain) else "A"


def _primary_uniprot(row):
    """The UniProt accession whose numbering and annotations belong to a structure: the
    model's own for AlphaFold DB models, the first one CAZy lists for PDB entries, and ''
    otherwise (CAZyme3D models, entries without UniProt ID). row: dict or pandas row."""
    struct_id = str(row.get("selected_pdb") or "")
    if row.get("source") == "AlphaFold" and struct_id.startswith("AF-"):
        return struct_id[3:]
    uids = row.get("uniprot_ids")
    if row.get("source") != "PDB" or not isinstance(uids, str):
        return ""
    return uids.split(";")[0].strip()


def _split_ids(value):
    """The accessions of a ';'-separated field (a dict value or pandas cell; NaN -> [])."""
    if not isinstance(value, str):
        return []
    return [v.strip() for v in value.split(";") if v.strip() and v.strip() != "nan"]


def _downloaded_structure(struct_id):
    """Path of the downloaded structure file (.pdb or .cif), or None."""
    for ext in (".pdb", ".cif"):
        path = os.path.join(STRUCTURES_DIR, f"{struct_id}{ext}")
        if os.path.exists(path):
            return path
    return None


def _remove_stale(directory, keep, suffix):
    """Remove the files ending in suffix whose names are not in keep (outputs of
    structures that are no longer selected, or failed this time)."""
    for name in sorted(os.listdir(directory)):
        if name.endswith(suffix) and name not in keep:
            os.remove(os.path.join(directory, name))
            print(f"  removed outdated {name}")


def _require_mkdssp():
    """Stop unless mkdssp (DSSP, step 5b) is installed."""
    from sasa_dssp import mkdssp_version

    version = mkdssp_version()
    if not version:
        raise StepError("mkdssp not found or not working: install it with 'conda create "
                        "-n dssp -c conda-forge dssp=4.4.11' or set MKDSSP to the binary")
    print(f"  DSSP: {version}")


def _read_alignment():
    """{struct_id: aligned sequence} from alignment.fasta (headers like '7OKR.pdb:A')."""
    from Bio import SeqIO

    if not os.path.exists(ALIGNMENT_FASTA):
        raise StepError(f"{ALIGNMENT_FASTA} not found. Run the 'align' step first.")
    return {record.id.split(".pdb")[0]: str(record.seq)
            for record in SeqIO.parse(ALIGNMENT_FASTA, "fasta")}


def _residues(lines):
    """ATOM residues in file order -- the order of the residues in the alignment:
    [(number, name, {atom name: (x, y, z)})], one per chain, number and insertion code."""
    residues, index = [], {}
    for line in lines:
        if not line.startswith("ATOM"):
            continue
        try:
            number = int(line[22:26])
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        key = (line[21], line[22:27])
        if key not in index:
            index[key] = len(residues)
            residues.append((number, line[17:20].strip(), {}))
        residues[index[key]][2][line[12:16].strip()] = xyz
    return residues


def _aligned_lines(struct_ids):
    """{struct_id: lines of 04_alignment/aligned/<id>.pdb}."""
    lines = {}
    for struct_id in struct_ids:
        path = os.path.join(ALIGNED_DIR, f"{struct_id}.pdb")
        if not os.path.exists(path):
            raise StepError(f"{path} not found. Run the 'align' step first.")
        with open(path) as f:
            lines[struct_id] = f.readlines()
    return lines


def _load_json(path, step_name):
    """Contents of a JSON file written by an earlier step."""
    if not os.path.exists(path):
        raise StepError(f"{path} not found. Run the '{step_name}' step first.")
    with open(path) as f:
        return json.load(f)


def _alignment_columns(seq, residues):
    """0-based alignment column of every residue, or None when the structure's residue count
    differs from its aligned sequence."""
    cols = [i for i, c in enumerate(seq) if c != "-"]
    return cols if len(cols) == len(residues) else None


# ---------------------------------------------------------------------------
# Step 1: List the LPMOs (CAZy and the other sources)
# ---------------------------------------------------------------------------

def step_scrape(families=None):
    """1a: scrape the CAZy characterized & structure pages -> cazy_lpmos.csv; 1b-1d: the
    characterized LPMOs of Swiss-Prot, BRENDA and curation/extra_proteins.csv that CAZy's lists
    miss -> other_lpmos.csv."""
    import scrape_cazy
    from activities import find_brenda_file
    from other_sources import collect as collect_other_lpmos

    if families:
        # Temporarily override the family list
        original = scrape_cazy.LPMO_FAMILIES
        scrape_cazy.LPMO_FAMILIES = [f for f in families if f in original]
        if not scrape_cazy.LPMO_FAMILIES:
            scrape_cazy.LPMO_FAMILIES = original
            raise StepError(f"no valid families in {families}")
        try:
            proteins = scrape_cazy.scrape_all_lpmo_families()
        finally:
            scrape_cazy.LPMO_FAMILIES = original
    else:
        proteins = scrape_cazy.scrape_all_lpmo_families()

    df = pd.DataFrame(proteins)
    df = df.sort_values(
        ["family", "organism", "protein_name"]).reset_index(drop=True)
    os.makedirs(os.path.dirname(CAZY_CSV), exist_ok=True)
    df.to_csv(CAZY_CSV, index=False)

    char_yes = df.in_cazy_characterized == "yes"
    struct_yes = df.in_cazy_structure == "yes"

    print(f"\nSaved {len(proteins)} proteins to {CAZY_CSV}")
    print(f"  Characterized: {char_yes.sum()}, Structure: {struct_yes.sum()}, "
          f"Both: {(char_yes & struct_yes).sum()}")

    print("\n1b-1d: characterized LPMOs outside CAZy's lists (Swiss-Prot, BRENDA, curated)...")
    other = collect_other_lpmos(pd.read_csv(CAZY_CSV, dtype=str), find_brenda_file(DATA_DIR))
    if families:
        other = other[other.family.isin(families)]
    other.to_csv(OTHER_CSV, index=False)
    print(f"Saved {len(other)} entries to {OTHER_CSV} ("
          + ", ".join(f"{source} {n}" for source, n in other.list_source.value_counts().items()) + ")")


# ---------------------------------------------------------------------------
# Step 2: Fetch data: structures, annotations, taxonomy and activities
# ---------------------------------------------------------------------------

_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "HIC": "H", "SEC": "U",
}

_SOURCE_RANK = {"PDB": 0, "AlphaFold": 1, "CAZyme3D": 2, "Local model": 3, "NONE": 4}


def _seq_from_pdb(pdb_path, chain=None):
    """Extract amino acid sequence from ATOM records, trimmed to the first HIS.

    LPMOs always start at a conserved N-terminal histidine after signal-peptide
    cleavage. Trimming to the first HIS normalises structures where the AF model
    still includes the signal peptide but the crystal structure does not, so that
    genuinely identical proteins compare equal regardless of source.

    If chain is given (e.g. 'A'), only that chain is read.
    Returns a sequence string, or '' on error.
    """
    seen = {}   # (chain_id, res_num) → res_name; first occurrence wins
    try:
        with open(pdb_path) as fh:
            for line in fh:
                if not line.startswith("ATOM"):
                    continue
                try:
                    res_name = line[17:20].strip()
                    ch = line[21]
                    res_num = int(line[22:26])
                except (ValueError, IndexError):
                    continue
                if chain and ch != chain:
                    continue
                if res_name in _THREE_TO_ONE and (ch, res_num) not in seen:
                    seen[(ch, res_num)] = res_name
    except OSError:
        return ""
    full = [(k, _THREE_TO_ONE[v]) for k, v in sorted(seen.items())]
    # Trim to first HIS (N-terminal brace histidine)
    for i, (_, aa) in enumerate(full):
        if aa == "H":
            return "".join(aa for _, aa in full[i:])
    return "".join(aa for _, aa in full)  # no HIS found — return as-is


def _dedup_by_sequence(selection_rows, structure_files, characterized=frozenset()):
    """2e: drop structures whose sequence is identical to a higher-priority one.

    Priority: PDB (any method) > AlphaFold > CAZyme3D. Among the same source, the lower
    resolution; then the entry CAZy lists as characterized (characterized: a set of
    (family, protein_name, organism)) -- e.g. 6TWE over 5LW4, two NMR structures of
    B. licheniformis LPMO10A without a resolution; then an entry of CAZy's lists over one of
    the other lists; then a UniProt ID listed by CAZy over one found through UniParc; then the
    first encountered.

    Returns (deduplicated_selection_rows, deduplicated_structure_files).
    """
    # Build a parallel list pairing each structure_file with its selection_row.
    # selection_rows includes NONE entries; structure_files does not.
    sf_iter = iter(structure_files)
    paired = []   # (sel_row, sf_tuple_or_None)
    for row in selection_rows:
        if row.get("source", "NONE") != "NONE":
            paired.append((row, next(sf_iter)))
        else:
            paired.append((row, None))

    # Extract sequences for non-NONE rows
    seq_map = {}   # index → sequence string
    for i, (row, sf) in enumerate(paired):
        if sf is None:
            continue
        pdb_path = sf[1]
        chain = row.get("selected_chain", "A") or "A"
        # chain field can be multi-char (e.g. "AAA" from RCSB); use first char
        chain = chain[0] if chain else "A"
        seq_map[i] = _seq_from_pdb(pdb_path, chain)

    # Group indices by sequence, then pick the best per group
    groups = defaultdict(list)
    for idx, seq in seq_map.items():
        if seq:
            groups[seq].append(idx)

    to_drop = set()
    for seq, idxs in groups.items():
        if len(idxs) < 2:
            continue
        # Sort by priority: source rank, then resolution (lower = better)
        def _key(i):
            row = paired[i][0]
            src_rank = _SOURCE_RANK.get(row.get("source", "NONE"), 9)
            try:
                res = float(row.get("resolution") or 99)
            except ValueError:
                res = 99
            entry = (row.get("family"), row.get("protein_name"), row.get("organism"))
            # a characterized CAZy entry beats a structure-only one, and a UniProt ID listed
            # by CAZy beats one found through UniParc
            return (src_rank, res, row.get("list_source", "CAZy") != "CAZy",
                    entry not in characterized, bool(row.get("uniprot_source")))
        idxs_sorted = sorted(idxs, key=_key)
        winner = idxs_sorted[0]
        for loser in idxs_sorted[1:]:
            to_drop.add(loser)
            row_l = paired[loser][0]
            row_w = paired[winner][0]
            print(f"  DEDUP drop {row_l['selected_pdb']!r} "
                  f"({row_l['source']}) — identical sequence to "
                  f"{row_w['selected_pdb']!r} ({row_w['source']})")

    new_sel_rows = [row for i, (row, _) in enumerate(paired) if i not in to_drop]
    new_struct_files = [sf for i, (_, sf) in enumerate(paired)
                        if sf is not None and i not in to_drop]
    return new_sel_rows, new_struct_files


def collect_uniprot_accessions(selection_rows, listed_df):
    """2f: every UniProt accession of the selected proteins (see
    fetch_structures.uniprot_identity): those of the entry (source: its list, "CAZy",
    "Swiss-Prot", "BRENDA" or "curated"), those of the other lists' strain variants dropped for
    it in step 2 ("strain variant"), the AlphaFold DB model's
    (found through UniParc when CAZy does not list it), the UniProt entry of a CAZyme3D model,
    and every active entry with the same sequence. Sets "uniprot_ids" (all), "uniprot_active"
    (not deleted) and "uniparc_ids" of the rows (dicts, modified in place).

    Returns (identity, references): {struct_id: {"checked", "uniparc", "accessions":
    {accession: {"source", "status", "organism", "ec_numbers"}}}} and {struct_id: [the papers
    the active entries cite, one per PubMed ID/DOI, with the citing "accessions"]}."""
    from fetch_structures import uniprot_identity

    listed_rows = {(r.family, r.protein_name, r.organism): r for r in listed_df.itertuples()}
    checked = datetime.now().strftime("%Y-%m-%d")
    identity, references = {}, {}
    for row in selection_rows:
        if row.get("source") == "NONE":
            continue
        entry = listed_rows.get((row["family"], row["protein_name"], row["organism"]))
        cazy_ids = _split_ids(entry.uniprot_ids) if entry is not None else []
        genbank_ids = _split_ids(entry.genbank_ids) if entry is not None else []
        list_source = entry.list_source if entry is not None else "CAZy"
        listed = [(acc, list_source) for acc in cazy_ids]
        if row["source"] == "AlphaFold" and row["selected_pdb"][3:] not in cazy_ids:
            listed.insert(0, (row["selected_pdb"][3:], "UniParc"))
        elif row["source"] == "CAZyme3D":
            model = re.match(r"AlphaFold DB model of (\S+)", str(row.get("method") or ""))
            if model and model.group(1) not in cazy_ids:
                listed.append((model.group(1), "CAZyme3D model"))

        info = uniprot_identity(listed, genbank_ids, [(acc, "strain variant")
                                                      for acc in _split_ids(row.get("strain_variants"))])
        papers = {}
        for acc, entry in info["accessions"].items():
            for ref in entry.pop("references"):
                key = ref["pmid"] or ref["doi"].lower()
                papers.setdefault(key, {**ref, "accessions": []})["accessions"].append(acc)
        active = [acc for acc, entry in info["accessions"].items() if entry["status"] != "deleted"]
        row["uniprot_ids"] = "; ".join(info["accessions"])
        row["uniprot_active"] = "; ".join(active)
        row["uniparc_ids"] = "; ".join(info["uniparc"])
        identity[row["selected_pdb"]] = {"checked": checked, **info}
        if papers:
            references[row["selected_pdb"]] = list(papers.values())
    return identity, references


def add_ec_numbers(selection_rows, cazy_df, identity=None, other_df=None):
    """2h: EC numbers of every selection row (dicts, modified in place): those of its entry
    (CAZy's, or those the other list shows experimentally -- other_df, other_lpmos.csv); else
    those of its reviewed (Swiss-Prot) UniProt entries (identity from 2f); else the family
    default -- the most common EC number string among the family's CAZy entries that list one
    (CAZy gives none for structure-only entries). Sets "ec_numbers" and "ec_source" ("CAZy",
    "Swiss-Prot", "BRENDA", "UniProt <accessions>", "family default" or ""). Returns the
    family defaults."""
    listed = {}
    for source, frame in (("CAZy", cazy_df), (None, other_df)):
        for _, r in (frame.iterrows() if frame is not None else ()):
            ec = str(r.ec_numbers).strip() if pd.notna(r.ec_numbers) else ""
            listed[(r.family, r.protein_name, r.organism)] = (ec, source or r.list_source)
    defaults = {}
    for family, group in cazy_df.groupby("family"):
        ecs = [str(v).strip() for v in group.ec_numbers if pd.notna(v) and str(v).strip()]
        if ecs:
            defaults[family] = Counter(ecs).most_common(1)[0][0]

    for row in selection_rows:
        ec, source = listed.get((row["family"], row["protein_name"], row["organism"]), ("", ""))
        accessions = (identity or {}).get(row.get("selected_pdb"), {}).get("accessions", {})
        reviewed = {acc: entry["ec_numbers"] for acc, entry in accessions.items()
                    if entry["status"] == "reviewed" and entry["ec_numbers"]}
        if ec:
            row["ec_numbers"], row["ec_source"] = ec, source
        elif reviewed:
            row["ec_numbers"] = "; ".join(sorted({e for ecs in reviewed.values() for e in ecs}))
            row["ec_source"] = "UniProt " + ", ".join(reviewed)
        elif row["family"] in defaults:
            row["ec_numbers"], row["ec_source"] = defaults[row["family"]], "family default"
        else:
            row["ec_numbers"], row["ec_source"] = "", ""
    return defaults


def step_fetch(families=None):
    """2a-2j: select and download one structure per entry of step 1 (CAZy's lists, then the
    other lists), then collect the UniProt accessions, domain annotations, EC numbers,
    taxonomy and activities of the selected proteins."""
    from fetch_structures import (
        fetch_pdb_metadata, select_best_pdb, download_pdb,
        download_alphafold, check_alphafold_available,
        uniparc_uniprot_ids, download_cazyme3d, cazyme3d_struct_id, import_local_model,
        fetch_all_domain_annotations, pinned_structures,
    )
    # PDB entries chosen by hand (curation/structures.csv), which win in 2a
    pinned = pinned_structures()

    if not os.path.exists(CAZY_CSV):
        raise StepError(f"{CAZY_CSV} not found. Run the 'scrape' step first.")

    cazy_df = pd.read_csv(CAZY_CSV)
    if not os.path.exists(OTHER_CSV):
        print(f"No {OTHER_CSV}: only CAZy's entries (rerun step 1, or lpmocat_pipeline/steps/other_sources.py)")
    df = listed_df = _listed_entries()
    other_df = df[df.list_source != "CAZy"]
    # the CAZy entry of each UniProt accession and GenBank protein CAZy lists, for the strain
    # variants of the other lists
    cazy_entry_of = {}
    for r in cazy_df.itertuples():
        for db_id in _split_ids(r.uniprot_ids) + _split_ids(r.genbank_ids):
            cazy_entry_of.setdefault(db_id, (r.family, r.protein_name, r.organism))
    # UniProt accession -> CAZy entry listing it. A UniParc match to one of these means an
    # identical-sequence duplicate of that entry (e.g. AfAA9B of A. fumigatus var. niveus)
    cazy_uniprot_owner = {}
    for _, r in cazy_df.iterrows():
        for u in (str(r.uniprot_ids).split(";") if pd.notna(r.uniprot_ids) else []):
            if u.strip():
                cazy_uniprot_owner.setdefault(u.strip(), f"{r.family} {r.protein_name}")
    if families:
        df = df[df.family.isin(families)]
    print(f"Processing {len(df)} proteins (" + ", ".join(
        f"{source} {n}" for source, n in df.list_source.value_counts().items()) + ")")

    os.makedirs(STRUCTURES_DIR, exist_ok=True)

    # Collect all PDB codes and fetch metadata
    all_pdb_codes = set()
    for _, row in df.iterrows():
        if pd.notna(row.pdb_codes) and str(row.pdb_codes).strip():
            for code in str(row.pdb_codes).split(","):
                code = code.strip().upper()
                if code:
                    all_pdb_codes.add(code)

    pdb_metadata = {}
    if all_pdb_codes:
        print(f"Fetching metadata for {len(all_pdb_codes)} PDB entries...")
        pdb_metadata = fetch_pdb_metadata(list(all_pdb_codes))

    # Select best PDB and download structures
    selection_rows = []
    structure_files = []

    for _, row in df.iterrows():
        protein_name = row.protein_name
        family = row.family
        organism = row.organism
        pdb_codes = (str(row.pdb_codes).strip()
                     if pd.notna(row.pdb_codes) else "")
        uniprot_ids = (str(row.uniprot_ids).strip()
                       if pd.notna(row.uniprot_ids) else "")

        label = f"{family}_{organism}_{protein_name}"[:60]

        # 1e: a model file made for the catalogue (the protein is in no public database)
        local_file = row.get("local_file")
        if isinstance(local_file, str) and local_file:
            filepath = import_local_model(os.path.join(PIPELINE_DIR, local_file), row.struct_id, STRUCTURES_DIR)
            structure_files.append((label, filepath, "", "A"))
            print(f"  LOC  {row.struct_id} ({local_file}) - {protein_name[:40]}")
            selection_rows.append({
                "family": family, "protein_name": protein_name, "organism": organism,
                "source": "Local model", "selected_pdb": row.struct_id, "selected_chain": "A",
                # 1e may name the UniProt entry of the protein the model stands in for, so the
                # entry keeps its EC numbers, papers and activities (2f-2j)
                "all_pdb_codes": "", "uniprot_ids": uniprot_ids, "all_af_models": "", "pdb_code": "",
                "chain_id": "A", "resolution": "",
                # AlphaFold Server names its files fold_<job>_model_<rank>.cif
                "method": ("AlphaFold 3 (AlphaFold Server)"
                           if os.path.basename(local_file).startswith("fold_") else "local model"),
                "seq_length": "", "modeled_count": "", "coverage": 1.0,
            })
            continue

        # a strain variant from another list stands in only for a CAZy entry without a structure
        if pd.notna(row.variant_of) and row.variant_of:
            owner = cazy_entry_of.get(row.variant_of)
            owner_row = next((r for r in selection_rows
                              if (r["family"], r["protein_name"], r["organism"]) == owner
                              and (r["source"] != "NONE" or r.get("duplicate_of"))), None)
            if owner_row:
                print(f"  VAR  {uniprot_ids} - {protein_name[:40]}: strain variant of "
                      f"{owner[1][:40]} ({row.variant_of}), which has a structure")
                # its UniProt entry joins the CAZy entry's accessions (2f), with its papers
                if owner_row["source"] != "NONE":
                    owner_row["strain_variants"] = "; ".join(
                        _split_ids(owner_row.get("strain_variants")) + _split_ids(uniprot_ids))
                selection_rows.append({
                    "family": family, "protein_name": protein_name, "organism": organism,
                    "source": "NONE", "selected_pdb": "", "selected_chain": "",
                    "all_pdb_codes": pdb_codes, "uniprot_ids": uniprot_ids, "all_af_models": "",
                    "duplicate_of": f"{owner[0]} {owner[1]}",
                })
                continue

        all_uids = []
        if uniprot_ids and uniprot_ids != "nan":
            all_uids = [u.strip() for u in uniprot_ids.split(";")
                        if u.strip()]

        # 2a: best PDB entry
        pdb_selected = False
        if pdb_codes and pdb_codes != "nan":
            pdb_code, chain_id, sel_info = select_best_pdb(
                pdb_codes, pdb_metadata, pinned=pinned)
            if pdb_code in pinned:
                print(f"  {pdb_code} chosen by curation/structures.csv")

            if pdb_code:
                filepath = download_pdb(pdb_code, STRUCTURES_DIR)
                if filepath:
                    af_models = [u for u in all_uids
                                 if check_alphafold_available(u)]
                    print(f"  PDB  {pdb_code} chain {chain_id}"
                          f" ({sel_info['resolution']}A,"
                          f" {sel_info['method']},"
                          f" coverage={sel_info['coverage']:.0%})"
                          f" - {protein_name[:40]}")
                    primary_uid = all_uids[0] if all_uids else ""
                    structure_files.append((label, filepath, primary_uid,
                                            chain_id))
                    selection_rows.append({
                        "family": family,
                        "protein_name": protein_name,
                        "organism": organism,
                        "source": "PDB",
                        "selected_pdb": pdb_code,
                        "selected_chain": chain_id,
                        "all_pdb_codes": pdb_codes,
                        "uniprot_ids": uniprot_ids,
                        "all_af_models": "; ".join(af_models),
                        **sel_info,
                    })
                    pdb_selected = True

        if pdb_selected:
            continue

        # 2b: AlphaFold DB -- try all UniProt IDs CAZy lists
        af_found_uids = []
        selected_path = None
        selected_uid = None
        for uid in all_uids:
            filepath = download_alphafold(uid, STRUCTURES_DIR)
            if filepath:
                af_found_uids.append(uid)
                if selected_path is None:
                    selected_uid = uid
                    selected_path = filepath

        if selected_path:
            structure_files.append((label, selected_path, selected_uid, "A"))
            if len(af_found_uids) > 1:
                print(f"  AF   {selected_uid} - {protein_name[:40]}"
                      f" (also found: {', '.join(af_found_uids[1:])})")
            else:
                print(f"  AF   {selected_uid} - {protein_name[:40]}")
            selection_rows.append({
                "family": family,
                "protein_name": protein_name,
                "organism": organism,
                "source": "AlphaFold",
                "selected_pdb": f"AF-{selected_uid}",
                "selected_chain": "A",
                "all_pdb_codes": pdb_codes,
                "uniprot_ids": uniprot_ids,
                "all_af_models": "; ".join(af_found_uids),
                "pdb_code": "",
                "chain_id": "A",
                "resolution": "",
                "method": "AlphaFold",
                "seq_length": "",
                "modeled_count": "",
                "coverage": 1.0,
            })
            continue

        genbank_ids = ([g.strip() for g in str(row.genbank_ids).split(";") if g.strip()]
                       if pd.notna(row.genbank_ids) else [])

        # 2c: UniParc -- UniProt entries with a sequence identical to the GenBank protein,
        # which CAZy does not always list, can have an AlphaFold DB model
        uniparc_path = uniparc_uid = duplicate_uid = None
        for genbank_id in genbank_ids:
            upi, uniparc_uids = uniparc_uniprot_ids(genbank_id, organism)
            duplicate_uid = next((u for u in uniparc_uids
                                  if u in cazy_uniprot_owner and u not in all_uids), None)
            if duplicate_uid:
                break
            for uid in [u for u in uniparc_uids if u not in all_uids][:20]:
                uniparc_path = download_alphafold(uid, STRUCTURES_DIR)
                if uniparc_path:
                    uniparc_uid = uid
                    break
            if uniparc_path:
                break

        if duplicate_uid:
            print(f"  DUP  {genbank_id} - {protein_name[:40]}: identical sequence to "
                  f"{cazy_uniprot_owner[duplicate_uid][:50]} (UniProt {duplicate_uid})")
            selection_rows.append({
                "family": family,
                "protein_name": protein_name,
                "organism": organism,
                "source": "NONE",
                "selected_pdb": "",
                "selected_chain": "",
                "all_pdb_codes": pdb_codes,
                "uniprot_ids": duplicate_uid,
                "all_af_models": "",
                "genbank_id": genbank_id,
                "uniprot_source": f"UniParc {upi}",
                "duplicate_of": cazy_uniprot_owner[duplicate_uid],
            })
            continue

        if uniparc_path:
            structure_files.append((label, uniparc_path, uniparc_uid, "A"))
            print(f"  AF   {uniparc_uid} via UniParc {upi} ({genbank_id}) - {protein_name[:40]}")
            selection_rows.append({
                "family": family,
                "protein_name": protein_name,
                "organism": organism,
                "source": "AlphaFold",
                "selected_pdb": f"AF-{uniparc_uid}",
                "selected_chain": "A",
                "all_pdb_codes": pdb_codes,
                "uniprot_ids": "; ".join([uniparc_uid] + all_uids),
                "all_af_models": uniparc_uid,
                "pdb_code": "",
                "chain_id": "A",
                "resolution": "",
                "method": "AlphaFold",
                "seq_length": "",
                "modeled_count": "",
                "coverage": 1.0,
                "genbank_id": genbank_id,
                "uniprot_source": f"UniParc {upi}",
            })
            continue

        # 2d: CAZyme3D, keyed by CAZy's GenBank accessions
        c3d_path = c3d_info = None
        for genbank_id in genbank_ids:
            c3d_path, c3d_info = download_cazyme3d(genbank_id, STRUCTURES_DIR)
            if c3d_path:
                break

        if c3d_path:
            structure_files.append((label, c3d_path, "", "A"))
            origin = (f"AlphaFold DB model of {c3d_info['uniprot']}" if c3d_info["uniprot"]
                      else "CAZyme3D AlphaFold2 prediction")
            print(f"  C3D  {genbank_id} ({origin}) - {protein_name[:40]}")
            selection_rows.append({
                "family": family,
                "protein_name": protein_name,
                "organism": organism,
                "source": "CAZyme3D",
                "selected_pdb": cazyme3d_struct_id(genbank_id),
                "selected_chain": "A",
                "all_pdb_codes": pdb_codes,
                "uniprot_ids": uniprot_ids,
                "all_af_models": "",
                "pdb_code": "",
                "chain_id": "A",
                "resolution": "",
                "method": origin,
                "seq_length": "",
                "modeled_count": "",
                "coverage": 1.0,
                "genbank_id": genbank_id,
                "plddt": c3d_info["plddt"],
            })
            continue

        print(f"  SKIP {family} {protein_name[:40]} - no PDB, AlphaFold DB or CAZyme3D model")
        selection_rows.append({
            "family": family,
            "protein_name": protein_name,
            "organism": organism,
            "source": "NONE",
            "selected_pdb": "",
            "selected_chain": "",
            "all_pdb_codes": pdb_codes,
            "uniprot_ids": uniprot_ids,
            "all_af_models": "",
        })

    source_of = {(r.family, r.protein_name, r.organism): (r.list_source, r.variant_of)
                 for r in df.itertuples()}
    for row in selection_rows:
        row["list_source"], variant_of = source_of[(row["family"], row["protein_name"], row["organism"])]
        row["variant_of"] = variant_of if pd.notna(variant_of) else ""

    # 2e: remove duplicate sequences (same protein from different CAZy strains)
    print("\n2e: checking for sequence duplicates across selected structures...")
    characterized = {(r.family, r.protein_name, r.organism) for r in cazy_df.itertuples()
                     if str(r.in_cazy_characterized).strip() == "yes"}
    selection_rows, structure_files = _dedup_by_sequence(selection_rows, structure_files,
                                                         characterized)

    # 2f: UniProt accessions, with their status, EC numbers and papers
    print("\n2f: collecting the UniProt accessions of the selected proteins (UniParc, UniProt)...")
    identity, references = collect_uniprot_accessions(selection_rows, listed_df)
    with open(UNIPROT_IDS_JSON, "w") as f:
        json.dump(identity, f, indent=1)
    with open(UNIPROT_REFERENCES_JSON, "w") as f:
        json.dump(references, f, indent=1)
    entries = [entry for info in identity.values() for entry in info["accessions"].values()]
    print(f"  {len(entries)} accessions for {len(identity)} proteins ("
          + ", ".join(f"{n} {status}" for status, n in Counter(e["status"] for e in entries).most_common())
          + f"), {sum(e['source'] == 'identical sequence' for e in entries)} of them for an identical sequence")
    print(f"  {sum(len(r) for r in references.values())} papers for {len(references)} proteins "
          f"-> {UNIPROT_REFERENCES_JSON}")
    without_active = sorted(s for s, info in identity.items()
                            if all(e["status"] == "deleted" for e in info["accessions"].values()))
    if without_active:
        print(f"  no active UniProt entry ({len(without_active)}): {', '.join(without_active)}")

    # 2g: domain annotations of the proteins' UniProt entries -- for trimming the N-terminus
    # (3h), the catalytic domain end (6b) and the associated CBMs
    uniprot_ids = sorted(({_primary_uniprot(row) for row in selection_rows}
                          | {acc for row in selection_rows
                             for acc in _split_ids(row.get("uniprot_active"))}) - {""})
    print(f"\n2g: fetching domain annotations for {len(uniprot_ids)} UniProt IDs...")
    domain_annotations = fetch_all_domain_annotations(uniprot_ids)
    with open(DOMAIN_ANNOTATIONS_JSON, "w") as f:
        json.dump(domain_annotations, f, indent=2)
    print(f"Saved domain annotations to {DOMAIN_ANNOTATIONS_JSON}")

    with_cat = sum(1 for a in domain_annotations.values() if a.get("catalytic_domain"))
    print(f"  With LPMO catalytic domain: {with_cat}")
    via_uniparc = {u: a["annotation_from"] for u, a in domain_annotations.items()
                   if a.get("annotation_from")}
    if via_uniparc:
        print("  Deleted UniProt entries annotated from an active entry with the same sequence: "
              + ", ".join(f"{u} <- {other}" for u, other in sorted(via_uniparc.items())))
    no_annotation = sorted(u for u, a in domain_annotations.items()
                           if not (a.get("catalytic_domain") or a.get("chain") or a.get("domains")))
    if no_annotation:
        print(f"  WARNING: no annotation at all (deleted UniProt entry without an active copy, or "
              f"a failed lookup) - 3h trims these models at the SignalP 6 signal peptide end, "
              f"else at the first histidine: {', '.join(no_annotation)}")
    if len(domain_annotations) > with_cat:
        print(f"  Without LPMO Pfam domain:   {len(domain_annotations) - with_cat}")
        for uid, a in domain_annotations.items():
            if a.get("catalytic_domain"):
                continue
            info = f"    {uid}:"
            if a.get("signal_peptide"):
                info += f" signal={a['signal_peptide']['start']}-{a['signal_peptide']['end']}"
            if a.get("chain"):
                info += f" chain={a['chain']['start']}-{a['chain']['end']}"
            pfam_domains = [d for d in a.get("domains", []) if d["source"] == "pfam"]
            if pfam_domains:
                info += " pfam=" + ", ".join(
                    f"{d['accession']}({d['start']}-{d['end']})" for d in pfam_domains)
            print(info)

    # 2h: EC numbers
    defaults = add_ec_numbers(selection_rows, cazy_df, identity, other_df)
    selected = [row for row in selection_rows if row["source"] != "NONE"]
    kinds = Counter("UniProt" if row["ec_source"].startswith("UniProt") else row["ec_source"] or "none"
                    for row in selected)
    print("\n2h: EC numbers of the selected proteins from "
          + ", ".join(f"{kind} {n}" for kind, n in kinds.most_common()))
    for label in ("UniProt", "Swiss-Prot", "BRENDA", "family default"):
        rows = [row for row in selected if row["ec_source"].startswith(label)]
        if rows:
            print(f"  {label}: " + ", ".join(f"{row['selected_pdb']} {row['ec_numbers']}" for row in rows))
    print("  family defaults: " + ", ".join(f"{family} {ec}" for family, ec in sorted(defaults.items())))

    # 2i: NCBI Taxonomy lineage of the source organisms (grouping by taxon in Django)
    from taxonomy import fetch_lineages, report as report_taxonomy
    organisms = {row["organism"] for row in selected}
    print(f"\n2i: NCBI Taxonomy lineage of {len(organisms)} organisms...")
    taxonomy = fetch_lineages(organisms)
    report_taxonomy(taxonomy)
    with open(TAXONOMY_JSON, "w") as f:
        json.dump(taxonomy, f, indent=1)
    print(f"Saved {TAXONOMY_JSON}")

    sel_df = pd.DataFrame(selection_rows).fillna("")

    # 2j: substrates and oxidised positions, each with the papers showing it (CAZy EC numbers,
    # Swiss-Prot, BRENDA); no source is preferred here, populate_db decides what to show
    from activities import collect_activities, find_brenda_file, report as report_activities
    brenda = find_brenda_file(DATA_DIR)
    print(f"\n2j: substrates and oxidised positions from CAZy, Swiss-Prot and BRENDA ({brenda})...")
    activities = collect_activities(sel_df, cazy_df, identity, references, taxonomy, brenda)
    report_activities(activities)
    with open(ACTIVITIES_JSON, "w") as f:
        json.dump(activities, f, indent=1)
    print(f"Saved {ACTIVITIES_JSON}")

    sel_df.to_csv(PDB_SELECTION_CSV, index=False)
    print(f"\nSaved PDB selection to {PDB_SELECTION_CSV}")

    sources = sel_df.source.value_counts()
    print("\nStructure sources:")
    for src, count in sources.items():
        print(f"  {src}: {count}")
    selected_df = sel_df[sel_df.source != "NONE"]
    print("Selected structures per list: " + ", ".join(
        f"{source} {n}" for source, n in selected_df.list_source.value_counts().items()))


# ---------------------------------------------------------------------------
# Step 3: Prepare structures
# ---------------------------------------------------------------------------

def _rebuild_loops(rows):
    """3a: rebuild unresolved internal loops of the selected crystal structures with MODELLER
    (fill_gaps.py). Returns {struct_id: path of the rebuilt structure}.

    Without a working MODELLER, the rebuilt files of an earlier run are reused where they
    fill exactly the gaps found now; the other structures continue unrebuilt."""
    from fill_gaps import MAX_GAP_LENGTH, fill_structure_gaps, modeller_version, plan_gaps

    plans = []
    for row in rows:
        if row.source != "PDB":
            continue
        path = _downloaded_structure(row.selected_pdb)
        if path:
            chain = _chain(row)
            plans.append((row.selected_pdb, path, chain, plan_gaps(path, chain)))

    for struct_id, _, chain, plan in plans:
        if plan.get("error"):
            print(f"  {struct_id}: {plan['error']}")
        for gap in plan["fill"]:
            print(f"  {struct_id} {chain}: gap {gap['start']}-{gap['end']} "
                  f"({gap['length']} aa {gap['sequence']}) -> rebuild")
        for gap in plan["skip"]:
            print(f"  {struct_id} {chain}: gap {gap['start']}-{gap['end']} "
                  f"({gap['length']} aa) left unfilled: {gap['reason']}")
    to_fill = [p for p in plans if p[3]["fill"] and not p[3].get("error")]
    print(f"\n  {len(to_fill)}/{len(plans)} PDB structures have internal gaps "
          f"<= {MAX_GAP_LENGTH} aa to rebuild")

    records_path = os.path.join(FILLED_DIR, "modelled_residues.json")
    version, error = modeller_version()
    if not version:
        previous = {}
        if os.path.exists(records_path):
            with open(records_path) as f:
                previous = json.load(f)
        reused = {}
        for struct_id, _, chain, plan in to_fill:
            rec = previous.get(struct_id, {})
            path = os.path.join(FILLED_DIR, f"{struct_id}.pdb")
            if rec.get("filled") == plan["fill"] and rec.get("chain") == chain \
                    and os.path.exists(path):
                reused[struct_id] = path
        print(f"  MODELLER not available ({error}): reusing {len(reused)} rebuilt structures "
              f"of an earlier run, {len(to_fill) - len(reused)} structures stay unrebuilt")
        return reused
    print(f"  MODELLER {version}")

    os.makedirs(FILLED_DIR, exist_ok=True)
    fill_ids = {struct_id for struct_id, *_ in to_fill}
    for struct_id, *_ in plans:   # drop outputs of structures that no longer need filling
        stale = os.path.join(FILLED_DIR, f"{struct_id}.pdb")
        if struct_id not in fill_ids and os.path.exists(stale):
            os.remove(stale)

    records = {}
    for struct_id, path, chain, plan in to_fill:
        try:
            rec = fill_structure_gaps(path, chain, FILLED_DIR, version, plan=plan)
        except Exception as e:
            print(f"  {struct_id}: FAILED - {e}")
            continue
        records[struct_id] = rec
        print(f"  {struct_id}: rebuilt {sum(g['length'] for g in rec['filled'])} residues, "
              f"DOPE {rec['dope']}, junction C-N {rec['junctions']} A, "
              f"closest contact {rec['min_contact']} A, "
              f"fixed atoms max shift {rec['fixed_max_dev']} A")

    with open(records_path, "w") as f:
        json.dump(records, f, indent=1)
    print(f"  Saved {records_path} ({len(records)} structures)")
    return {struct_id: os.path.join(FILLED_DIR, rec["output"])
            for struct_id, rec in records.items()}


def _curated_n_termini():
    """3h: {struct_id: first residue number or "keep"} from curation/n_terminus.csv."""
    if not os.path.exists(N_TERMINUS_CSV):
        return {}
    curated = pd.read_csv(N_TERMINUS_CSV, dtype=str, keep_default_na=False)
    return {r.struct_id.strip(): r.first_residue.strip() for r in curated.itertuples()
            if r.struct_id.strip()}


def _curated_c_termini():
    """3i: {struct_id: last residue number or "keep"} from curation/c_terminus.csv."""
    if not os.path.exists(C_TERMINUS_CSV):
        return {}
    curated = pd.read_csv(C_TERMINUS_CSV, dtype=str, keep_default_na=False)
    return {r.struct_id.strip(): r.last_residue.strip() for r in curated.itertuples()
            if r.struct_id.strip()}


def _signal_peptide_starts(rows, domain_annotations, manual):
    """3h: SignalP 6 for the models (AlphaFold DB, CAZyme3D; full precursor sequences) that
    have neither an annotated start nor a curated N-terminus. CAZy's kingdom sets the
    organism group: "other" for bacteria and archaea, "eukarya" for the rest (the viral
    fusolins are secreted by their eukaryotic hosts).

    Returns ({struct_id: residue number after the predicted signal peptide}, {struct_id of
    every model sent to SignalP}); the predictions go to signal_peptides.json and SignalP's
    own output to 03_preparation/signalp/<group>/."""
    from preprocess_structures import annotated_start, predict_signal_peptides, structure_sequence

    kingdom = {}
    if os.path.exists(CAZY_CSV):
        kingdom = {(r.family, r.protein_name, r.organism): r.kingdom for r in
                   _listed_entries().fillna("").itertuples()}
    groups, numbers = {"eukarya": {}, "other": {}}, {}
    for row in rows:
        struct_id = row.selected_pdb
        path = _downloaded_structure(struct_id)
        if (row.source not in MODEL_SOURCES or struct_id in manual or not path
                or annotated_start(domain_annotations.get(_primary_uniprot(row)))[0] is not None):
            continue
        numbers[struct_id], sequence = structure_sequence(path, _chain(row))
        group = ("other" if kingdom.get((row.family, row.protein_name, row.organism))
                 in ("Bacteria", "Archaea") else "eukarya")
        groups[group][struct_id] = sequence

    predictions, starts = {}, {}
    if numbers:
        print(f"\n3h: SignalP 6 for {len(numbers)} models without an annotated N-terminus")
    for group, sequences in groups.items():
        if not sequences:
            continue
        result, error = predict_signal_peptides(sequences, group,
                                                os.path.join(PREPROCESSING_DIR, "signalp", group))
        if error:
            print(f"  WARNING: SignalP 6 failed for the {group} group ({error}); these models "
                  f"fall back to the first histidine")
            continue
        for struct_id, p in result.items():
            end = p["sp_end"]
            after = numbers[struct_id][end] if end is not None and end < len(numbers[struct_id]) else None
            predictions[struct_id] = {**p, "organism_group": group, "first_residue_after": after}
            if after is not None:
                starts[struct_id] = after
            print(f"  {struct_id}: {p['prediction']}"
                  + (f", signal peptide 1-{end} (cleavage site Pr {p['cleavage_probability']}), "
                     f"next residue {after}" if end is not None else ""))
    if numbers:
        os.makedirs(PREPROCESSING_DIR, exist_ok=True)
        with open(SIGNAL_PEPTIDES_JSON, "w") as f:
            json.dump(predictions, f, indent=1)
    return starts, set(numbers)


def step_prepare():
    """3a-3i: one clean single-chain PDB file per selected structure, made from the
    MODELLER-rebuilt structure where there is one (see preprocess_structures.py)."""
    from preprocess_structures import preprocess_structure, tagged_residues

    rows = _selected_structures()
    domain_annotations = {}
    if os.path.exists(DOMAIN_ANNOTATIONS_JSON):
        with open(DOMAIN_ANNOTATIONS_JSON) as f:
            domain_annotations = json.load(f)
    else:
        print(f"  WARNING: {DOMAIN_ANNOTATIONS_JSON} not found (fetched in step 2g): "
              f"N-termini are trimmed by SignalP 6 (models) or to the first histidine")

    print("3a: rebuilding unresolved loops (MODELLER)")
    rebuilt = _rebuild_loops(rows)

    manual = _curated_n_termini()
    manual_c = _curated_c_termini()
    predicted, signalp_models = _signal_peptide_starts(rows, domain_annotations, manual)

    print(f"\n3b-3i: preprocessing {len(rows)} structures")
    os.makedirs(PREPROCESSED_DIR, exist_ok=True)
    done, n_metal = set(), 0
    for row in rows:
        struct_id = row.selected_pdb
        src_path = rebuilt.get(struct_id) or _downloaded_structure(struct_id)
        if not src_path:
            print(f"  {struct_id}: no structure file in {STRUCTURES_DIR} - skipped")
            continue
        out_path = os.path.join(PREPROCESSED_DIR, f"{struct_id}.pdb")
        # 3i: the tag annotation is in the downloaded entry; a MODELLER-rebuilt file (3a) has
        # only its ATOM records, so the tags are always read from the download.
        downloaded = _downloaded_structure(struct_id)
        tagged = tagged_residues(downloaded, _chain(row)) if downloaded else set()
        messages, metal = preprocess_structure(
            src_path, out_path, chain_id=_chain(row),
            annotation=domain_annotations.get(_primary_uniprot(row)),
            predicted_start=predicted.get(struct_id), manual=manual.get(struct_id),
            tagged=tagged, manual_c=manual_c.get(struct_id))
        if struct_id in rebuilt:
            messages.insert(0, "input: structure with rebuilt loops (3a)")
        if struct_id in signalp_models and not any("SignalP" in m for m in messages):
            with open(out_path) as f:
                first = next((line[17:20] for line in f if line.startswith("ATOM")), "")
            if first != "HIS" or any("first HIS" in m for m in messages):
                messages.append("WARNING: no annotation and no SignalP signal peptide ending "
                                "before a histidine - check the N-terminus (curation/n_terminus.csv)")
        for message in messages:
            print(f"  {struct_id}: {message}")
        done.add(f"{struct_id}.pdb")
        n_metal += bool(metal)

    _remove_stale(PREPROCESSED_DIR, done, ".pdb")
    print(f"\n  {len(done)} structures preprocessed ({len(rebuilt)} with rebuilt loops, "
          f"{n_metal} with an active-site metal) -> {PREPROCESSED_DIR}/")
    if len(done) < len(rows):
        raise StepError(f"{len(rows) - len(done)} selected structures could not be preprocessed")


# ---------------------------------------------------------------------------
# Step 4: Structural alignment
# ---------------------------------------------------------------------------

def step_align():
    """4: US-align -mm 4 on all preprocessed structures (see align_structures.py)."""
    from align_structures import copy_aligned_structures, run_usalign

    files = []
    for row in _selected_structures():
        path = os.path.join(PREPROCESSED_DIR, f"{row.selected_pdb}.pdb")
        if not os.path.exists(path):
            raise StepError(f"{path} not found. Run the 'prepare' step first.")
        files.append((_label(row), path))
    if len(files) < 2:
        raise StepError("need at least 2 structures for an alignment")

    results = run_usalign(files, ALIGNMENT_DIR)
    copied = copy_aligned_structures(ALIGNMENT_DIR, results["matrices"])

    print(f"\n  Alignment: {ALIGNMENT_FASTA}")
    print(f"  Matrices:  {os.path.join(ALIGNMENT_DIR, 'matrices.json')}")
    print(f"  {len(results['matrices'])} structures aligned, {copied} saved to {ALIGNED_DIR}/")
    if copied < len(files):
        raise StepError(f"US-align returned {copied} of {len(files)} structures")


# ---------------------------------------------------------------------------
# Step 5: Solvent accessibility and secondary structure (DSSP)
# ---------------------------------------------------------------------------

def step_dssp_sasa():
    """5a-5b: per-residue SASA and DSSP secondary structure of every aligned structure (see
    sasa_dssp.py).

    Output: 05_dssp_sasa/sasa.json
              {struct_id: {"residues": [numbers], "chains": [...], "names": [...],
                           "sasa": [A^2], "rsa": [relative SASA, %, null for unknown types]}}
            05_dssp_sasa/secondary_structure.json   (3-state, loaded by Django)
              {struct_id: ["H" | "E" | "-", ...]}
            05_dssp_sasa/dssp.json
              {struct_id: {"residues": [numbers], "ss8": "DSSP states, - for none",
                           "sheet": "sheet labels, - for none",
                           "bridges": [[residue, partner, ladder label], ...]}}
                (lowercase label: parallel ladder, uppercase: antiparallel)
    """
    from sasa_dssp import residue_sasa, run_dssp, three_state

    _require_mkdssp()
    sasa, ss, dssp, failed = {}, {}, {}, []
    for struct_id in _read_alignment():
        path = os.path.join(ALIGNED_DIR, f"{struct_id}.pdb")
        if not os.path.exists(path):
            raise StepError(f"{path} not found. Run the 'align' step first.")
        try:
            rows = residue_sasa(path)
            result = run_dssp(path)
        except Exception as e:
            failed.append(f"{struct_id} ({e})")
            continue
        sasa[struct_id] = {"residues": [r[0] for r in rows], "chains": [r[1] for r in rows],
                           "names": [r[2] for r in rows], "sasa": [r[3] for r in rows],
                           "rsa": [r[4] for r in rows]}
        ss[struct_id] = three_state(result["ss8"])
        dssp[struct_id] = {"residues": [r[1] for r in result["residues"]],
                           "ss8": result["ss8"], "sheet": result["sheet"],
                           "bridges": result["bridges"]}

    os.makedirs(os.path.dirname(SASA_JSON), exist_ok=True)
    for path, data in ((SASA_JSON, sasa), (SECONDARY_STRUCTURE_JSON, ss), (DSSP_JSON, dssp)):
        with open(path, "w") as f:
            json.dump(data, f)
    counts = Counter(c for states in ss.values() for c in states)
    rsa = [v for entry in sasa.values() for v in entry["rsa"] if v is not None]
    exposed = sum(v >= SURFACE_RSA for v in rsa)
    print(f"  SASA of {len(sasa)} structures -> {SASA_JSON} ({exposed} of {len(rsa)} residues "
          f"surface exposed, relative SASA >= {SURFACE_RSA} %)")
    print(f"  Secondary structure -> {SECONDARY_STRUCTURE_JSON}, {DSSP_JSON} (residues: helix "
          f"{counts['H']}, strand {counts['E']}, coil {counts['-']})")
    if failed:
        raise StepError(f"failed for {', '.join(failed)}")


# ---------------------------------------------------------------------------
# Step 6: Catalytic domain: disulfide bonds, domain end and histidine brace
# ---------------------------------------------------------------------------

def _column_of_position(seq, position):
    """1-based alignment column of the position-th residue (1-based) of an aligned sequence."""
    count = 0
    for i, c in enumerate(seq):
        if c != "-":
            count += 1
            if count == position:
                return i + 1
    return None


def _position_at_column(seq, column):
    """Residue position (1-based) at a 1-based alignment column; 0 if the column is a gap."""
    if column > len(seq) or seq[column - 1] == "-":
        return 0
    return sum(1 for c in seq[:column] if c != "-")


def _nearest_position(seq, column):
    """Residue position nearest to a 1-based alignment column (searching +-50 columns)."""
    for delta in range(0, 51):
        for sign in (-1, 1):
            col = column + sign * delta
            if 1 <= col <= len(seq) and seq[col - 1] != "-":
                return sum(1 for c in seq[:col] if c != "-")
    return None


def _disulfide_bonds(residues):
    """6a: disulfide bonds of one structure, [(index, index, SG-SG distance)] with 0-based
    residue indices in file order: the Cys pairs whose SG atoms lie within
    DISULFIDE_MAX_DISTANCE. From the coordinates, not from SSBOND records: models and some
    experimental files have none, and the records also list bonds to symmetry mates."""
    sulfurs = [(k, r[2]["SG"]) for k, r in enumerate(residues) if r[1] == "CYS" and "SG" in r[2]]
    bonds = []
    for (i, a), (j, b) in combinations(sulfurs, 2):
        distance = math.dist(a, b)
        if distance <= DISULFIDE_MAX_DISTANCE:
            bonds.append((i, j, round(distance, 2)))
    return bonds


def _consensus_domain_ends(aligned, residues, families, domain_annotations):
    """6b: catalytic domain end residue of every aligned structure.

    AlphaFold models (residue number = UniProt position) whose UniProt entry has an
    annotated catalytic domain vote with the alignment column of its end; the median column
    per family (clade) is the consensus. Each structure then gets
      annotation  its own annotated end (AlphaFold models with an annotation),
      consensus   the residue at, or nearest to, its family's consensus column -- or its
                  last residue when that is at most 15 residues further,
      fallback    its last residue (no consensus for the family).
    Returns {struct_id: {"end_residue", "family", "source"}}."""
    import statistics

    votes, own_end = {}, {}
    for struct_id, seq in aligned.items():
        if not struct_id.startswith("AF-") or not residues[struct_id]:
            continue
        cat = domain_annotations.get(struct_id[3:], {}).get("catalytic_domain")
        if not cat:
            continue
        numbers = [r[0] for r in residues[struct_id]]
        end = cat["end"]
        if end in numbers:
            position = numbers.index(end) + 1
        else:   # the annotated end lies beyond the model: clamp to its last residue
            end, position = numbers[-1], len(numbers)
        column = _column_of_position(seq, position)
        if column:
            family = families.get(struct_id, "unknown")
            votes.setdefault(family, []).append(column)
            own_end[struct_id] = end

    if not votes:
        raise StepError("no AlphaFold models with an annotated catalytic domain to vote")

    consensus = {}
    for family, columns in sorted(votes.items()):
        consensus[family] = int(statistics.median(columns))
        stdev = statistics.stdev(columns) if len(columns) > 1 else 0
        print(f"  {family:6s}: {len(columns):2d} votes, consensus col={consensus[family]}, "
              f"range={min(columns)}–{max(columns)}, stdev={stdev:.1f}")

    boundaries = {}
    for struct_id, seq in aligned.items():
        numbers = [r[0] for r in residues[struct_id]]
        if not numbers:
            continue
        family = families.get(struct_id, "unknown")
        end, source = numbers[-1], "fallback"
        if struct_id in own_end:
            end, source = own_end[struct_id], "annotation"
        elif family in consensus:
            position = (_position_at_column(seq, consensus[family])
                        or _nearest_position(seq, consensus[family]))
            if position and position <= len(numbers):
                # within 15 residues of the C-terminus: the structure is the catalytic domain
                end = numbers[-1] if len(numbers) - position <= 15 else numbers[position - 1]
                source = "consensus"
        boundaries[struct_id] = {"end_residue": end, "family": family, "source": source}

    counts = Counter(b["source"] for b in boundaries.values())
    print(f"  from annotation: {counts['annotation']}, from consensus: {counts['consensus']}, "
          f"fallback (last residue): {counts['fallback']}")
    return boundaries


def _end_column(seq, numbers, end_residue):
    """0-based alignment column of the last residue numbered <= end_residue."""
    last, k = None, 0
    for col, ch in enumerate(seq):
        if ch == "-":
            continue
        number = numbers[k] if k < len(numbers) else k + 1
        k += 1
        if end_residue is None or number <= end_residue:
            last = col
    return last


def _residue_at_or_before_column(seq, numbers, cap_col):
    """Number of the last residue at a 0-based alignment column <= cap_col."""
    last, k = None, 0
    for col, ch in enumerate(seq):
        if col > cap_col:
            break
        if ch == "-":
            continue
        last = numbers[k] if k < len(numbers) else k + 1
        k += 1
    return last


def _linker_snap(ss, numbers, end_residue, min_coil):
    """6d: if end_residue lies in a coil run of at least min_coil residues -- a domain
    linker -- the last helix/strand residue before that run, else end_residue."""
    if not min_coil or end_residue is None or not ss or len(ss) != len(numbers):
        return end_residue
    idx = max((i for i, n in enumerate(numbers) if n <= end_residue), default=None)
    if idx is None or ss[idx] not in ("U", "-"):
        return end_residue
    lo = hi = idx
    while lo >= 0 and ss[lo] in ("U", "-"):
        lo -= 1
    while hi < len(ss) and ss[hi] in ("U", "-"):
        hi += 1
    if (hi - lo - 1) < min_coil or lo < 0:
        return end_residue
    return numbers[lo]


def _cap_domain_ends(boundaries, aligned, residues):
    """6c-6d: tighten catalytic domain ends that run into a CBM, linker or disordered tail
    -- a consensus column landing in the CBM of a multi-domain protein, an over-long
    annotation, or a fallback to the last residue.

    6c  A family's consensus end column is the median end column of its members with an
        own annotation (else of all members), ignoring columns beyond 1.5x the median of
        all members. An end more than CAP_MARGIN columns past it moves to the last residue
        at or before that column.
    6d  An end lying in a coil run of at least LINKER_MIN_COIL residues moves back to the
        last helix/strand residue before the run (the 3-state secondary structure of step 5b).
    Changed entries get source "capped" and keep "uncapped_end_residue" and
    "original_source". Modifies boundaries in place."""
    import statistics

    secondary = _load_json(SECONDARY_STRUCTURE_JSON, "dssp_sasa")
    end_columns, family_columns, annotated_columns = {}, defaultdict(list), defaultdict(list)
    for struct_id, b in boundaries.items():
        numbers = [r[0] for r in residues[struct_id]]
        column = _end_column(aligned[struct_id], numbers, b["end_residue"])
        if column is None:
            continue
        end_columns[struct_id] = column
        family_columns[b["family"]].append(column)
        if b["source"] == "annotation":
            annotated_columns[b["family"]].append(column)

    consensus = {}
    for family, columns in family_columns.items():
        threshold = 1.5 * statistics.median(columns)
        pool = ([c for c in annotated_columns.get(family, []) if c <= threshold]
                or [c for c in columns if c <= threshold] or columns)
        consensus[family] = int(statistics.median(pool))

    n_ss = 0
    changed = []
    for struct_id, column in end_columns.items():
        b = boundaries[struct_id]
        numbers = [r[0] for r in residues[struct_id]]
        end = b["end_residue"]
        cap_column = consensus[b["family"]] + CAP_MARGIN
        new_end = (_residue_at_or_before_column(aligned[struct_id], numbers, cap_column)
                   if column > cap_column else end)
        if new_end is None:
            continue
        ss = secondary.get(struct_id)
        n_ss += ss is not None
        snapped = _linker_snap(ss, numbers, new_end, LINKER_MIN_COIL)
        if snapped is not None and snapped < new_end:
            new_end = snapped
        if new_end != end:
            changed.append((struct_id, b["family"], end, new_end, column))
            b["uncapped_end_residue"] = end
            b["original_source"] = b["source"]
            b["end_residue"] = new_end
            b["source"] = "capped"

    print("  Family consensus (median end column) and cap:")
    for family in sorted(consensus):
        print(f"    {family:6s} median col={consensus[family]:5d}  ->  cap col="
              f"{consensus[family] + CAP_MARGIN:5d}  (n={len(family_columns[family])})")
    if n_ss < len(end_columns):
        print(f"  WARNING: no secondary structure for {len(end_columns) - n_ss} structures in "
              f"{SECONDARY_STRUCTURE_JSON} (mkdssp failed in step 5) - linker check skipped "
              f"for them")
    print(f"  Capped {len(changed)} structure(s):")
    for struct_id, family, old, new, column in sorted(changed, key=lambda c: -c[4]):
        print(f"    {struct_id:<16} {family:<6} end {old} -> {new}  (end col was {column})")


def _keep_disulfide_bonds(boundaries, residues, bonds):
    """6e: an end lying inside a disulfide bond -- one Cys up to the end, the other after it
    -- moves to the second Cys, again while the new end cuts another bond, so that no bond
    ties the catalytic domain to the part left out (the UniProt domain annotations of AA9
    models end up to 25 residues before a Cys bonded to the domain). Applied after 6b-6d, so
    it overrides them. Changed entries get source "disulfide" and keep "end_before_disulfide"
    and "original_source" (unless 6c-6d set it). bonds: {struct_id: _disulfide_bonds result}.
    Modifies boundaries in place."""
    changed = []
    for struct_id, b in boundaries.items():
        numbers = [r[0] for r in residues[struct_id]]
        last = max((k for k, n in enumerate(numbers) if n <= b["end_residue"]), default=None)
        if last is None:
            continue
        new_last, crossing = last, []
        while True:     # each move goes further; ends when no bond crosses the end
            bond = next(((i, j) for i, j, _ in bonds.get(struct_id, []) if i <= new_last < j),
                        None)
            if bond is None:
                break
            new_last = bond[1]
            crossing.append(bond)
        if new_last == last:
            continue
        changed.append((struct_id, b["family"], b["end_residue"], numbers[new_last],
                        ", ".join(f"Cys{numbers[i]}-Cys{numbers[j]}" for i, j in crossing)))
        b["end_before_disulfide"] = b["end_residue"]
        b.setdefault("original_source", b["source"])
        b["source"] = "disulfide"
        b["end_residue"] = numbers[new_last]

    print(f"  Moved {len(changed)} end(s) to the second Cys of a bond:")
    for struct_id, family, old, new, crossing in changed:
        print(f"    {struct_id:<16} {family:<6} end {old} -> {new}  ({crossing})")


def _candidate_columns(counter):
    """0-based alignment columns at which to look for a residue in structures without a
    metal: the columns where metal-bound structures have it (most frequent first), then the
    columns up to BRACE_COLUMN_WINDOW next to them (nearest first)."""
    observed = [c for c, _ in counter.most_common()]
    neighbours = {c + d for c in observed
                  for d in range(-BRACE_COLUMN_WINDOW, BRACE_COLUMN_WINDOW + 1)}
    return observed + sorted(neighbours - set(observed),
                             key=lambda c: (min(abs(c - o) for o in observed), c))


def _his_nitrogen_distance(residue_a, residue_b):
    """Smallest distance (A) between the nitrogens of two residues."""
    n_a = [xyz for name, xyz in residue_a[2].items() if name.startswith("N")]
    n_b = [xyz for name, xyz in residue_b[2].items() if name.startswith("N")]
    return min((math.dist(a, b) for a in n_a for b in n_b), default=math.inf)


def _pseudo_cu(his1_atoms, his2_atoms):
    """6f: the pseudo-Cu, the centroid of the nitrogens that hold the copper (BRACE_NITROGENS:
    His-1 N, ND1, NE2 and His-2 ND1, NE2; 0.16 A from the Cu, median of the Cu structures), as
    [x, y, z] -- or None with fewer than two of them. Six decimals: rounded to the 0.001 A
    of the coordinates, a few of 7c's copper depths would change in their last digit. his1_atoms,
    his2_atoms: {atom name: (x, y, z)} of the two residues, {} for a missing one."""
    atoms = {"his1": his1_atoms, "his2": his2_atoms}
    coords = [atoms[his][name] for his, name in BRACE_NITROGENS if name in atoms[his]]
    if len(coords) < 2:
        return None
    return [round(sum(xyz[i] for xyz in coords) / len(coords), 6) for i in range(3)]


def _copper_position(metal_xyz, his1_atoms, his2_atoms, pseudo_cu):
    """6f: the copper position every later step uses (7c, 8c-8d, 12), as ([x, y, z], source):
    the structure's metal ("metal") when a ring nitrogen of each brace histidine lies within
    COPPER_MAX_RING_N_DISTANCE of it, else the pseudo-Cu ("pseudo-Cu"), else (None, None).
    his1_atoms, his2_atoms: {atom name: (x, y, z)} of the two residues."""
    if metal_xyz is not None:
        def holds(atoms):
            return any(name in atoms and math.dist(metal_xyz, atoms[name]) <= COPPER_MAX_RING_N_DISTANCE
                       for name in ("ND1", "NE2"))
        if holds(his1_atoms) and holds(his2_atoms):
            return [round(float(v), 3) for v in metal_xyz], "metal"
    if pseudo_cu is not None:
        return pseudo_cu, "pseudo-Cu"
    return None, None


def _histidine_braces(aligned, residues, metal_lines):
    """6f: His-1 and His-2 of the histidine brace of every aligned structure.

    metal      the two histidines with nitrogens closest to the active-site metal (kept in
               3g), in sequence order
    alignment  (no metal: apo crystal structures, models) His-1 is the first residue when it
               is a histidine, else the histidine at a His-1 column of the metal-bound
               structures; His-2 is the first histidine at their His-2 columns (most
               frequent first), then at the columns up to BRACE_COLUMN_WINDOW next to them
               (nearest first)
    Braces whose His-1 and His-2 nitrogens are more than BRACE_WARN_N_DIST apart are
    reported. Returns {struct_id: {"his1", "his2", "his1_column", "his2_column",
    "his_n_distance", "method", "metal", "pseudo_cu", "copper", "copper_source"}} with residue
    numbers, 1-based columns, the smallest His-1/His-2 nitrogen distance (A), the pseudo-Cu
    (_pseudo_cu, for every structure, also those with a metal) and the copper position the
    later steps use (_copper_position); None where not found."""
    columns = {}
    for struct_id, seq in aligned.items():
        cols = [i for i, c in enumerate(seq) if c != "-"]
        if len(cols) != len(residues[struct_id]):
            print(f"  WARNING: {struct_id}: {len(residues[struct_id])} residues in aligned/ "
                  f"but {len(cols)} in the alignment - skipped")
            continue
        columns[struct_id] = cols

    braces, his1_columns, his2_columns = {}, Counter(), Counter()
    for struct_id, cols in columns.items():
        metal_line = metal_lines.get(struct_id)
        if not metal_line:
            continue
        res = residues[struct_id]
        metal_xyz = (float(metal_line[30:38]), float(metal_line[38:46]),
                     float(metal_line[46:54]))
        his = sorted((min(math.dist(metal_xyz, xyz) for name, xyz in r[2].items()
                          if name.startswith("N")), k)
                     for k, r in enumerate(res)
                     if r[1] == "HIS" and any(name.startswith("N") for name in r[2]))
        if len(his) < 2:
            continue
        k1, k2 = sorted(k for _, k in his[:2])
        his1_columns[cols[k1]] += 1
        his2_columns[cols[k2]] += 1
        braces[struct_id] = {"his1": res[k1][0], "his2": res[k2][0],
                             "his1_column": cols[k1] + 1, "his2_column": cols[k2] + 1,
                             "his_n_distance": round(_his_nitrogen_distance(res[k1], res[k2]), 2),
                             "method": "metal",
                             "metal": (metal_line[76:78].strip().upper()      # element, e.g. CU
                                       or metal_line[17:20].strip().rstrip("0123456789")),
                             "pseudo_cu": _pseudo_cu(res[k1][2], res[k2][2])}
        braces[struct_id]["copper"], braces[struct_id]["copper_source"] = _copper_position(
            metal_xyz, res[k1][2], res[k2][2], braces[struct_id]["pseudo_cu"])

    his2_candidates = _candidate_columns(his2_columns)

    for struct_id, cols in columns.items():
        if struct_id in braces:
            continue
        res = residues[struct_id]
        at_column = {col: k for k, col in enumerate(cols)}
        k1 = 0 if res and res[0][1] == "HIS" else next(
            (at_column[c] for c, _ in his1_columns.most_common()
             if c in at_column and res[at_column[c]][1] == "HIS"), None)
        k2 = None
        if k1 is not None:
            k2 = next((at_column[c] for c in his2_candidates
                       if c in at_column and at_column[c] != k1
                       and res[at_column[c]][1] == "HIS"), None)
        braces[struct_id] = {
            "his1": res[k1][0] if k1 is not None else None,
            "his2": res[k2][0] if k2 is not None else None,
            "his1_column": cols[k1] + 1 if k1 is not None else None,
            "his2_column": cols[k2] + 1 if k2 is not None else None,
            "his_n_distance": (round(_his_nitrogen_distance(res[k1], res[k2]), 2)
                               if k2 is not None else None),
            "method": "alignment", "metal": None,
            "pseudo_cu": _pseudo_cu(res[k1][2] if k1 is not None else {},
                                    res[k2][2] if k2 is not None else {})}
        braces[struct_id]["copper"], braces[struct_id]["copper_source"] = _copper_position(
            None, {}, {}, braces[struct_id]["pseudo_cu"])

    n_metal = sum(b["method"] == "metal" for b in braces.values())
    print(f"  pseudo-Cu (centroid of the brace nitrogens): "
          f"{sum(b['pseudo_cu'] is not None for b in braces.values())} structures")
    sources = Counter(b["copper_source"] for b in braces.values())
    print(f"  copper position for steps 7c, 8 and 12: the metal in {sources['metal']}, the pseudo-Cu "
          f"in {sources['pseudo-Cu']}, none in {sources[None]}")
    for struct_id, b in sorted(braces.items()):
        if b["metal"] and b["copper_source"] != "metal":
            print(f"  NOTE: {struct_id}: its {b['metal']} is not held by the brace rings - pseudo-Cu used")
    incomplete = sorted(s for s, b in braces.items() if b["his1"] is None or b["his2"] is None)
    print(f"  from the metal: {n_metal}, from the alignment: {len(braces) - n_metal}, "
          f"incomplete: {len(incomplete)}")
    for method in ("metal", "alignment"):
        distances = [b["his_n_distance"] for b in braces.values()
                     if b["method"] == method and b["his_n_distance"] is not None]
        if distances:
            print(f"  His-1/His-2 nitrogen distance ({method}): "
                  f"{min(distances):.2f}–{max(distances):.2f} A")
    if his2_columns:
        print("  His-2 columns of metal-bound structures: "
              + ", ".join(f"{c + 1} ({n})" for c, n in his2_columns.most_common()))
    for struct_id, b in sorted(braces.items()):
        if b["his_n_distance"] is not None and b["his_n_distance"] > BRACE_WARN_N_DIST:
            print(f"  WARNING: {struct_id}: His-1 {b['his1']} and His-2 {b['his2']} nitrogens "
                  f"are {b['his_n_distance']} A apart - check the structure")
    for struct_id in incomplete:
        print(f"  no complete brace: {struct_id} {braces[struct_id]}")
    return braces


def step_catalytic_domain():
    """6a-6f: disulfide bonds, catalytic domain end and histidine brace of every aligned
    structure.

    Output: 06_catalytic_domain/disulfides.json
              {struct_id: [[residue, residue, SG-SG distance], ...]}
            06_catalytic_domain/catalytic_domain_boundaries.json
              {struct_id: {"end_residue", "family", "source", ...}}
            06_catalytic_domain/histidine_brace.json   (see _histidine_braces)
    """
    from preprocess_structures import active_site_metal

    aligned = _read_alignment()
    families = {row.selected_pdb: row.family for row in _selected_structures()}
    domain_annotations = {}
    if os.path.exists(DOMAIN_ANNOTATIONS_JSON):
        with open(DOMAIN_ANNOTATIONS_JSON) as f:
            domain_annotations = json.load(f)

    residues, metal_lines = {}, {}
    for struct_id in aligned:
        path = os.path.join(ALIGNED_DIR, f"{struct_id}.pdb")
        if not os.path.exists(path):
            raise StepError(f"{path} not found. Run the 'align' step first.")
        with open(path) as f:
            lines = f.readlines()
        residues[struct_id] = _residues(lines)
        metal_lines[struct_id] = active_site_metal(lines)

    print("6a: disulfide bonds")
    bonds = {struct_id: _disulfide_bonds(res) for struct_id, res in residues.items()}
    os.makedirs(ANNOTATION_DIR, exist_ok=True)
    with open(DISULFIDES_JSON, "w") as f:
        json.dump({struct_id: [[residues[struct_id][i][0], residues[struct_id][j][0], distance]
                               for i, j, distance in b]
                   for struct_id, b in bonds.items()}, f)
    print(f"  {sum(map(len, bonds.values()))} bonds in {sum(1 for b in bonds.values() if b)} "
          f"of {len(bonds)} structures -> {DISULFIDES_JSON}")

    print("\n6b: catalytic domain end (family consensus)")
    boundaries = _consensus_domain_ends(aligned, residues, families, domain_annotations)
    print("\n6c-6d: capping tails and linkers")
    _cap_domain_ends(boundaries, aligned, residues)
    print("\n6e: disulfide bonds across the end")
    _keep_disulfide_bonds(boundaries, residues, bonds)
    with open(BOUNDARIES_JSON, "w") as f:
        json.dump(boundaries, f, indent=2)
    print(f"  Saved {BOUNDARIES_JSON} ({len(boundaries)} structures)")

    print("\n6f: histidine brace")
    braces = _histidine_braces(aligned, residues, metal_lines)
    with open(BRACE_JSON, "w") as f:
        json.dump(braces, f, indent=1)
    print(f"  Saved {BRACE_JSON} ({len(braces)} structures)")


# ---------------------------------------------------------------------------
# Step 7: Beta-sandwich core, loop regions and binding surface
# ---------------------------------------------------------------------------

def _beta_core(aligned, families, n_residues):
    """7a: core strands of every structure (see beta_core.py); prints the strands found per
    family and the sheets they form. Returns ({struct_id: {name: (first, last, sheet)}} with
    0-based residue indices, {name: [first column, last column]} for the block strands)."""
    from beta_core import CORE_STRANDS, assign_strands, core_ladders, find_blocks, strand_name

    dssp = _load_json(DSSP_JSON, "dssp_sasa")
    ss8 = {s: d["ss8"] for s, d in dssp.items() if s in aligned}
    blocks = find_blocks(aligned, ss8, families, n_residues)
    if len(blocks) != CORE_STRANDS - 1:
        raise StepError(f"{len(blocks)} core strand blocks instead of {CORE_STRANDS - 1}: the loop "
                        f"definitions (beta_core.LOOPS) assume the strands of the 2026-09 data - "
                        f"check the blocks and update the definitions")
    columns = {strand_name(k): [b[0] + 1, b[1] + 1] for k, b in enumerate(blocks)}
    bridges = {}
    for s in ss8:
        index = {number: i for i, number in enumerate(dssp[s]["residues"])}
        bridges[s] = [(index[a], index[b], label) for a, b, label in dssp[s]["bridges"]]
    # first by the alignment, then checked against the ladders this first assignment shows
    by_columns = {s: assign_strands(aligned[s], ss8[s], dssp[s]["sheet"], n_residues[s], blocks,
                                    bridges[s]) for s in ss8}
    ladders = core_ladders(by_columns, bridges)
    strands = {s: assign_strands(aligned[s], ss8[s], dssp[s]["sheet"], n_residues[s], blocks,
                                 bridges[s], ladders) for s in ss8}

    names = [f"b{k}" for k in range(1, CORE_STRANDS + 1)]
    print(f"  7a beta-sandwich core: {CORE_STRANDS} strands (alignment columns): "
          + ", ".join(f"{n} {c[0]}-{c[1]}" for n, c in columns.items())
          + "; b2 = the hairpin partner of b1")
    members = Counter(families.get(s) for s in strands)
    print(f"    {'family':8s}{'n':>4s}  " + " ".join(f"{n:>4s}" for n in names))
    for family in sorted(members, key=lambda f: (len(f or ""), f or "")):
        found = Counter(n for s, st in strands.items() if families.get(s) == family for n in st)
        print(f"    {family or '?':8s}{members[family]:4d}  "
              + " ".join(f"{found[n]:4d}" for n in names))
    sheets = Counter()
    for st in strands.values():
        groups = defaultdict(list)
        for name, (_, _, label, _) in st.items():
            groups[label].append(name[1:])
        sheets[" | ".join(sorted(",".join(g) for g in groups.values()))] += 1
    print("    sheets (DSSP) of the core strands: "
          + "; ".join(f"{k} ({n})" for k, n in sheets.most_common(5)))
    print("    ladders: " + ", ".join(f"{a}{'||' if parallel else '-'}{b}"
                                      for (a, b), parallel in sorted(ladders.items())))
    for s in sorted(strands):
        changes = [f"{name} {_range(by_columns[s].get(name))} -> {_range(strands[s].get(name))}"
                   for name in sorted(set(by_columns[s]) | set(strands[s]), key=lambda n: int(n[1:]))
                   if by_columns[s].get(name) != strands[s].get(name)]
        if changes:
            print(f"    by pairing: {s} ({families.get(s)}) {', '.join(changes)}")
    return strands, columns


def _range(strand):
    """'first-last' (1-based) of a (first, last, ...) tuple of 0-based indices, or 'none'."""
    return f"{strand[0] + 1}-{strand[1] + 1}" if strand else "none"


def step_core_fold():
    """7a-7c: beta-sandwich core strands, loop regions and binding surface of every aligned
    structure (see beta_core.py and binding_surface.py). 7b the loops are the segments
    between given core strands (beta_core.LOOPS), within the catalytic domain; 7c the binding
    surface is the residues a plane tangent to the copper face touches -- the stand-in for a
    flat crystalline substrate (see binding_surface.py; it replaced the family templates of
    loops.py on 2026-09-23).

    Output: 07_core_fold/beta_core.json
      {"columns": {strand: [first, last]} (alignment columns; b2 has none),
       "structures": {struct_id: {strand: {"range": "start-end", "sheet": DSSP sheet label,
                                           "pieces": ["start-end", ...]}}}}
            07_core_fold/loop_regions.json
      {struct_id: {"loop_regions": {loop: ["start-end"]},
                   "binding_surface": {region: ["start-end", ...]} or null}}
      with ungapped sequence positions (1 = the first residue of the aligned structure);
      binding-surface residues outside every loop go under "other"
            07_core_fold/binding_surface.json  the plane, the copper depth, how evenly the
      contacts lie around it, and every residue with its distance to the plane, RSA and
      region (see binding_surface.py)
    """
    import statistics

    from beta_core import LOOPS, loop_regions
    from binding_surface import MAX_CU_DEPTH, assign_binding_surface, region_ranges

    aligned = _read_alignment()
    families = {row.selected_pdb: row.family for row in _selected_structures()}
    boundaries = _load_json(BOUNDARIES_JSON, "catalytic_domain")
    lines = _aligned_lines(aligned)

    # residues up to the catalytic domain end (numbered up to end_residue), else all
    n_residues = {}
    for struct_id, seq in aligned.items():
        end = boundaries.get(struct_id, {}).get("end_residue")
        n_residues[struct_id] = (len(seq.replace("-", "")) if end is None else
                                 sum(1 for r in _residues(lines[struct_id]) if r[0] <= end))

    strands, columns = _beta_core(aligned, families, n_residues)
    os.makedirs(os.path.dirname(BETA_CORE_JSON), exist_ok=True)
    with open(BETA_CORE_JSON, "w") as f:
        json.dump({"columns": columns,
                   "structures": {s: {name: {"range": f"{v[0] + 1}-{v[1] + 1}", "sheet": v[2],
                                             "pieces": [f"{a + 1}-{b + 1}" for a, b in v[3]]}
                                      for name, v in st.items()}
                                  for s, st in strands.items()}}, f, indent=1)
    print(f"    -> {BETA_CORE_JSON}")

    loops = {s: loop_regions(strands.get(s, {}), n_residues[s]) for s in aligned}

    # 7c: the region each sequence position belongs to, for the binding-surface records
    labels = {}
    for struct_id in aligned:
        for name, value in strands.get(struct_id, {}).items():
            for a, b in value[3]:
                for position in range(a + 1, b + 2):
                    labels.setdefault(struct_id, {})[position] = name.replace("b", "β")
        for loop, (a, b) in loops[struct_id].items():
            for position in range(a + 1, b + 2):
                labels.setdefault(struct_id, {})[position] = loop
    surfaces, no_surface = assign_binding_surface(
        aligned, ALIGNED_DIR, boundaries, _load_json(BRACE_JSON, "catalytic_domain"),
        _load_json(SASA_JSON, "dssp_sasa"), labels)

    regions = {s: {"loop_regions": {loop: [f"{a + 1}-{b + 1}"] for loop, (a, b) in loops[s].items()},
                   "binding_surface": (region_ranges(surfaces[s]) if s in surfaces else None)}
               for s in aligned}
    with open(LOOPS_JSON, "w") as f:
        json.dump(regions, f, indent=1)
    with open(BINDING_SURFACE_JSON, "w") as f:
        json.dump(surfaces, f, indent=1)

    print("  7b loops between core strands ("
          + ", ".join(f"{loop} {a}-{b or 'end'}" for loop, (a, b) in LOOPS.items())
          + "): structures with the loop, length range")
    members = Counter(families.get(s) for s in aligned)
    print(f"    {'family':8s}{'n':>4s}  " + "".join(f"{loop:>13s}" for loop in LOOPS))
    for family in sorted(members, key=lambda f: (len(f or ""), f or "")):
        cells = []
        for loop in LOOPS:
            lengths = [loops[s][loop][1] - loops[s][loop][0] + 1 for s in aligned
                       if families.get(s) == family and loop in loops[s]]
            cells.append(f"{len(lengths):3d} {min(lengths):3d}-{max(lengths):<3d}" if lengths
                         else f"{0:3d} -      ")
        print(f"    {family or '?':8s}{members[family]:4d}  " + "".join(f"{c:>13s}" for c in cells))
    print(f"  7c binding surface of {len(surfaces)} of {len(aligned)} structures "
          f"-> {BINDING_SURFACE_JSON}, {LOOPS_JSON}")
    print(f"    {'family':8s}{'n':>4s}{'residues':>12s}{'Cu below the plane':>21s}"
          f"{'around the Cu':>15s}{'recessed':>10s}   where the residues are")
    for family in sorted({families.get(s) for s in surfaces}, key=lambda f: (len(f or ""), f or "")):
        rows = [r for s, r in surfaces.items() if families.get(s) == family]
        counts = [r["n_residues"] for r in rows]
        depths = [r["copper"]["depth_below_plane"] for r in rows]
        spread = [r["contacts_around_copper"] for r in rows]
        where = Counter(residue["region"] or "other" for r in rows for residue in r["residues"])
        total = sum(where.values()) or 1
        print(f"    {family or '?':8s}{len(rows):4d}"
              f"{statistics.median(counts):7.0f} ({min(counts)}-{max(counts)})"
              f"{statistics.median(depths):14.1f} A ({min(depths):.1f}-{max(depths):.1f})"
              f"{statistics.median(spread):15.2f}"
              f"{sum(1 for r in rows if r['copper']['recessed']):10d}   "
              + " ".join(f"{k} {100 * v // total}%" for k, v in where.most_common(4)))
    recessed = [s for s, r in surfaces.items() if r["copper"]["recessed"]]
    if recessed:
        print(f"    a plane fits the copper face of {len(recessed)} structures poorly (a groove "
              f"or a bowl; the copper stays more than {MAX_CU_DEPTH} A below "
              f"any plane): " + ", ".join(sorted(recessed)))
    if no_surface:
        print(f"  WARNING: no binding surface for {', '.join(no_surface)}")


# ---------------------------------------------------------------------------
# Step 8: Second-sphere aromatic, brace interactions and motifs
# ---------------------------------------------------------------------------

def step_second_sphere():
    """8a-8e: second-sphere aromatic, histidine brace interactions, second-sphere motif, the
    solvent-facing axial residue and the motifs after the brace histidines of every aligned
    structure (see brace_interactions.py).

    8a The Phe/Tyr with an atom closest to the active-site metal; without a metal, the
        Phe/Tyr at the alignment columns of the metal-bound structures' aromatic (most
        frequent first, then up to BRACE_COLUMN_WINDOW columns next to them).
    8b Hydrogen bonds, pi-pi stacking and cation-pi contacts of His-1 and His-2
        (histidine_brace.json) for every structure with a complete brace.
    8c The second-sphere hydrogen-bonding motif (H-X8-Q-X-Y, E…N-X-F, ...) ending at the
        aromatic, for every structure with a brace and an aromatic.
    8d The residue on the copper's solvent-facing axial position (the Ala of AA10s).
    8c and 8d use the active-site metal, or the pseudo-Cu of 6f ("copper").
    8e The sequence motifs after the N-terminal His (N-His: H-G-X4-P-X2-R, H-G-X4-P, H-G) and
        the internal His (Internal-His: H-X-G-P, H-X2-G-P, H-X-G), for every structure with a
        complete brace.

    Output: 08_second_sphere/brace_interactions.json
      {struct_id: {"second_sphere_aromatic": {"name", "residue", "column", "method"} or null,
                   "interactions": {"hydrogen_bonds", "pi_pi_stacking", "cation_pi"} or null,
                   "copper": "metal" | "brace nitrogens" | null,
                   "second_sphere_motif": {"motif", "short", "residues": [{"role", "name",
                                           "residue", ...}]} or null,
                   "distal_axial": {"name", "residue", "atom", "distance", "angle"},
                                   {"name": null} (open) or null,
                   "brace_motifs": {"n_his": {"motif": "H-G-X4-P-X2-R" or null, "sequence",
                                              "residues": [{"role", "name", "residue",
                                                            "offset"}]},
                                    "internal_his": {...}} or null}}
    """
    from brace_interactions import (brace_motifs, detect_interactions,
                                    distal_axial_residue, nearest_aromatic, parse_residues,
                                    second_sphere_motif)
    import numpy as np
    from preprocess_structures import active_site_metal

    aligned = _read_alignment()
    braces = _load_json(BRACE_JSON, "catalytic_domain")
    families = {s: b.get("family")
                for s, b in _load_json(BOUNDARIES_JSON, "catalytic_domain").items()
                if isinstance(b, dict)}
    lines = _aligned_lines(aligned)

    results = {s: {"second_sphere_aromatic": None, "interactions": None, "copper": None,
                   "second_sphere_motif": None, "distal_axial": None, "brace_motifs": None}
               for s in aligned}
    ordered, parsed, columns = {}, {}, {}
    aromatic_columns = Counter()
    for struct_id, seq in aligned.items():
        ordered[struct_id] = _residues(lines[struct_id])
        parsed[struct_id] = parse_residues(os.path.join(ALIGNED_DIR, f"{struct_id}.pdb"))
        columns[struct_id] = _alignment_columns(seq, ordered[struct_id])
        metal_line = active_site_metal(lines[struct_id])
        if not metal_line:
            continue
        aromatic = nearest_aromatic(parsed[struct_id], (float(metal_line[30:38]),
                                                        float(metal_line[38:46]),
                                                        float(metal_line[46:54])))
        if not aromatic:
            continue
        number, name = aromatic
        k = next((k for k, r in enumerate(ordered[struct_id]) if r[0] == number), None)
        column = columns[struct_id][k] if columns[struct_id] and k is not None else None
        if column is not None:
            aromatic_columns[column] += 1
        results[struct_id]["second_sphere_aromatic"] = {
            "name": name, "residue": number,
            "column": column + 1 if column is not None else None, "method": "metal"}

    candidates = _candidate_columns(aromatic_columns)
    for struct_id, cols in columns.items():
        if results[struct_id]["second_sphere_aromatic"] or not cols:
            continue
        at_column = {col: k for k, col in enumerate(cols)}
        k = next((at_column[c] for c in candidates if c in at_column
                  and ordered[struct_id][at_column[c]][1] in ("PHE", "TYR")), None)
        if k is not None:
            number, name = ordered[struct_id][k][:2]
            results[struct_id]["second_sphere_aromatic"] = {
                "name": name, "residue": number, "column": cols[k] + 1, "method": "alignment"}

    for struct_id in aligned:
        brace = braces.get(struct_id, {})
        if brace.get("his1") is not None and brace.get("his2") is not None:
            results[struct_id]["interactions"] = detect_interactions(
                parsed[struct_id], brace["his1"], brace["his2"])
            results[struct_id]["brace_motifs"] = brace_motifs(
                ordered[struct_id], brace["his1"], brace["his2"])

    for struct_id in aligned:
        site, brace = results[struct_id], braces.get(struct_id, {})
        aromatic = site["second_sphere_aromatic"]
        if site["interactions"] is None or not aromatic:
            continue
        residues = ordered[struct_id]
        # the copper position decided in step 6f
        cu = np.array(brace["copper"], dtype=float) if brace.get("copper") else None
        site["copper"] = {"metal": "metal", "pseudo-Cu": "brace nitrogens"}.get(brace.get("copper_source"))
        index = next((k for k, r in enumerate(residues) if r[0] == aromatic["residue"]), None)
        stacked = {s["other_res_id"] for s in site["interactions"]["pi_pi_stacking"]
                   if s["brace_res_id"] == brace["his2"] and s["other_res_name"] == "HIS"}
        site["second_sphere_motif"] = second_sphere_motif(residues, index, stacked, cu)
        site["distal_axial"] = distal_axial_residue(residues, brace["his1"], brace["his2"],
                                                    aromatic["residue"], cu)

    os.makedirs(os.path.dirname(BRACE_INTERACTIONS_JSON), exist_ok=True)
    with open(BRACE_INTERACTIONS_JSON, "w") as f:
        json.dump(results, f, indent=1)

    methods = Counter((r["second_sphere_aromatic"] or {}).get("method") for r in results.values())
    print(f"  8a second-sphere aromatic: from the metal {methods['metal']}, from the alignment "
          f"{methods['alignment']}, none {methods[None]}")
    if aromatic_columns:
        print("     aromatic columns of metal-bound structures: "
              + ", ".join(f"{c + 1} ({n})" for c, n in aromatic_columns.most_common()))
    with_interactions = [r["interactions"] for r in results.values() if r["interactions"]]
    print(f"  8b brace interactions of {len(with_interactions)} structures: "
          + ", ".join(f"{sum(len(i[key]) for i in with_interactions)} {key}"
                      for key in ("hydrogen_bonds", "pi_pi_stacking", "cation_pi")))
    copper = Counter(r["copper"] for r in results.values() if r["copper"])
    print(f"  8c-8d copper: the metal in {copper['metal']}, the brace nitrogens' centroid in "
          f"{copper['brace nitrogens']}")

    def per_family(value_of, label):
        counts = defaultdict(Counter)
        for struct_id, r in results.items():
            value = value_of(r)
            if value is not None:
                counts[families.get(struct_id) or "?"][value] += 1
        print(f"  {label}:")
        for family in sorted(counts, key=lambda f: int(re.sub(r"\D", "", f) or 0)):
            print(f"     {family}: " + ", ".join(f"{v} {n}" for v, n in counts[family].most_common()))

    per_family(lambda r: (r["second_sphere_motif"] or {}).get("motif"), "8c second-sphere motifs")
    per_family(lambda r: None if r["distal_axial"] is None else (r["distal_axial"]["name"] or "open"),
               "8d residue on the solvent-facing axial position")
    for which, label in (("n_his", "N-His"), ("internal_his", "Internal-His")):
        per_family(lambda r: None if not (r["brace_motifs"] or {}).get(which)
                   else (r["brace_motifs"][which]["motif"] or "other"),
                   f"8e motifs after the {label}")
    print(f"  Saved {BRACE_INTERACTIONS_JSON}")


# ---------------------------------------------------------------------------
# Step 9: Distance matrix of the catalytic domains
# ---------------------------------------------------------------------------

def step_distance_matrix():
    """9: pairwise distances between the catalytic domains, from which the catalogue draws
    its dendrogram.

    Every alignment row is masked past the structure's catalytic domain end (6b-6e), so CBMs,
    linkers and tails do not pull otherwise similar catalytic domains apart. Biopython's
    DistanceCalculator (BLOSUM62) then compares each pair over the columns where neither
    sequence has a gap.

    Output: distance_matrix.csv (square matrix labelled by struct_id, in alignment order)
    """
    from Bio.Align import MultipleSeqAlignment
    from Bio.Phylo.TreeConstruction import DistanceCalculator
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    aligned = _read_alignment()
    boundaries = _load_json(BOUNDARIES_JSON, "catalytic_domain")
    lines = _aligned_lines(aligned)

    records, unmasked = [], []
    for struct_id, seq in aligned.items():
        end = boundaries.get(struct_id, {}).get("end_residue")
        if end is None:
            unmasked.append(struct_id)
            records.append(SeqRecord(Seq(seq), id=struct_id))
            continue
        numbers = [number for number, _, _ in _residues(lines[struct_id])]
        n_residues = sum(1 for c in seq if c != "-")
        if len(numbers) != n_residues:
            print(f"  WARNING: {struct_id}: {len(numbers)} residues in the structure, "
                  f"{n_residues} in the alignment - residues numbered 1-{n_residues}")
            numbers = list(range(1, n_residues + 1))
        masked, k = [], 0
        for c in seq:
            if c != "-":
                c = c if numbers[k] <= end else "-"
                k += 1
            masked.append(c)
        records.append(SeqRecord(Seq("".join(masked)), id=struct_id))
    if unmasked:
        print(f"  no catalytic domain end, whole chain compared: {', '.join(unmasked)}")

    print(f"  BLOSUM62 distances between {len(records)} catalytic domains...")
    matrix = DistanceCalculator("blosum62").get_distance(MultipleSeqAlignment(records))
    names = matrix.names
    table = pd.DataFrame([[matrix[a, b] for b in names] for a in names],
                         index=names, columns=names)
    os.makedirs(os.path.dirname(DISTANCE_MATRIX_CSV), exist_ok=True)
    table.to_csv(DISTANCE_MATRIX_CSV)
    print(f"  Saved {DISTANCE_MATRIX_CSV} ({len(names)} x {len(names)})")


# ---------------------------------------------------------------------------
# Step 10: Per-residue flexibility
# ---------------------------------------------------------------------------

def step_flexibility(run_cabsflex=False):
    """10a-10c: per-residue flexibility of every preprocessed structure.

    10a  bfactor: crystallographic B-factors of the CA atoms (PDB entries); residues rebuilt
        by MODELLER (occupancy 0) get none (null). plddt: AlphaFold pLDDT (AlphaFold DB and
        CAZyme3D models).
    10b  gnm: squared fluctuations from a ProDy Gaussian network model (20 modes).
    10c  cabsflex: RMSF of a CABSflex coarse-grained simulation, run in cabsflex_results/<id>/.
        A run is reused while it simulated the current structure (same input fingerprint in
        cabsflex_results/<id>/input.sha1); missing or outdated runs are (re)run in parallel
        only with --cabsflex (hours for the whole set).

    Output: flexibility_raw.json
      {struct_id: {"residues": [numbers], "bfactor" or "plddt": [...], "gnm": [...],
                   "cabsflex": [...]}}
    """
    try:
        import prody
        prody.confProDy(verbosity="none")
    except ImportError:
        raise StepError("ProDy not installed. Run: pip install prody")

    rows = _selected_structures()
    results = {}
    for i, row in enumerate(rows, 1):
        struct_id = row.selected_pdb
        pdb_file = os.path.join(PREPROCESSED_DIR, f"{struct_id}.pdb")
        if not os.path.exists(pdb_file):
            raise StepError(f"{pdb_file} not found. Run the 'prepare' step first.")

        print(f"  [{i}/{len(rows)}] {struct_id} ({row.source})", end="", flush=True)
        ca = prody.parsePDB(pdb_file).select("calpha")
        if ca is None:
            print(" - no CA atoms, skipped")
            continue
        entry = {"residues": ca.getResnums().tolist()}

        # 10a: B-factor / pLDDT
        betas = ca.getBetas().tolist()
        if row.source in MODEL_SOURCES:
            entry["plddt"] = [round(v, 2) for v in betas]
        else:
            entry["bfactor"] = [None if occupancy == 0 else round(v, 2)
                                for v, occupancy in zip(betas, ca.getOccupancies().tolist())]

        # 10b: GNM
        try:
            gnm = prody.GNM(struct_id)
            gnm.buildKirchhoff(ca)
            gnm.calcModes(n_modes=min(20, len(ca) - 1))
            entry["gnm"] = [round(v, 6) for v in prody.calcSqFlucts(gnm).tolist()]
            print(" gnm✓", end="", flush=True)
        except Exception as e:
            print(f"\n    WARNING: GNM failed: {e}", end="")

        results[struct_id] = entry
        print()

    # 10c: CABSflex
    todo = []
    for struct_id, entry in results.items():
        work_dir = os.path.join(CABSFLEX_DIR, struct_id)
        pdb_file = os.path.join(PREPROCESSED_DIR, f"{struct_id}.pdb")
        rmsf = (_parse_cabsflex_rmsf(work_dir) if _cabsflex_run_is_current(work_dir, pdb_file)
                else None)
        if rmsf and len(rmsf) == len(entry["residues"]):
            entry["cabsflex"] = [round(v, 4) for v in rmsf]
        else:
            todo.append(struct_id)
    print(f"\n  10c CABSflex: {len(results) - len(todo)} runs reused, "
          f"{len(todo)} missing or outdated")
    if todo and run_cabsflex:
        import multiprocessing
        n_workers = max(1, multiprocessing.cpu_count() - 1)
        print(f"  running {len(todo)} on {n_workers} workers")
        targets = [(struct_id, os.path.join(PREPROCESSED_DIR, f"{struct_id}.pdb"),
                    os.path.join(CABSFLEX_DIR, struct_id)) for struct_id in todo]
        with multiprocessing.Pool(n_workers) as pool:
            for struct_id, rmsf in pool.imap_unordered(_run_cabsflex_worker, targets):
                if rmsf and len(rmsf) == len(results[struct_id]["residues"]):
                    results[struct_id]["cabsflex"] = [round(v, 4) for v in rmsf]
                    print(f"    {struct_id}: cabsflex✓ ({len(rmsf)} residues)")
                else:
                    print(f"    {struct_id}: cabsflex failed")
    elif todo:
        print(f"  not run (add --cabsflex): {', '.join(todo)}")

    os.makedirs(os.path.dirname(FLEXIBILITY_JSON), exist_ok=True)
    with open(FLEXIBILITY_JSON, "w") as f:
        json.dump(results, f, indent=2)

    counts = {m: sum(1 for v in results.values() if m in v)
              for m in ("bfactor", "plddt", "gnm", "cabsflex")}
    print(f"\n  Done: {len(results)} structures  "
          + "  ".join(f"{m}: {n}" for m, n in counts.items()))
    print(f"  Saved to {FLEXIBILITY_JSON}")


def _cabsflex_input_key(pdb_file):
    """Fingerprint of the protein CABSflex simulates: atom names, residue names and numbers,
    and coordinates of the ATOM records (serial numbers, B-factors, the metal left out)."""
    digest = hashlib.sha1()
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM  "):
                digest.update(line[12:54].encode())
    return digest.hexdigest()


def _cabsflex_run_is_current(work_dir, pdb_file):
    """True if the CABSflex run in work_dir simulated pdb_file as it is now."""
    key_path = os.path.join(work_dir, CABSFLEX_INPUT_KEY)
    if not os.path.exists(key_path):
        return False
    with open(key_path) as f:
        return f.read().strip() == _cabsflex_input_key(pdb_file)


def _parse_cabsflex_rmsf(work_dir):
    """Compute per-residue RMSF from the CABSflex replica trajectory.

    CABSflex saves the full CA trajectory as replica.pdb (multi-model).
    We superpose all conformations and compute RMSF with ProDy —
    equivalent to what CABSflex computes internally.

    Returns a list of floats (Å) or None if parsing fails.
    """
    replica_pdb = os.path.join(work_dir, "output_pdbs", "replica.pdb")
    if not os.path.exists(replica_pdb):
        return None
    try:
        import prody
        ag = prody.parsePDB(replica_pdb, subset="ca")
        if ag is None or ag.numCoordsets() < 2:
            return None
        ens = prody.Ensemble("cabsflex")
        ens.setCoords(ag.getCoordsets()[0])
        ens.addCoordset(ag.getCoordsets())
        ens.superpose()
        return prody.calcRMSF(ens).tolist()
    except Exception:
        return None


def _run_cabsflex_worker(args):
    """Worker function for parallel CABSflex execution via multiprocessing.Pool.

    Args: (struct_id, pdb_file, out_dir); an earlier run in out_dir is replaced.
    Returns: (struct_id, rmsf_list_or_None)
    """
    import subprocess
    struct_id, pdb_file, out_dir = args
    shutil.rmtree(out_dir, ignore_errors=True)
    os.makedirs(out_dir)
    try:
        cmd = [
            "python3", "-m", "CABS", "flex",
            "-i", pdb_file,
            "--work-dir", out_dir,
            "--csv-output", "A",
            "--no-progress-bar",
        ]
        # check=False: CG2ALL reconstruction fails (missing binary) but
        # replica.pdb is written before that step — RMSF is extracted from it.
        subprocess.run(cmd, check=False, capture_output=True, timeout=3600)
        rmsf = _parse_cabsflex_rmsf(out_dir)
        if rmsf:
            with open(os.path.join(out_dir, CABSFLEX_INPUT_KEY), "w") as f:
                f.write(_cabsflex_input_key(pdb_file) + "\n")
        return (struct_id, rmsf)
    except Exception:
        return (struct_id, None)


# ---------------------------------------------------------------------------
# Step 11: Electrostatics
# ---------------------------------------------------------------------------

def _pdb2pqr_input(lines, his1, his2):
    """Structure lines prepared for pdb2pqr: brace His-1 renamed HIE and His-2 HID, and the
    methyl carbon (CZ) of N-methylated His-1 (HIC, renamed HIS in 3f) and the HIC header
    records removed -- the AMBER force field has no parameters for them."""
    out = []
    for line in lines:
        if line.startswith(("MODRES", "HET ", "HETNAM", "FORMUL")) and "HIC" in line:
            continue
        if line.startswith(("ATOM  ", "HETATM")) and line[17:20] == "HIS":
            if line[12:16].strip() == "CZ":
                continue
            try:
                number = int(line[22:26])
            except ValueError:
                number = None
            if number is not None and number == his1:
                line = line[:17] + "HIE" + line[20:]
            elif number is not None and number == his2:
                line = line[:17] + "HID" + line[20:]
        out.append(line)
    return out


def step_electrostatics():
    """11: PQR files (per-atom charges and radii) of every aligned structure with pdb2pqr
    and the AMBER force field.

    The brace histidines from histidine_brace.json (6f) get fixed protonation states:
    His-1 binds Cu through ND1 (and its amine), so NE2 carries the proton (HIE); His-2 binds
    through NE2, so ND1 carries it (HID). pdb2pqr leaves out the copper (no AMBER
    parameters), optimises the rest of the hydrogen-bond network and writes a chain column.
    All PQR files are regenerated.
    """
    import subprocess
    import tempfile

    struct_ids = list(_read_alignment())
    if not os.path.exists(BRACE_JSON):
        raise StepError(f"{BRACE_JSON} not found. Run the 'catalytic_domain' step first.")
    with open(BRACE_JSON) as f:
        braces = json.load(f)

    os.makedirs(PQR_DIR, exist_ok=True)
    done, failed, no_brace = set(), [], []
    for i, struct_id in enumerate(struct_ids, 1):
        pdb_path = os.path.join(ALIGNED_DIR, f"{struct_id}.pdb")
        pqr_path = os.path.join(PQR_DIR, f"{struct_id}.pqr")
        if os.path.exists(pqr_path):
            os.remove(pqr_path)
        brace = braces.get(struct_id, {})
        his1, his2 = brace.get("his1"), brace.get("his2")
        if his1 is None or his2 is None:
            no_brace.append(struct_id)

        print(f"  [{i}/{len(struct_ids)}] {struct_id}", end="", flush=True)
        with open(pdb_path) as f:
            lines = _pdb2pqr_input(f.readlines(), his1, his2)
        with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as tmp:
            tmp.writelines(lines)
        try:
            result = subprocess.run(
                ["pdb2pqr", "--ff=amber", "--chain", tmp.name, pqr_path],
                capture_output=True, text=True, timeout=600)
        finally:
            os.unlink(tmp.name)

        if result.returncode == 0 and os.path.exists(pqr_path) and os.path.getsize(pqr_path):
            done.add(f"{struct_id}.pqr")
            print(f" ✓ HIE:{his1}/HID:{his2}")
        else:
            failed.append(struct_id)
            print(f" FAILED: {(result.stderr or result.stdout).strip()[-200:]}")

    _remove_stale(PQR_DIR, done, ".pqr")
    print(f"\n  Done: {len(done)} PQR files in {PQR_DIR}/")
    if no_brace:
        print(f"  WARNING: no complete brace, His left to pdb2pqr: {', '.join(no_brace)}")
    if failed:
        raise StepError(f"pdb2pqr failed for {', '.join(failed)}")


# ---------------------------------------------------------------------------
# Step 12: Hole hopping
# ---------------------------------------------------------------------------

def step_hole_hopping():
    """12: hole hopping pathways from the protein surface to the Cu of every aligned
    structure with a histidine brace, with EHPath (see hole_hopping_ehpath.py).

    A structure is recalculated only when an input changed since the last run -- the aligned
    structure, the brace, the catalytic domain end, its relative SASA (step 5), its SASA
    threshold or the hole hopping scripts; the others keep their paths. Structures run in parallel.

    Output: 12_hole_hopping/hole_hopping.json
              {struct_id: {"inputs": hash, "paths": [{"rank", "path", "mean_residence_time",
                           "donor_residue", "donor_distance"}] or null, "note": str}}
            12_hole_hopping/<id>_paths.txt, EH_files/, EH_results/
    """
    import multiprocessing

    import hole_hopping_ehpath as hh

    aligned = _read_alignment()
    braces = _load_json(BRACE_JSON, "catalytic_domain")
    boundaries = _load_json(BOUNDARIES_JSON, "catalytic_domain")
    sasa = _load_json(SASA_JSON, "dssp_sasa")
    previous = {}
    if os.path.exists(HOLE_HOPPING_JSON):
        with open(HOLE_HOPPING_JSON) as f:
            previous = json.load(f)
    scripts_sha = hh.scripts_fingerprint()
    os.makedirs(HOLE_HOPPING_DIR, exist_ok=True)

    results, todo = {}, []
    for struct_id in aligned:
        pdb_path = os.path.join(ALIGNED_DIR, f"{struct_id}.pdb")
        if not os.path.exists(pdb_path):
            raise StepError(f"{pdb_path} not found. Run the 'align' step first.")
        brace = braces.get(struct_id, {})
        note = (f"excluded: {hh.EXCLUDED_STRUCT_IDS[struct_id]}" if struct_id in hh.EXCLUDED_STRUCT_IDS
                else "no histidine brace" if brace.get("his1") is None or brace.get("his2") is None
                else None)
        if note:
            hh.clear_outputs(struct_id, HOLE_HOPPING_DIR)
            results[struct_id] = {"inputs": None, "paths": None, "note": note}
            continue
        end_residue = boundaries.get(struct_id, {}).get("end_residue")
        threshold = hh.SASA_THRESHOLD_OVERRIDES.get(struct_id)
        if struct_id not in sasa:
            raise StepError(f"no SASA for {struct_id}. Run the 'dssp_sasa' step first.")
        rsa = (sasa[struct_id]["residues"], sasa[struct_id]["rsa"])
        inputs = hh.input_fingerprint(pdb_path, brace["his1"], brace["his2"], end_residue,
                                      threshold, scripts_sha, rsa)
        if previous.get(struct_id, {}).get("inputs") == inputs:
            results[struct_id] = previous[struct_id]
            continue
        results[struct_id] = {"inputs": inputs, "paths": None, "note": ""}
        todo.append((struct_id, pdb_path, brace["his1"], brace["his2"], end_residue,
                     HOLE_HOPPING_DIR, threshold, rsa, brace.get("copper")))

    for name in sorted(os.listdir(HOLE_HOPPING_DIR)):   # structures no longer aligned
        if name.endswith("_paths.txt") and name[:-len("_paths.txt")] not in aligned:
            hh.clear_outputs(name[:-len("_paths.txt")], HOLE_HOPPING_DIR)
            print(f"  removed outdated {name}")

    n_workers = max(1, min(len(todo), multiprocessing.cpu_count() - 1))
    print(f"  {len(aligned) - len(todo)} structures reused or skipped, {len(todo)} to run "
          f"on {n_workers} workers")
    if todo:
        with multiprocessing.Pool(n_workers) as pool:
            for i, (struct_id, paths, log) in enumerate(
                    pool.imap_unordered(hh.run_structure, todo), 1):
                notes = [line.strip() for line in log.splitlines() if line.strip()
                         and not line.startswith(("Processing protein", "Written:"))]
                results[struct_id]["paths"] = paths or None
                results[struct_id]["note"] = "; ".join(notes)[:300]
                print(f"  [{i}/{len(todo)}] {struct_id}: {len(paths) if paths else 'no'} paths"
                      + (f" ({results[struct_id]['note']})" if notes else ""))

    with open(HOLE_HOPPING_JSON, "w") as f:
        json.dump(dict(sorted(results.items())), f, indent=1)
    without = sorted(s for s, r in results.items() if not r.get("paths"))
    print(f"\n  {len(results) - len(without)}/{len(results)} structures with hole hopping paths "
          f"-> {HOLE_HOPPING_JSON}")
    for struct_id in without:
        print(f"  no paths: {struct_id} ({results[struct_id].get('note') or 'none found'})")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

STEPS = [   # (name, title, function) in pipeline order; the step number is the position
    ("scrape", "Listing the LPMOs (CAZy and the other sources)", step_scrape),
    ("fetch", "Fetching data", step_fetch),
    ("prepare", "Preparing structures", step_prepare),
    ("align", "Structural alignment (US-align)", step_align),
    ("dssp_sasa", "Secondary structure and solvent accessibility", step_dssp_sasa),
    ("catalytic_domain", "Catalytic domain, disulfide bonds and histidine brace",
     step_catalytic_domain),
    ("core_fold", "Beta-sandwich core, loop regions and binding surface", step_core_fold),
    ("second_sphere", "Second-sphere aromatic, brace interactions and motifs",
     step_second_sphere),
    ("distance_matrix", "Distance matrix of the catalytic domains", step_distance_matrix),
    ("flexibility", "Per-residue flexibility", step_flexibility),
    ("electrostatics", "Electrostatics (pdb2pqr)", step_electrostatics),
    ("hole_hopping", "Hole hopping pathways (EHPath)", step_hole_hopping),
]


def parse_steps(spec):
    """Selected step names in pipeline order from a comma-separated list of step names and
    numbers, numbers also as ranges: 'fetch,align', '3-7', '2,5-6'."""
    names = [name for name, _, _ in STEPS]
    selected = set()
    for token in (t.strip().lower() for t in spec.split(",")):
        if not token:
            continue
        if token in names:
            selected.add(token)
            continue
        first, dash, last = token.partition("-")
        try:
            first = int(first)
            last = int(last) if dash else first
        except ValueError:
            raise ValueError(f"unknown step '{token}'")
        if not 1 <= first <= last <= len(STEPS):
            raise ValueError(f"no step(s) '{token}' - steps are numbered 1-{len(STEPS)}")
        selected.update(names[first - 1:last])
    return [name for name in names if name in selected]


def main():
    step_list = "\n".join(f"  {number:>2}  {name:<17} {title}"
                          for number, (name, title, _) in enumerate(STEPS, 1))
    parser = argparse.ArgumentParser(
        description="LPMOcat data pipeline (see the module docstring for the sub-steps "
                    "and output files)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
Steps:
{step_list}

Examples:
  python pipeline.py                        # all steps
  python pipeline.py --steps scrape,fetch   # by name
  python pipeline.py --steps 3-12           # from the preparation on
  python pipeline.py --steps 10 --cabsflex  # flexibility, running missing CABSflex
  python pipeline.py --steps 1-2 --families AA9,AA10
""")
    parser.add_argument(
        "--steps", type=str, default=None,
        help="Comma-separated step names or numbers, ranges allowed, e.g. 'fetch,align' "
             "or '5-12' (default: all)")
    parser.add_argument(
        "--families", type=str, default=None,
        help="Comma-separated families for the scrape and fetch steps (e.g. AA15,AA16); "
             "their output files then list only these families. Default: all families.")
    parser.add_argument(
        "--cabsflex", action="store_true", default=False,
        help="Run CABSflex for structures without a current run in the flexibility step "
             "(slow: ~3 min/protein). Requires the CABS package.")
    parser.add_argument(
        "--data-dir", type=str, default=DEFAULT_DATA_DIR,
        help="Directory holding the step folders 01_proteins ... 12_hole_hopping and logs "
             "(default: data/ next to lpmocat_pipeline/)")
    args = parser.parse_args()

    try:
        steps = parse_steps(args.steps) if args.steps else [name for name, _, _ in STEPS]
    except ValueError as e:
        parser.error(str(e))
    families = [f.strip().upper() for f in args.families.split(",")] if args.families else None
    set_data_dir(os.path.abspath(args.data_dir))

    # Log to the terminal and a log file
    os.makedirs(LOG_DIR, exist_ok=True)
    log_path = os.path.join(LOG_DIR, f"pipeline_{datetime.now():%Y%m%d_%H%M%S}.log")
    log_file = open(log_path, "w", encoding="utf-8")
    sys.stdout = TeeStream(log_file, sys.__stdout__)
    sys.stderr = TeeStream(log_file, sys.__stderr__)

    step_kwargs = {"scrape": {"families": families}, "fetch": {"families": families},
                   "flexibility": {"run_cabsflex": args.cabsflex}}
    exit_code = 0
    try:
        print(f"Pipeline started at {datetime.now():%Y-%m-%d %H:%M:%S}")
        print(f"Data directory: {DATA_DIR}")
        print(f"Steps: {', '.join(steps)}")
        if families:
            print(f"Filtering to families: {families}")

        pipeline_start = datetime.now()
        for number, (name, title, function) in enumerate(STEPS, 1):
            if name not in steps:
                continue
            t0 = datetime.now()
            print(f"\n{'=' * 60}\nSTEP {number}: {title}\n{'=' * 60}")
            print(f"[{t0:%H:%M:%S}] Starting step {number} ({name})")
            try:
                function(**step_kwargs.get(name, {}))
            except StepError as e:
                print(f"\nERROR in step {number} ({name}): {e}")
                print("Pipeline stopped.")
                exit_code = 1
                break
            print(f"[{datetime.now():%H:%M:%S}] Step {number} ({name}) done in "
                  f"{datetime.now() - t0}")

        print("\n" + "=" * 60)
        if not exit_code:
            print(f"Pipeline complete. Total time: {datetime.now() - pipeline_start}")
        print(f"Log saved to {log_path}")
        print("=" * 60)
    finally:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        log_file.close()
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
