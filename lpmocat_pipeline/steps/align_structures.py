"""
Step 4 of the LPMOcat pipeline: structural alignment with US-align.

US-align -mm 4 aligns all preprocessed structures at once (multiple structural alignment)
and writes the sequence alignment implied by the superposition, one transformation matrix
per structure and the superposed structures (result.<n>.pdb), saved as aligned/<id>.pdb.
"""

import glob
import json
import os
import shutil
import subprocess

from paths import TOOLS_DIR

USALIGN_BIN = os.path.join(TOOLS_DIR, "USalign", "USalign")


def run_usalign(structure_files, outdir):
    """Run US-align multiple structural alignment.

    Args:
        structure_files: list of (label, filepath) tuples, all in one directory
        outdir: directory for results

    Returns:
        dict with keys: fasta_text, matrices (list of dicts), stdout
    """
    os.makedirs(outdir, exist_ok=True)
    for stale in glob.glob(os.path.join(outdir, "result.*.pdb")):
        os.remove(stale)

    # US-align -mm 4 requires: USalign -dir folder/ list_file -mm 4
    # The list file must contain just filenames relative to the -dir folder.
    struct_dir = os.path.dirname(structure_files[0][1])
    list_path = os.path.join(outdir, "structure_list.txt")
    with open(list_path, "w") as f:
        for label, filepath in structure_files:
            f.write(os.path.basename(filepath) + "\n")

    matrix_path = os.path.join(outdir, "matrices.txt")
    result_prefix = os.path.join(outdir, "result")

    cmd = [
        USALIGN_BIN,
        "-dir", struct_dir + "/",
        list_path,
        "-mm", "4",
        "-mol", "prot",
        "-m", matrix_path,
        "-o", result_prefix,
    ]

    print(f"Running US-align on {len(structure_files)} structures...")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

    stdout_text = result.stdout + result.stderr
    with open(os.path.join(outdir, "usalign_stdout.txt"), "w") as f:
        f.write(stdout_text)
    if result.returncode != 0:
        raise RuntimeError(f"US-align exited with code {result.returncode}: "
                           f"{stdout_text.strip()[-500:]}")

    # Extract FASTA from stdout
    fasta_lines = []
    in_fasta = False
    for line in stdout_text.splitlines():
        if line.startswith(">"):
            in_fasta = True
        if in_fasta:
            fasta_lines.append(line)
            if not line.startswith(">") and line.strip():
                in_fasta = False  # sequence line after header
    if not fasta_lines:
        raise RuntimeError("US-align wrote no alignment - see usalign_stdout.txt")

    fasta_text = "\n".join(fasta_lines) + "\n"
    fasta_path = os.path.join(outdir, "alignment.fasta")
    with open(fasta_path, "w") as f:
        f.write(fasta_text)

    # Parse transformation matrices
    matrices = parse_usalign_matrices(matrix_path, structure_files)

    matrices_path = os.path.join(outdir, "matrices.json")
    with open(matrices_path, "w") as f:
        json.dump(matrices, f, indent=2)

    return {
        "fasta_text": fasta_text,
        "matrices": matrices,
        "stdout": stdout_text,
    }


def parse_usalign_matrices(matrix_path, structure_files):
    """Parse US-align matrix output into Mol*-compatible 4x4 matrices.

    Returns list of dicts: [{label, source_file, matrix_4x4}, ...]
    The reference structure (last one, marked with * by US-align) gets
    an identity matrix.
    """
    with open(matrix_path) as f:
        content = f.read()

    blocks = content.split("The rotation matrix to rotate ")[1:]

    result = {}
    reference_target = None

    for block in blocks:
        header = block.split("\n")[0]
        source = header.split(" to ")[0].strip()
        target = header.split(" to ")[1].split(" -")[0].strip()
        reference_target = target

        t = [0.0, 0.0, 0.0]
        u = [[0.0] * 3 for _ in range(3)]

        for line in block.split("\n"):
            parts = line.split()
            if len(parts) == 5 and parts[0] in ("0", "1", "2"):
                i = int(parts[0])
                t[i] = float(parts[1])
                u[i] = [float(parts[2]), float(parts[3]), float(parts[4])]

        # Row-major 4x4 matrix for Mol* TransformStructureConformation
        mat4x4 = [
            u[0][0], u[0][1], u[0][2], t[0],
            u[1][0], u[1][1], u[1][2], t[1],
            u[2][0], u[2][1], u[2][2], t[2],
            0, 0, 0, 1,
        ]

        result[source] = mat4x4

    # Add identity matrix for reference structure
    if reference_target:
        result[reference_target] = [
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1,
        ]

    # Map back to structure labels
    matrices = []
    for label, filepath in structure_files:
        basename = os.path.basename(filepath)
        # US-align uses "filename:chain" as key
        matched_key = None
        for key in result:
            if key.startswith(basename):
                matched_key = key
                break

        if matched_key:
            matrices.append({
                "label": label,
                "source_file": basename,
                "usalign_key": matched_key,
                "matrix_4x4": result[matched_key],
            })
        else:
            print(f"  WARNING: no matrix found for {label} ({basename})")

    return matrices


def copy_aligned_structures(outdir, matrices):
    """Save the superposed structures result.<n>.pdb as aligned/<id>.pdb.

    Earlier versions are replaced, and aligned files of structures that are no longer part
    of the alignment are removed, so aligned/ always matches alignment.fasta.
    Returns the number of structures copied."""
    aligned_dir = os.path.join(outdir, "aligned")
    os.makedirs(aligned_dir, exist_ok=True)

    copied, current = 0, set()
    for i, entry in enumerate(matrices):
        src = os.path.join(outdir, f"result.{i}.pdb")
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(aligned_dir, entry["source_file"]))
            current.add(entry["source_file"])
            copied += 1

    for name in os.listdir(aligned_dir):
        if name.endswith(".pdb") and name not in current:
            os.remove(os.path.join(aligned_dir, name))
            print(f"  removed {name} from aligned/ (not in the alignment)")
    return copied
