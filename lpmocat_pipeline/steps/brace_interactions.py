"""
Step 8 of the LPMOcat pipeline: the surroundings of the histidine brace in the aligned
structures, using His-1/His-2 from step 6f (histidine_brace.json).

  8a second-sphere aromatic: the Phe/Tyr with an atom closest to the active-site metal; in
      structures without a metal, the Phe/Tyr at the alignment column where the metal-bound
      structures have theirs (chosen in pipeline.step_second_sphere)
  8b brace interactions of His-1 and His-2:
        hydrogen bonds  brace ND1/NE2 to side-chain N/O acceptors of other residues, and N/O
                        donors of other residues to brace ND1/NE2, heavy atoms closer than
                        3.5 A
        pi-pi stacking  histidine ring centroids closer than 6.5 A with ring planes within
                        10 degrees of parallel or perpendicular
        cation-pi       Lys NZ / Arg NH1, NH2, NE centre closer than 6.0 A to a brace ring
                        centroid, projecting onto the ring normal (|cos| > 0.5)
  8c second-sphere hydrogen-bonding motif (second_sphere_motif), e.g. H-X8-Q-X-Y: the Gln/Glu
      two residues before the axial aromatic and the residue before it that holds it -- a His
      stacked on His-2, else a published sequence motif (KNOWN_MOTIFS), else a side-chain
      hydrogen bond. Without a Gln/Glu there (chitin-active AA10s: N-X-F), the Gln/Glu closest
      to the copper stands first: E…N-X-F.
  8d the residue on the copper's solvent-facing axial position (distal_axial_residue), the
      conserved Ala of AA10s (Hemsworth et al. 2013, doi:10.1021/ja402106e)
  8e sequence motifs after the two brace histidines (brace_motifs): the N-terminal His (N-His,
      His-1) followed by H-G-X4-P-X2-R, H-G-X4-P or H-G, the internal His (Internal-His, His-2)
      by H-X-G-P, H-X2-G-P or H-X-G (N_HIS_MOTIFS, INTERNAL_HIS_MOTIFS)

The copper position comes from step 6f (histidine_brace.json "copper"): the structure's metal
when the rings of both brace histidines hold it, else the pseudo-Cu (the centroid of the brace
nitrogens).
"""

import math

import numpy as np
from Bio.SeqUtils import seq1

AROMATIC_RINGS = {"HIS": ["CG", "CD2", "NE2", "CE1", "ND1"]}
CATION_ATOMS = {"LYS": ["NZ"], "ARG": ["NH1", "NH2", "NE"]}

# 8c: side-chain nitrogens and oxygens
POLAR_ATOMS = {"HIS": ("ND1", "NE2"), "ARG": ("NE", "NH1", "NH2"), "LYS": ("NZ",),
               "ASN": ("OD1", "ND2"), "GLN": ("OE1", "NE2"), "ASP": ("OD1", "OD2"),
               "GLU": ("OE1", "OE2")}
GLN_GLU = ("GLN", "GLU")
# the motif's first residue lies at most this many residues before the Gln/Glu
MOTIF_WINDOW = 16
# A, side-chain N/O distance of a first residue found by its hydrogen bond
MOTIF_HBOND = 3.5
# A, side-chain N/O of the Gln/Glu that heads an E…N-X-F motif to the copper
MOTIF_GLN_GLU_TO_CU = 8.0
# Published sequence motifs, (first residue, residues between, Gln/Glu, aromatic), used before
# the hydrogen-bond search since the first residue is not always in contact in the structures:
# R-X4-E-X-F and R-X5-E-X-F of cellulose-active and viral AA10s, N-X-E-X-Y of AA11s,
# Q-X2-Q-X-Y of AA13s and H-X5-8-Q-X-Y of AA9s (Vu et al. 2019, doi:10.1074/jbc.RA119.009509,
# and the works cited there), and Q-X8-Q-X-Y, which all AA17 structures have. AA11s also have
# an Arg four residues before their Glu; the aromatic keeps them N-X-E-X-Y. The AA9 His is
# mostly found by its stacking; in MtAA9 (5UFV) His161 lies 7.1 A from His-2 (H-X5-Q-X-Y).
KNOWN_MOTIFS = (("ARG", 4, "GLU", "PHE"), ("ARG", 5, "GLU", "PHE"), ("ASN", 1, "GLU", "TYR"),
                ("GLN", 2, "GLN", "TYR"), ("GLN", 8, "GLN", "TYR"),
                *(("HIS", between, "GLN", "TYR") for between in (5, 6, 7, 8)))
# 8d: an atom covers the solvent-facing axial position when it lies within this distance (A)
# of the copper and this angle (degrees) of the brace plane's normal, on the side away from the
# axial aromatic (the Ala CB of the H-X2-Q-X-Y AA10s lies 53-61 degrees off the normal, that of
# the other AA10s mostly within 35; equatorial atoms lie near 90)
DISTAL_MAX_DISTANCE = 4.5
DISTAL_MAX_ANGLE = 65
# 8e: motifs after the brace histidines, (name, {offset from the His: residue}), the first that
# fits wins. N-His: Gly next to it in AA10, AA13, AA15, AA17 and some AA9s and AA16s (its CA
# packs on b4, too close for a side chain); AA10, AA13, AA15 and AA16 also have a Pro and an Arg
# after it, the Arg in a salt bridge to an Asp of b9. Internal-His: H-X-G-P in AA9 and AA17
# (Laurent et al. 2019, doi:10.3390/ijms20246219, from TaAA9A; "in all LPMOs" it is not), with
# one more residue in some AA9s (H-X2-G-P), and the Gly alone (H-X-G) in AA11, AA13, AA15 and
# some AA10s; its CA lies 3.0-3.6 A from a backbone O of b8, too close for a side chain.
N_HIS_MOTIFS = (("H-G-X4-P-X2-R", {1: "GLY", 6: "PRO", 9: "ARG"}),
                ("H-G-X4-P", {1: "GLY", 6: "PRO"}),
                ("H-G", {1: "GLY"}))
INTERNAL_HIS_MOTIFS = (("H-X-G-P", {2: "GLY", 3: "PRO"}),
                       ("H-X2-G-P", {3: "GLY", 4: "PRO"}),
                       ("H-X-G", {2: "GLY"}))
# residues kept after each His (the His included): the longest motif's span
N_HIS_WINDOW, INTERNAL_HIS_WINDOW = 10, 5
# A, longest C-N distance of a peptide bond; a longer one is a chain break
PEPTIDE_BOND_MAX = 2.0


def parse_residues(pdb_path):
    """{residue number: {"res_name", "atoms": {atom name: (x, y, z)}}} of all ATOM and HETATM
    records."""
    residues = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                res_id = int(line[22:26].strip())
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except (ValueError, IndexError):
                continue
            residues.setdefault(res_id, {"res_name": res_name, "atoms": {}})
            residues[res_id]["atoms"][atom_name] = (x, y, z)
    return residues


def nearest_aromatic(residues, metal_xyz):
    """8a: (residue number, residue name) of the Phe/Tyr with an atom closest to the metal,
    or None."""
    best = None
    for number in sorted(residues):
        residue = residues[number]
        if residue["res_name"] not in ("PHE", "TYR"):
            continue
        distance = min(math.dist(metal_xyz, xyz) for xyz in residue["atoms"].values())
        if best is None or distance < best[0]:
            best = (distance, number, residue["res_name"])
    return (best[1], best[2]) if best else None


def _ring_centroid(res_name, atoms):
    names = AROMATIC_RINGS.get(res_name, [])
    coords = [atoms[n] for n in names if n in atoms]
    return np.mean(coords, axis=0) if len(coords) >= 3 else None


def _ring_normal(res_name, atoms):
    names = AROMATIC_RINGS.get(res_name, [])
    pts = [atoms[n] for n in names[:3] if n in atoms]
    if len(pts) < 3:
        return None
    v1 = np.array(pts[1]) - np.array(pts[0])
    v2 = np.array(pts[2]) - np.array(pts[0])
    n = np.cross(v1, v2)
    norm = np.linalg.norm(n)
    return n / norm if norm > 0 else None


def detect_interactions(all_res, brace_res1, brace_res2):
    """8b: {'hydrogen_bonds': [...], 'pi_pi_stacking': [...], 'cation_pi': [...]} of the two
    brace histidines, from parse_residues() output."""
    brace_ids = (brace_res1, brace_res2)
    brace_atoms = {rid: all_res[rid] for rid in brace_ids if rid in all_res}

    # ── 1. Hydrogen bonds ─────────────────────────────────────────────────────
    hbonds = []
    HIS_DONORS = ["ND1", "NE2"]
    HIS_ACCEPTORS = ["ND1", "NE2"]
    OTHER_ACCEPTORS = {"ND1", "NE2", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH"}
    OTHER_DONORS = {"ND1", "NE2", "N", "OH", "OG", "OG1"}
    CUTOFF = 3.5

    for bid in brace_ids:
        if bid not in brace_atoms:
            continue
        ba = brace_atoms[bid]["atoms"]
        for oid, oi in all_res.items():
            if oid in brace_ids:
                continue
            oa = oi["atoms"]
            rn = oi["res_name"]
            # brace → other
            for dn in HIS_DONORS:
                if dn not in ba:
                    continue
                for an, ac in oa.items():
                    if an in OTHER_ACCEPTORS:
                        d = math.dist(ba[dn], ac)
                        if d < CUTOFF:
                            hbonds.append({"brace_res_id": bid, "other_res_name": rn,
                                           "other_res_id": oid, "donor_atom": dn,
                                           "acceptor_atom": an, "distance": round(d, 2)})
            # other → brace
            for an, ac in oa.items():
                if an in OTHER_DONORS:
                    for acc in HIS_ACCEPTORS:
                        if acc in ba:
                            d = math.dist(ac, ba[acc])
                            if d < CUTOFF:
                                hbonds.append({"brace_res_id": bid, "other_res_name": rn,
                                               "other_res_id": oid, "donor_atom": an,
                                               "acceptor_atom": acc, "distance": round(d, 2)})

    # De-duplicate
    seen, unique_hbonds = set(), []
    for hb in hbonds:
        key = (min(hb["brace_res_id"], hb["other_res_id"]),
               max(hb["brace_res_id"], hb["other_res_id"]),
               round(hb["distance"], 1))
        if key not in seen:
            seen.add(key)
            unique_hbonds.append(hb)

    # ── 2. Pi-pi stacking ─────────────────────────────────────────────────────
    pi_stacking = []
    brace_centroids = []
    for bid in brace_ids:
        if bid not in brace_atoms:
            continue
        rn = brace_atoms[bid]["res_name"]
        c = _ring_centroid(rn, brace_atoms[bid]["atoms"])
        nb = _ring_normal(rn, brace_atoms[bid]["atoms"])
        if c is not None and nb is not None:
            brace_centroids.append((bid, rn, c, nb))

    for bid, brn, bc, nrm in brace_centroids:
        for oid, oi in all_res.items():
            if oid in brace_ids or oi["res_name"] not in AROMATIC_RINGS:
                continue
            oc = _ring_centroid(oi["res_name"], oi["atoms"])
            on = _ring_normal(oi["res_name"], oi["atoms"])
            if oc is None or on is None:
                continue
            dist = math.dist(bc, oc)
            if dist < 6.5:
                angle = np.arccos(np.clip(abs(np.dot(nrm, on)), 0, 1)) * 180 / np.pi
                if angle < 80 or angle > 100:
                    pi_stacking.append({"brace_res_id": bid, "brace_res_name": brn,
                                        "other_res_name": oi["res_name"], "other_res_id": oid,
                                        "distance": round(dist, 2), "angle": round(angle, 1)})

    # ── 3. Cation-pi interactions ─────────────────────────────────────────────
    cation_pi = []
    for bid, brn, bc, nrm in brace_centroids:
        for oid, oi in all_res.items():
            if oid in brace_ids or oi["res_name"] not in CATION_ATOMS:
                continue
            cat_names = CATION_ATOMS[oi["res_name"]]
            cat_coords = [np.array(oi["atoms"][n]) for n in cat_names if n in oi["atoms"]]
            if not cat_coords:
                continue
            cat_center = np.mean(cat_coords, axis=0)
            dist = math.dist(bc, cat_center)
            if dist < 6.0:
                vec = cat_center - np.array(bc)
                vn = np.linalg.norm(vec)
                if vn > 0 and abs(np.dot(vec / vn, nrm)) > 0.5:
                    cation_pi.append({
                        "brace_res_id": bid, "brace_res_name": brn,
                        "other_res_name": oi["res_name"], "other_res_id": oid,
                        "cation_atom": ", ".join(n for n in cat_names if n in oi["atoms"]),
                        "distance": round(dist, 2),
                        "projection": round(abs(np.dot(vec / vn, nrm)), 2),
                    })

    return {"hydrogen_bonds": unique_hbonds, "pi_pi_stacking": pi_stacking, "cation_pi": cation_pi}


def _atoms_of(residues, number):
    """{atom name: (x, y, z)} of the first residue numbered `number` in `residues`, else {}."""
    return next((atoms for n, _, atoms in residues if n == number), {})


def _is_hydrogen(atom_name):
    return atom_name.startswith("H") or (atom_name[:1].isdigit() and atom_name[1:2] == "H")


def _polar_distance(residue_a, residue_b):
    """Smallest distance (A) between the side-chain nitrogens and oxygens of two residues."""
    atoms_a = [residue_a[2][a] for a in POLAR_ATOMS.get(residue_a[1], ()) if a in residue_a[2]]
    atoms_b = [residue_b[2][b] for b in POLAR_ATOMS.get(residue_b[1], ()) if b in residue_b[2]]
    return min((math.dist(a, b) for a in atoms_a for b in atoms_b), default=math.inf)


def _spacer(between):
    """'X8-' for eight residues between two motif residues, 'X-' for one, '' for none."""
    return f"X{between if between > 1 else ''}-" if between else ""


def _motif_residue(role, residue, **extra):
    return {"role": role, "name": residue[1], "residue": residue[0], **extra}


def second_sphere_motif(residues, aromatic_index, stacked_his, cu):
    """8c: the second-sphere hydrogen-bonding motif of a structure.
    residues        [(number, name, {atom: xyz})] in sequence order
    aromatic_index  index of the axial aromatic (8a) in residues
    stacked_his     numbers of the His residues pi-stacked on His-2 (8b)
    cu              the copper position (step 6f, histidine_brace.json "copper"), or None
    The motif runs from its first residue to the aromatic, one-letter codes with the number of
    residues between them (H-X8-Q-X-Y); `short` leaves the spacing out (H-Q-Y). With a Gln/Glu
    two residues before the aromatic, the first residue is the one that holds it, searched up
    to MOTIF_WINDOW residues before it: a His stacked on His-2 (AA9, AA14 and the H-X2-Q-X-Y
    AA10s), else a KNOWN_MOTIFS residue, else the closest side-chain hydrogen bond partner
    (MOTIF_HBOND); without one the motif starts at the Gln/Glu (Q-X-Y). Otherwise the Gln/Glu
    whose side-chain N/O is closest to the copper, within MOTIF_GLN_GLU_TO_CU, heads it: E…N-X-F
    (the L2 Glu of chitin-active AA10s). Returns {"motif", "short", "residues": [{"role", "name",
    "residue", ...}]} with the roles first, gln_glu, middle (a residue two before the aromatic
    that is no Gln/Glu) and aromatic, or None."""
    if aromatic_index is None or aromatic_index < 2:
        return None
    aromatic, middle = residues[aromatic_index], residues[aromatic_index - 2]
    m = aromatic_index - 2
    aro, mid = seq1(aromatic[1]), seq1(middle[1])
    if middle[1] in GLN_GLU:
        window = range(max(0, m - MOTIF_WINDOW), m)
        first, found_by = None, None
        stacked = [i for i in window if residues[i][1] == "HIS" and residues[i][0] in stacked_his]
        if stacked:
            first, found_by = stacked[-1], "stacked on His-2"
        for name, between, gln_glu, aromatic_name in (KNOWN_MOTIFS if first is None else ()):
            i = m - between - 1
            if ((middle[1], aromatic[1]) == (gln_glu, aromatic_name) and i >= 0
                    and residues[i][1] == name):
                first, found_by = i, "sequence motif"
                break
        if first is None:
            contacts = [(d, i) for i in window if residues[i][1] in POLAR_ATOMS
                        for d in [_polar_distance(residues[i], middle)] if d <= MOTIF_HBOND]
            if contacts:
                first, found_by = min(contacts)[1], "hydrogen bond"
        entries = [_motif_residue("gln_glu", middle), _motif_residue("aromatic", aromatic)]
        if first is None:
            return {"motif": f"{mid}-X-{aro}", "short": f"{mid}-{aro}", "residues": entries}
        head = seq1(residues[first][1])
        return {"motif": f"{head}-{_spacer(m - first - 1)}{mid}-X-{aro}",
                "short": f"{head}-{mid}-{aro}",
                "residues": [_motif_residue("first", residues[first], found_by=found_by)] + entries}

    best = None
    for i, (_, name, atoms) in enumerate(residues if cu is not None else ()):
        if name in GLN_GLU:
            d = min((math.dist(cu, atoms[a]) for a in POLAR_ATOMS[name] if a in atoms),
                    default=math.inf)
            if d <= MOTIF_GLN_GLU_TO_CU and (best is None or d < best[0]):
                best = (d, i)
    entries = [_motif_residue("middle", middle), _motif_residue("aromatic", aromatic)]
    if best is None:
        return {"motif": f"{mid}-X-{aro}", "short": f"{mid}-{aro}", "residues": entries}
    head_residue = residues[best[1]]
    head = seq1(head_residue[1])
    return {"motif": f"{head}…{mid}-X-{aro}", "short": f"{head}-{mid}-{aro}",
            "residues": [_motif_residue("gln_glu", head_residue, cu_distance=round(best[0], 1))]
            + entries}


def distal_axial_residue(residues, his1, his2, aromatic_number, cu):
    """8d: the residue on the copper's solvent-facing axial position, opposite the axial
    aromatic: the one with a heavy side-chain atom (Gly: CA) closest to the copper among those
    within DISTAL_MAX_DISTANCE of it and DISTAL_MAX_ANGLE of the normal of the plane through
    His-1 N, His-1 ND1 and His-2 NE2. Returns {"name", "residue", "atom", "distance", "angle"},
    {"name": None} when no residue covers the position (open), or None without the brace atoms,
    the aromatic or the copper."""
    his1_atoms, his2_atoms = _atoms_of(residues, his1), _atoms_of(residues, his2)
    aromatic = _atoms_of(residues, aromatic_number)
    if cu is None or "N" not in his1_atoms or "ND1" not in his1_atoms or "NE2" not in his2_atoms:
        return None
    normal = np.cross(np.subtract(his1_atoms["N"], his2_atoms["NE2"]),
                      np.subtract(his1_atoms["ND1"], his2_atoms["NE2"]))
    ring = [aromatic[a] for a in ("CZ", "CE1", "CE2") if a in aromatic]
    reference = aromatic.get("OH") or (tuple(np.mean(ring, axis=0)) if ring else None)
    if reference is None or not np.linalg.norm(normal):
        return None
    normal = normal / np.linalg.norm(normal)
    if np.dot(np.subtract(reference, cu), normal) > 0:
        normal = -normal
    best = None
    for number, name, atoms in residues:
        if number in (his1, his2):
            continue
        for atom, xyz in atoms.items():
            if (_is_hydrogen(atom) or atom in ("N", "C", "O", "OXT")
                    or (atom == "CA" and name != "GLY")):
                continue
            vector = np.subtract(xyz, cu)
            distance = float(np.linalg.norm(vector))
            if not 0 < distance <= DISTAL_MAX_DISTANCE:
                continue
            cosine = max(-1.0, min(1.0, float(np.dot(vector, normal)) / distance))
            angle = math.degrees(math.acos(cosine))
            if angle <= DISTAL_MAX_ANGLE and (best is None or distance < best[0]):
                best = (distance, angle, number, name, atom)
    if best is None:
        return {"name": None}
    distance, angle, number, name, atom = best
    return {"name": name, "residue": number, "atom": atom, "distance": round(distance, 1),
            "angle": round(angle)}


def _following(residues, index, length):
    """The residue at index and the ones after it, up to length residues in all, ending before
    a chain break."""
    window = [residues[index]]
    for k in range(index + 1, min(index + length, len(residues))):
        c, n = residues[k - 1][2].get("C"), residues[k][2].get("N")
        if c is None or n is None or math.dist(c, n) > PEPTIDE_BOND_MAX:
            break
        window.append(residues[k])
    return window


def _his_motif(residues, number, motifs, length):
    """The motif after the His with this residue number: {"motif": name or None, "sequence":
    one-letter codes of the His and the residues after it, "residues": [{"role", "name",
    "residue", "offset"}] of the motif residues after the His}, or None without that His."""
    index = next((k for k, r in enumerate(residues) if r[0] == number and r[1] == "HIS"), None)
    if index is None:
        return None
    window = _following(residues, index, length)
    sequence = "".join(seq1(r[1]) for r in window)
    for motif, positions in motifs:
        if all(k < len(window) and window[k][1] == name for k, name in positions.items()):
            return {"motif": motif, "sequence": sequence,
                    "residues": [_motif_residue(name.lower(), window[k], offset=k)
                                 for k, name in sorted(positions.items())]}
    return {"motif": None, "sequence": sequence, "residues": []}


def brace_motifs(residues, his1, his2):
    """8e: the motifs after the N-terminal His (N_HIS_MOTIFS, over N_HIS_WINDOW residues) and
    the internal His (INTERNAL_HIS_MOTIFS, INTERNAL_HIS_WINDOW), residues in sequence order.
    Returns {"n_his": ..., "internal_his": ...} (see _his_motif); positions after a chain break
    count as missing."""
    return {"n_his": _his_motif(residues, his1, N_HIS_MOTIFS, N_HIS_WINDOW),
            "internal_his": _his_motif(residues, his2, INTERNAL_HIS_MOTIFS, INTERNAL_HIS_WINDOW)}
