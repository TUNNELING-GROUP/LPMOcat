#!/usr/bin/env python
"""
Step 7c of the LPMOcat pipeline: the binding surface of every structure, from a TANGENT
PLANE at the copper face -- the stand-in for a flat crystalline substrate. It replaced the
hand-picked family templates of loops.py on 2026-09-23.

Per structure, on its own coordinates (the alignment is not used; that the residues land in
similar alignment columns is something to look at afterwards, not an input):

  1. the catalytic domain only (step 6), so a CBM or a linker cannot decide where the plane
     lies;
  2. the copper: step 6f's copper position (histidine_brace.json "copper": the active-site
     metal when the rings of both brace histidines hold it, else the pseudo-Cu, the centroid
     of the brace nitrogens);
  3. plane orientations sampled on a sphere. For each direction the plane is pushed in until
     it touches the outermost residue, so it never cuts the protein (ALLOW_ABOVE residues may
     stick through it, for the odd protruding side chain);
  4. of the directions holding the copper at most MAX_CU_DEPTH below the plane, the one with
     the best score wins, then it is refined locally. Score = the residues it touches (their
     outermost heavy atom within SLAB of the plane) taken down by BALANCE when they all lie on
     ONE SIDE of the copper instead of around it -- a plane resting on one shoulder of the
     protein loses to one a substrate could really lie on;
  5. binding surface = the touching residues that have any solvent accessibility at all
     (step 5a; no threshold -- a quarter of them are below the 25 % the app calls exposed).

Structures whose copper never comes within MAX_CU_DEPTH of any tangent plane are marked
`recessed` and fitted with the limit opened to the shallowest depth they allow: their site is
a groove or a bowl (AA13, AA14 and above all AA17, whose copper sits 7-8 A down behind a rim
of loops), and a plane describes it poorly. The copper depth and the spread of the contacts
are reported per structure to tell those apart.

Nothing is filtered by residue type or by loop membership: chemistry in the search would make
the composition of the result circular, and the loop overlap (91 % of the residues on the
2026-09-22 data) is a result, not an input. The residue type, the loop or core strand it sits
in, its RSA and its distance to the plane are all in binding_surface.json for analysis.

Output (written by pipeline.step_core_fold):
    07_core_fold/binding_surface.json  {struct_id: {"plane": {"normal", "offset"},
        "copper": {"source", "coord", "depth_below_plane", "shallowest_possible", "recessed"},
        "contacts_around_copper": 0-1, "n_residues": int, "settings": {...},
        "residues": [{"position", "residue_number", "insertion_code", "name", "distance",
                      "rsa", "region"}]}}
    and the "binding_surface" key of 07_core_fold/loop_regions.json, {region: ["start-end", ...]}
    in sequence positions, as the Django import and the app read it (residues outside every
    loop -- the N-His, core-strand edges -- go under "other").

Standalone (PyMOL sessions to check a structure by eye: the plane as a transparent square,
the binding-surface residues as orange sticks, "side" and "face" views):
    python lpmocat_pipeline/steps/binding_surface.py 8RRY 7PU1 4B5Q --pymol
    python lpmocat_pipeline/steps/binding_surface.py --all
"""
import os
import json
import math
import argparse

import numpy as np

from paths import DEFAULT_DATA_DIR

SLAB = 3.0              # a residue touches the plane when its outermost atom is this close (A)
MAX_CU_DEPTH = 3.0      # how far above the copper the plane may lie (A)
DEPTH_MARGIN = 1.0      # extra depth allowed for a recessed copper, over the shallowest possible
BALANCE = 1.0           # 0-1: penalty for contacts all on one side of the copper
ALLOW_ABOVE = 0         # residues allowed to stick through the plane
DIRECTIONS = 20000      # plane orientations in the first search

METALS = ("CU", "ZN", "NI", "CO", "MN", "FE", "MO", "MG")   # = step 3g's BRACE_METALS

SETTINGS = {"slab": SLAB, "max_cu_depth": MAX_CU_DEPTH, "depth_margin": DEPTH_MARGIN,
            "balance": BALANCE, "allow_above": ALLOW_ABOVE, "directions": DIRECTIONS}


# --- structure ------------------------------------------------------------------------
def read_structure(path):
    """Residues of a PDB file in file order, and the active-site metal.

    A residue is {"number", "icode", "name", "atoms": [names], "coords": (n, 3)} with its
    heavy atoms; the metal is {"name", "coord"} or None."""
    residues, metal, seen = [], None, {}
    with open(path) as handle:
        for line in handle:
            record = line[:6]
            if record not in ("ATOM  ", "HETATM"):
                continue
            element = line[76:78].strip() or line[12:16].strip()[0]
            if element in ("H", "D"):
                continue
            coord = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            if record == "HETATM":
                if metal is None and element.upper() in METALS:
                    metal = {"name": line[17:20].strip(), "coord": np.array(coord)}
                continue
            key = (int(line[22:26]), line[26])
            if key not in seen:
                seen[key] = len(residues)
                residues.append({"number": key[0], "icode": key[1].strip(),
                                 "name": line[17:20].strip(), "atoms": [], "coords": []})
            residues[seen[key]]["atoms"].append(line[12:16].strip())
            residues[seen[key]]["coords"].append(coord)
    for residue in residues:
        residue["coords"] = np.array(residue["coords"])
    return residues, metal


def catalytic_domain(residues, end_residue):
    """The residues up to the catalytic domain end (step 6), in file order."""
    if end_residue is None:
        return residues
    last = max((i for i, r in enumerate(residues) if r["number"] <= end_residue),
               default=len(residues) - 1)
    return residues[:last + 1]


# --- the plane ------------------------------------------------------------------------
def sphere_directions(count):
    """`count` roughly evenly spread unit vectors (Fibonacci sphere)."""
    index = np.arange(count) + 0.5
    phi = np.arccos(1 - 2 * index / count)
    theta = math.pi * (1 + 5 ** 0.5) * index
    return np.column_stack((np.cos(theta) * np.sin(phi),
                            np.sin(theta) * np.sin(phi),
                            np.cos(phi)))


def cone_directions(axis, count, max_angle):
    """`count` unit vectors within `max_angle` degrees of `axis` (deterministic)."""
    axis = axis / np.linalg.norm(axis)
    helper = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(helper, axis)) > 0.9:
        helper = np.array([0.0, 1.0, 0.0])
    u = np.cross(axis, helper)
    u /= np.linalg.norm(u)
    v = np.cross(axis, u)
    index = np.arange(count) + 0.5
    cos_max = math.cos(math.radians(max_angle))          # equal-area sampling of the cap
    cos_t = 1 - (index / count) * (1 - cos_max)
    sin_t = np.sqrt(np.clip(1 - cos_t ** 2, 0, 1))
    phi = math.pi * (1 + 5 ** 0.5) * index
    return (cos_t[:, None] * axis
            + (sin_t * np.cos(phi))[:, None] * u
            + (sin_t * np.sin(phi))[:, None] * v)


def residue_heights(coords, starts, directions, chunk=512):
    """(residues, directions) matrix of the largest projection of each residue's atoms."""
    out = np.empty((len(starts), len(directions)))
    for first in range(0, len(directions), chunk):
        block = directions[first:first + chunk]
        projection = coords @ block.T
        out[:, first:first + len(block)] = np.maximum.reduceat(projection, starts, axis=0)
    return out


def contact_spread(keep, offsets_from_cu, directions, chunk=256):
    """How evenly the touching residues lie around the copper, per direction: 1 = all around
    it, 0 = all on one side.

    Every touching residue is reduced to the direction of its offset from the copper within
    the plane; the mean of those unit vectors is short when they point every way and long
    when they crowd one side."""
    spread = np.ones(len(directions))
    for first in range(0, len(directions), chunk):
        block = directions[first:first + chunk]
        along = offsets_from_cu @ block.T                              # (residues, block)
        lateral = offsets_from_cu[:, None, :] - along[:, :, None] * block[None, :, :]
        length = np.linalg.norm(lateral, axis=2)
        unit = np.divide(lateral, np.where(length[:, :, None] > 1e-6, length[:, :, None], 1))
        mask = keep[:, first:first + len(block)]
        counts = mask.sum(axis=0)
        resultant = np.linalg.norm((unit * mask[:, :, None]).sum(axis=0), axis=1)
        spread[first:first + len(block)] = np.where(
            counts > 0, 1.0 - resultant / np.maximum(counts, 1), 1.0)
    return spread


def _score(heights, directions, cu, countable, centroids, settings):
    """Per direction: the plane offset, the copper's depth below it, the residues it touches
    and how evenly they lie around the copper."""
    allow_above = settings["allow_above"]
    if allow_above:
        offsets = np.partition(heights, -(allow_above + 1), axis=0)[-(allow_above + 1)]
    else:
        offsets = heights.max(axis=0)
    depth = offsets - directions @ cu
    touching = heights >= (offsets - settings["slab"])
    keep = touching & countable[:, None]
    counts = keep.sum(axis=0)
    # tie-break: residues sitting closer to the plane
    closeness = np.where(keep, settings["slab"] - (offsets - heights), 0.0).sum(axis=0)
    spread = (contact_spread(keep, centroids - cu, directions)
              if settings["balance"] else np.ones(len(directions)))
    return offsets, depth, counts, closeness, spread


def _best_direction(heights_of, cu, countable, centroids, settings, directions, cap):
    """The best (direction, offset, depth, count, heights, spread) of `directions` with the
    copper at most `cap` below the plane, plus the smallest copper depth any of them reach."""
    heights = heights_of(directions)
    offsets, depth, counts, closeness, spread = _score(
        heights, directions, cu, countable, centroids, settings)
    allowed = (depth >= -0.5) & (depth <= cap)
    if not allowed.any():
        return None, float(depth.min())
    balanced = counts * (1 - settings["balance"] * (1 - spread))
    score = np.where(allowed, balanced + closeness / (1e4 * max(len(countable), 1)), -1)
    winner = int(np.argmax(score))
    return ((directions[winner], float(offsets[winner]), float(depth[winner]),
             int(counts[winner]), heights[:, winner], float(spread[winner])),
            float(depth.min()))


def fit_plane(residues, cu, settings=None):
    """Fit the binding plane of one structure's catalytic domain. Returns
    {"normal", "offset", "cu_depth", "min_cu_depth", "recessed", "spread", "gaps", "picked",
     "above"} or None; `gaps` is each residue's distance below the plane and `picked` the
    indices of the binding-surface residues."""
    settings = {**SETTINGS, **(settings or {})}
    coords = np.vstack([r["coords"] for r in residues])
    starts = np.cumsum([0] + [len(r["coords"]) for r in residues])[:-1]
    centroids = np.array([r["coords"].mean(axis=0) for r in residues])
    countable = np.array([r["countable"] for r in residues])

    def heights_of(directions):
        return residue_heights(coords, starts, directions)

    sphere = sphere_directions(settings["directions"])
    cap = settings["max_cu_depth"]
    result, min_depth = _best_direction(heights_of, cu, countable, centroids, settings,
                                        sphere, cap)
    recessed = result is None
    if recessed:                       # a groove or a bowl: open the limit, flag the result
        cap = min_depth + settings["depth_margin"]
        result, _ = _best_direction(heights_of, cu, countable, centroids, settings, sphere, cap)
        if result is None:
            return None
    for angle, count in ((10.0, 2000), (3.0, 2000)):
        refined, _ = _best_direction(heights_of, cu, countable, centroids, settings,
                                     cone_directions(result[0], count, angle), cap)
        if refined is not None and refined[3] >= result[3]:
            result = refined
    normal, offset, depth, _, heights, spread = result

    gaps = offset - heights                    # 0 = touches the plane, < 0 = above it
    picked = [i for i in range(len(residues))
              if gaps[i] <= settings["slab"] and countable[i]]
    if picked:                                 # the spread of the plane that won
        lateral = centroids[picked] - cu
        lateral = lateral - np.outer(lateral @ normal, normal)
        length = np.linalg.norm(lateral, axis=1)
        unit = lateral / np.where(length[:, None] > 1e-6, length[:, None], 1)
        spread = float(1 - np.linalg.norm(unit.sum(axis=0)) / len(picked))
    return {"normal": normal, "offset": offset, "cu_depth": depth, "min_cu_depth": min_depth,
            "recessed": recessed, "spread": spread, "gaps": gaps, "picked": picked,
            "above": [i for i in range(len(residues)) if gaps[i] < -1e-6]}


# --- per structure --------------------------------------------------------------------
def binding_surface(pdb_path, end_residue=None, brace=None, sasa=None, labels=None,
                    settings=None):
    """The binding surface of one structure. `brace` is its histidine_brace.json entry (its
    "copper" is the copper position; files without it: the metal, else "pseudo_cu"), `sasa` its sasa.json entry, `labels` {position (1-based):
    loop or core strand}. Returns the binding_surface.json record, or None when the copper or
    the plane cannot be found."""
    settings = {**SETTINGS, **(settings or {})}
    brace, labels = brace or {}, labels or {}
    residues, metal = read_structure(pdb_path)
    domain = catalytic_domain(residues, end_residue)
    if not domain:
        return None

    if "copper" in brace:                                 # decided in step 6f
        if brace["copper"] is None:
            return None
        cu = np.array(brace["copper"], dtype=float)
        source = metal["name"] if brace.get("copper_source") == "metal" and metal else "pseudo-Cu"
    elif metal is not None:
        cu, source = metal["coord"], metal["name"]
    else:
        cu, source = brace.get("pseudo_cu"), "pseudo-Cu"
        if cu is None:
            return None
        cu = np.array(cu, dtype=float)

    rsa = {number: value for number, value in zip(sasa.get("residues", []),
                                                  sasa.get("rsa", []))} if sasa else {}
    for residue in domain:
        residue["rsa"] = rsa.get(residue["number"])
        # any solvent accessibility at all, no threshold; an unknown RSA counts as exposed
        residue["countable"] = residue["rsa"] is None or residue["rsa"] > 0

    plane = fit_plane(domain, cu, settings)
    if plane is None:
        return None
    picked = []
    for index in plane["picked"]:
        residue = domain[index]
        picked.append({"position": index + 1, "residue_number": residue["number"],
                       "insertion_code": residue["icode"], "name": residue["name"],
                       "distance": round(float(plane["gaps"][index]), 2),
                       "rsa": residue["rsa"], "region": labels.get(index + 1, "")})
    return {"plane": {"normal": [round(float(v), 4) for v in plane["normal"]],
                      "offset": round(plane["offset"], 3)},
            "copper": {"source": source, "coord": [round(float(v), 3) for v in cu],
                       "depth_below_plane": round(plane["cu_depth"], 2),
                       "shallowest_possible": round(plane["min_cu_depth"], 2),
                       "recessed": plane["recessed"]},
            "contacts_around_copper": round(plane["spread"], 2),
            "residues_above_plane": [f"{domain[i]['name']}{domain[i]['number']}"
                                     for i in plane["above"]],
            "n_residues": len(picked), "settings": settings, "residues": picked}


def region_ranges(record):
    """{region: ["start-end", ...]} of a binding_surface.json record: the runs of consecutive
    sequence positions per loop or core strand, everything else under "other". This is the
    shape the Django import and the app read from loop_regions.json."""
    by_region = {}
    for residue in record["residues"]:
        by_region.setdefault(residue["region"] or "other", []).append(residue["position"])
    ranges = {}
    for region, positions in by_region.items():
        runs, first, last = [], None, None
        for position in sorted(positions):
            if first is None:
                first = last = position
            elif position == last + 1:
                last = position
            else:
                runs.append(f"{first}-{last}")
                first = last = position
        if first is not None:
            runs.append(f"{first}-{last}")
        ranges[region] = runs
    return ranges


def assign_binding_surface(struct_ids, aligned_dir, boundaries, braces, sasa, labels=None,
                           settings=None):
    """7c: the binding surface of every structure. Returns
    ({struct_id: binding_surface.json record}, [struct_ids without one])."""
    records, failed = {}, []
    for struct_id in struct_ids:
        path = os.path.join(aligned_dir, f"{struct_id}.pdb")
        record = (binding_surface(path, boundaries.get(struct_id, {}).get("end_residue"),
                                  braces.get(struct_id), sasa.get(struct_id),
                                  (labels or {}).get(struct_id), settings)
                  if os.path.exists(path) else None)
        if record is None:
            failed.append(struct_id)
        else:
            records[struct_id] = record
    return records, failed


# --- standalone: PyMOL sessions --------------------------------------------------------
def write_pymol(struct_id, record, pdb_path, out_dir, size=24.0):
    """A .pml that loads the structure, draws the plane as a transparent square and shows the
    binding-surface residues as sticks, with "side" and "face" views."""
    normal = np.array(record["plane"]["normal"], dtype=float)
    cu = np.array(record["copper"]["coord"], dtype=float)
    point = cu + normal * record["copper"]["depth_below_plane"]
    helper = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(helper, normal)) > 0.9:
        helper = np.array([0.0, 1.0, 0.0])
    u = np.cross(normal, helper)
    u /= np.linalg.norm(u)
    v = np.cross(normal, u)
    corners = [point + size * (a * u + b * v)
               for a, b in ((-1, -1), (-1, 1), (1, -1), (1, 1))]     # triangle-strip order
    vertices = "".join(f"           VERTEX, {c[0]:.3f}, {c[1]:.3f}, {c[2]:.3f},\n"
                       for c in corners)
    selection = "+".join(str(r["residue_number"]) for r in record["residues"]) or "0"
    cu_block = ""
    if record["copper"]["source"] == "pseudo-Cu":
        cu_block = (f"pseudoatom pseudo_cu, pos=[{cu[0]}, {cu[1]}, {cu[2]}], color=orange\n"
                    "show spheres, pseudo_cu\nset sphere_scale, 0.4, pseudo_cu\n")

    def view_matrix(axes):
        """PyMOL set_view: the 3x3 part holds the camera axes as COLUMNS."""
        centre = point - normal * 12.0
        numbers = [value for row in np.array(axes).T for value in row]
        numbers += [0.0, 0.0, -120.0] + list(centre) + [60.0, 190.0, -20.0]
        return ", ".join(f"{value:.3f}" for value in numbers)

    path = os.path.join(out_dir, f"{struct_id}.pml")
    with open(path, "w") as handle:
        handle.write(f"""# Binding surface of {struct_id}: {record['n_residues']} residues touching the
# plane, copper {record['copper']['depth_below_plane']} A below it, contacts around it
# {record['contacts_around_copper']}{' (RECESSED: a groove, not a flat face)' if record['copper']['recessed'] else ''}
load {os.path.abspath(pdb_path)}, {struct_id}
hide everything
show cartoon, {struct_id}
color grey80, {struct_id}
show spheres, {struct_id} and not polymer
color orange, {struct_id} and not polymer
{cu_block}
select binding_surface, {struct_id} and resi {selection}
show sticks, binding_surface and not (name C+N+O)
color tv_orange, binding_surface
deselect
# label the binding-surface residues:
# label binding_surface and name CA, "%s%s" % (resn, resi)

python
from pymol.cgo import BEGIN, TRIANGLE_STRIP, COLOR, NORMAL, VERTEX, END
plane = [BEGIN, TRIANGLE_STRIP,
         COLOR, 0.25, 0.55, 0.90,
         NORMAL, {normal[0]:.4f}, {normal[1]:.4f}, {normal[2]:.4f},
{vertices}         END]
cmd.load_cgo(plane, "binding_plane")
cmd.set("cgo_transparency", 0.55, "binding_plane")
python end

set cartoon_transparency, 0.25, {struct_id}
bg_color white
set ray_shadows, 0

# "side" looks along the plane (is the copper face flat?), "face" looks onto it
set_view ({view_matrix((u, normal, np.cross(u, normal)))})
alias side, set_view ({view_matrix((u, normal, np.cross(u, normal)))})
alias face, set_view ({view_matrix((u, v, normal))})
""")
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("structures", nargs="*", help="struct_ids (default: all)")
    parser.add_argument("--all", action="store_true", help="every structure (the default)")
    parser.add_argument("--data-dir", default=DEFAULT_DATA_DIR)
    parser.add_argument("--out", help="where to write (default: <data-dir>/07_core_fold)")
    parser.add_argument("--pymol", action="store_true", help="also write a PyMOL session each")
    parser.add_argument("--slab", type=float, default=SLAB)
    parser.add_argument("--max-cu-depth", type=float, default=MAX_CU_DEPTH)
    parser.add_argument("--balance", type=float, default=BALANCE)
    parser.add_argument("--allow-above", type=int, default=ALLOW_ABOVE)
    args = parser.parse_args()

    data_dir = args.data_dir
    aligned_dir = os.path.join(data_dir, "04_alignment", "aligned")
    boundaries = json.load(open(os.path.join(data_dir, "06_catalytic_domain",
                                             "catalytic_domain_boundaries.json")))
    braces = json.load(open(os.path.join(data_dir, "06_catalytic_domain", "histidine_brace.json")))
    sasa_file = os.path.join(data_dir, "05_dssp_sasa", "sasa.json")
    sasa = json.load(open(sasa_file)) if os.path.exists(sasa_file) else {}
    loops_file = os.path.join(data_dir, "07_core_fold", "loop_regions.json")
    core_file = os.path.join(data_dir, "07_core_fold", "beta_core.json")
    labels = {}
    for path, key in ((loops_file, "loop_regions"), (core_file, "structures")):
        if not os.path.exists(path):
            continue
        content = json.load(open(path))
        content = content.get("structures", content) if key == "structures" else content
        for struct_id, entry in content.items():
            regions = entry if key == "structures" else entry.get(key, {})
            for name, value in regions.items():
                pieces = value.get("pieces", []) if isinstance(value, dict) else value
                label = name.replace("b", "β") if key == "structures" else name
                for piece in pieces:
                    start, end = (int(v) for v in piece.split("-"))
                    for position in range(start, end + 1):
                        labels.setdefault(struct_id, {}).setdefault(position, label)

    settings = {"slab": args.slab, "max_cu_depth": args.max_cu_depth, "balance": args.balance,
                "allow_above": args.allow_above}
    struct_ids = args.structures or sorted(boundaries)
    out_dir = args.out or os.path.join(data_dir, "07_core_fold")
    os.makedirs(out_dir, exist_ok=True)

    records, failed = assign_binding_surface(struct_ids, aligned_dir, boundaries, braces,
                                             sasa, labels, settings)
    for struct_id, record in records.items():
        copper = record["copper"]
        regions = {}
        for residue in record["residues"]:
            regions[residue["region"] or "other"] = regions.get(residue["region"] or "other",
                                                                0) + 1
        print(f"{struct_id:<16}{record['n_residues']:>3} residues, Cu "
              f"{copper['depth_below_plane']:>5} A below the plane"
              f"{' (RECESSED)' if copper['recessed'] else ''}, around it "
              f"{record['contacts_around_copper']}, "
              + " ".join(f"{k}:{v}" for k, v in sorted(regions.items())))
        if args.pymol:
            write_pymol(struct_id, record, os.path.join(aligned_dir, f"{struct_id}.pdb"),
                        out_dir)
    if failed:
        print(f"no binding surface for {len(failed)}: {', '.join(failed)}")
    with open(os.path.join(out_dir, "binding_surface.json"), "w") as handle:
        json.dump(records, handle, indent=1)
    print(f"\n{len(records)} structures -> {os.path.join(out_dir, 'binding_surface.json')}")


if __name__ == "__main__":
    main()
