#!/usr/bin/env python
"""
Offline exploration (NOT the web app): what the binding surfaces of pipeline step 7c are made
of, and whether the patterns follow the family, the substrate or the oxidised position.

Reads 07_core_fold/binding_surface.json (the residues a plane tangent to the copper face touches),
05_dssp_sasa/sasa.json (the exposed residues they are compared against), the structural
alignment (which columns the residues fall in) and the live Django DB (family, substrates,
oxidised position, whether a paper reports them).

Sections:
  1 composition      residue types on the surface against the rest of the exposed catalytic
                     domain, overall and per family (the brace histidines are left out: the
                     plane is anchored at the copper, so they are there by construction)
  2 where            which loop or core strand the residues sit in, per family
  3 substrate        composition and loop use of the proteins of one substrate class
  4 position         C1 / C4 / C1+C4 within AA9 and AA10: size, copper depth, composition
  5 columns          alignment columns where many structures put a surface residue, per family
                     and shared between families -- the structurally equivalent positions
  6 aromatics        aromatic and polar counts per protein, by substrate class

Usage:
    python lpmocat_pipeline/analysis/binding_surface_patterns.py
    python lpmocat_pipeline/analysis/binding_surface_patterns.py --section columns --family AA10
"""
import os
import sys
import json
import argparse
import statistics as st
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve()
DATA_DIR = HERE.parents[2] / "data"
DJANGO_DIR = HERE.parents[2] / "lpmocat_web"
sys.path.insert(0, str(DJANGO_DIR))
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "config.settings")
import django  # noqa: E402

django.setup()
from catalogue.models import Protein  # noqa: E402

GROUPS = (("aromatic F/Y/W", "PHE TYR TRP"), ("His", "HIS"), ("charged D/E/K/R", "ASP GLU LYS ARG"),
          ("polar N/Q/S/T", "ASN GLN SER THR"), ("Gly/Pro", "GLY PRO"),
          ("aliphatic", "ALA VAL LEU ILE MET CYS"))
AMINO = ("TYR", "TRP", "PHE", "HIS", "ASN", "GLN", "SER", "THR", "ASP", "GLU",
         "LYS", "ARG", "GLY", "PRO", "ALA", "VAL", "LEU", "ILE", "MET", "CYS")


def load():
    """(surfaces, proteins) keyed by struct_id; proteins carry family, substrates and position."""
    surfaces = json.load(open(DATA_DIR / "07_core_fold" / "binding_surface.json"))
    sasa = json.load(open(DATA_DIR / "05_dssp_sasa" / "sasa.json"))
    braces = json.load(open(DATA_DIR / "06_catalytic_domain" / "histidine_brace.json"))
    boundaries = json.load(open(DATA_DIR / "06_catalytic_domain" / "catalytic_domain_boundaries.json"))
    alignment, name = {}, None
    for line in open(DATA_DIR / "04_alignment" / "alignment.fasta"):
        if line.startswith(">"):
            name = line[1:].strip().split()[0].split(".pdb")[0]   # ">7OKR.pdb:A\tL=192 ..."
            alignment[name] = ""
        elif name:
            alignment[name] += line.strip()

    proteins = {}
    for protein in Protein.objects.prefetch_related("substrate", "regioselectivity").all():
        struct = protein.selected_structure
        if struct in surfaces:
            proteins[struct] = {
                "code": protein.protein_code,
                "family": boundaries.get(struct, {}).get("family", ""),
                "substrates": {s.name for s in protein.substrate.all()},
                "position": " and ".join(sorted(r.name.replace("-hydroxylating", "")
                                                 .replace("-dehydrogenating", "")
                                                for r in protein.regioselectivity.all())),
            }
    return surfaces, proteins, sasa, braces, boundaries, alignment


def background(struct_id, surfaces, sasa, braces, boundaries):
    """(surface residues, other exposed residues of the catalytic domain) as name lists, with
    the brace histidines left out of both."""
    entry = sasa.get(struct_id)
    if not entry:
        return [], []
    brace = {braces.get(struct_id, {}).get("his1"), braces.get(struct_id, {}).get("his2")}
    end = boundaries.get(struct_id, {}).get("end_residue")
    picked = {residue["residue_number"] for residue in surfaces[struct_id]["residues"]}
    on, off = [], []
    for number, name, rsa in zip(entry["residues"], entry["names"], entry["rsa"]):
        if (end is not None and number > end) or rsa <= 0 or number in brace:
            continue
        (on if number in picked else off).append(name)
    return on, off


def enrichment(on, off, rows=GROUPS):
    """[(label, % on the surface, % elsewhere, enrichment)]."""
    total_on, total_off = len(on) or 1, len(off) or 1
    counts_on, counts_off = Counter(on), Counter(off)
    out = []
    for label, names in rows:
        share_on = sum(counts_on[n] for n in names.split()) / total_on
        share_off = sum(counts_off[n] for n in names.split()) / total_off
        out.append((label, 100 * share_on, 100 * share_off,
                    share_on / share_off if share_off else float("inf")))
    return out


def section_composition(data, args):
    surfaces, proteins, sasa, braces, boundaries, _ = data
    print("1 COMPOSITION  what the binding surfaces are made of\n")
    everything = defaultdict(lambda: ([], []))
    for struct_id in surfaces:
        on, off = background(struct_id, surfaces, sasa, braces, boundaries)
        family = proteins.get(struct_id, {}).get("family", "?")
        for key in ("all", family):
            everything[key][0].extend(on)
            everything[key][1].extend(off)
    print(f"  {'group':<17}{'surface':>9}{'elsewhere':>11}{'enrichment':>12}")
    for label, share_on, share_off, ratio in enrichment(*everything["all"]):
        print(f"  {label:<17}{share_on:>8.1f}%{share_off:>10.1f}%{ratio:>11.2f}x")
    print("\n  per residue type (enrichment, share of the surface):")
    on, off = everything["all"]
    values = sorted(enrichment(on, off, [(n, n) for n in AMINO]), key=lambda r: -r[3])
    for index in range(0, len(values), 4):
        print("   " + "  ".join(f"{label[:3]} {ratio:>5.2f}x ({share:>4.1f}%)"
                                for label, share, _, ratio in values[index:index + 4]))
    print("\n  enrichment per family (surface share in brackets):")
    families = sorted(k for k in everything if k not in ("all", "?"))
    print(f"  {'group':<17}" + "".join(f"{f:>16}" for f in families))
    for label, names in GROUPS:
        cells = []
        for family in families:
            row = next(r for r in enrichment(*everything[family], rows=[(label, names)]))
            cells.append(f"{row[3]:>6.2f}x ({row[1]:>4.1f}%)")
        print(f"  {label:<17}" + "".join(f"{c:>16}" for c in cells))


def section_where(data, args):
    surfaces, proteins, *_ = data
    print("\n2 WHERE  which loop or core strand the residues sit in\n")
    per_family = defaultdict(Counter)
    sizes = defaultdict(list)
    for struct_id, record in surfaces.items():
        family = proteins.get(struct_id, {}).get("family", "?")
        sizes[family].append(record["n_residues"])
        for residue in record["residues"]:
            per_family[family][residue["region"] or "other"] += 1
    regions = ["L2", "L3", "LS", "L8", "LC", "other"]
    print(f"  {'family':<8}{'n':>4}{'residues':>10}   " + "".join(f"{r:>8}" for r in regions))
    for family in sorted(per_family, key=lambda f: -len(sizes[f])):
        total = sum(per_family[family].values()) or 1
        cells = "".join(f"{100 * per_family[family][r] // total:>7}%" for r in regions)
        print(f"  {family:<8}{len(sizes[family]):>4}{st.median(sizes[family]):>7.0f} med"
              f"   {cells}")


def _substrate_class(entry):
    """One substrate class per protein, for grouping: the specific ones win over cellulose."""
    substrates = entry["substrates"]
    for name in ("Starch", "Pectin", "Xylan", "Chitin"):
        if name in substrates:
            return name
    if substrates & {"Xyloglucan", "Mixed-linkage glucan", "Glucomannan"}:
        return "Cellulose + hemicellulose"
    if "Cellulose" in substrates or "Cello-oligosaccharides" in substrates:
        return "Cellulose only"
    return None


def section_substrate(data, args):
    surfaces, proteins, sasa, braces, boundaries, _ = data
    print("\n3 SUBSTRATE  composition and loop use by substrate class\n")
    groups = defaultdict(lambda: ([], []))
    regions = defaultdict(Counter)
    members = defaultdict(list)
    for struct_id in surfaces:
        entry = proteins.get(struct_id)
        label = _substrate_class(entry) if entry else None
        if not label:
            continue
        on, off = background(struct_id, surfaces, sasa, braces, boundaries)
        groups[label][0].extend(on)
        groups[label][1].extend(off)
        members[label].append(struct_id)
        for residue in surfaces[struct_id]["residues"]:
            regions[label][residue["region"] or "other"] += 1
    order = sorted(groups, key=lambda k: -len(members[k]))
    print(f"  {'substrate':<26}{'n':>4}{'residues':>10}" + "".join(f"{g[:9]:>12}" for g, _ in GROUPS))
    for label in order:
        sizes = [surfaces[s]["n_residues"] for s in members[label]]
        cells = "".join(f"{row[3]:>10.2f}x" for row in enrichment(*groups[label]))
        print(f"  {label:<26}{len(members[label]):>4}{st.median(sizes):>7.0f} med{cells}")
    print(f"\n  loops:  {'substrate':<24}" + "".join(f"{r:>8}" for r in
                                                     ("L2", "L3", "LS", "L8", "LC", "other")))
    for label in order:
        total = sum(regions[label].values()) or 1
        cells = "".join(f"{100 * regions[label][r] // total:>7}%"
                        for r in ("L2", "L3", "LS", "L8", "LC", "other"))
        print(f"          {label:<24}{cells}")


def section_position(data, args):
    surfaces, proteins, sasa, braces, boundaries, _ = data
    print("\n4 OXIDISED POSITION  within a family\n")
    for family in ("AA9", "AA10"):
        groups = defaultdict(list)
        for struct_id, record in surfaces.items():
            entry = proteins.get(struct_id)
            if not entry or entry["family"] != family or not entry["position"]:
                continue
            groups[entry["position"]].append(struct_id)
        print(f"  {family}:")
        print(f"    {'position':<12}{'n':>4}{'residues':>10}{'Cu depth':>10}{'spread':>8}"
              + "".join(f"{g[:9]:>11}" for g, _ in GROUPS))
        for position, members in sorted(groups.items(), key=lambda kv: -len(kv[1])):
            on, off = [], []
            for struct_id in members:
                a, b = background(struct_id, surfaces, sasa, braces, boundaries)
                on.extend(a)
                off.extend(b)
            sizes = [surfaces[s]["n_residues"] for s in members]
            depths = [surfaces[s]["copper"]["depth_below_plane"] for s in members]
            spread = [surfaces[s]["contacts_around_copper"] for s in members]
            cells = "".join(f"{row[3]:>10.2f}x" for row in enrichment(on, off))
            print(f"    {position:<12}{len(members):>4}{st.median(sizes):>7.0f} med"
                  f"{st.median(depths):>9.1f}A{st.median(spread):>8.2f}{cells}")


def column_roles(surfaces, braces, alignment, data_dir=DATA_DIR):
    """{column: role} for the active-site positions the pipeline knows: the two brace
    histidines (step 6f), the second-sphere aromatic and the motif residues (step 8)."""
    interactions = json.load(open(data_dir / "08_second_sphere" / "brace_interactions.json"))
    votes = defaultdict(Counter)
    for struct_id in surfaces:
        sequence = alignment.get(struct_id)
        if not sequence:
            continue
        number_to_column, position = {}, 0
        numbers = [residue["residue_number"] for residue in surfaces[struct_id]["residues"]]
        # every residue of the structure, not only the surface ones
        record = interactions.get(struct_id, {})
        marks = {}
        brace = braces.get(struct_id, {})
        if brace.get("his1_column"):
            votes[brace["his1_column"]]["N-His"] += 1
        if brace.get("his2_column"):
            votes[brace["his2_column"]]["internal His"] += 1
        aromatic = record.get("second_sphere_aromatic") or {}
        if aromatic.get("column"):
            votes[aromatic["column"]]["second-sphere aromatic"] += 1
        motif = record.get("second_sphere_motif") or {}
        for residue in motif.get("residues", []):
            marks[residue["residue"]] = f"motif {residue['role']}"
        distal = record.get("distal_axial") or {}
        if distal.get("residue"):
            marks[distal["residue"]] = "shielding residue"
        if marks:
            column, position = 0, 0
            for index, letter in enumerate(sequence, start=1):
                if letter == "-":
                    continue
                position += 1
            # map residue numbers through the structure's own order
            order = json.load(open(data_dir / "05_dssp_sasa" / "sasa.json")).get(struct_id, {})
            numbers_in_order = order.get("residues", [])
            position_of = {number: i + 1 for i, number in enumerate(numbers_in_order)}
            column_of, position = {}, 0
            for index, letter in enumerate(sequence, start=1):
                if letter != "-":
                    position += 1
                    column_of[position] = index
            for number, role in marks.items():
                where = column_of.get(position_of.get(number, 0))
                if where:
                    votes[where][role] += 1
    return {column: counter.most_common(1)[0][0] for column, counter in votes.items()
            if counter.most_common(1)[0][1] >= 5}


def section_columns(data, args):
    surfaces, proteins, sasa, braces, boundaries, alignment = data
    roles = column_roles(surfaces, braces, alignment)
    print("\n5 COLUMNS  alignment columns many structures put a surface residue in")
    print("   (a role in brackets = what the pipeline calls that position)\n")
    columns = defaultdict(lambda: defaultdict(list))     # family -> column -> [residue names]
    members = Counter()
    for struct_id, record in surfaces.items():
        family = proteins.get(struct_id, {}).get("family", "?")
        sequence = alignment.get(struct_id)
        if not sequence:
            continue
        members[family] += 1
        position_to_column, position = {}, 0
        for index, letter in enumerate(sequence, start=1):
            if letter != "-":
                position += 1
                position_to_column[position] = index
        for residue in record["residues"]:
            column = position_to_column.get(residue["position"])
            if column:
                columns[family][column].append(residue["name"])
    families = [f for f in sorted(members) if members[f] >= (2 if args.family else 3)]
    for family in families:
        if args.family and family != args.family:
            continue
        total = members[family]
        common = sorted(((len(names), column, names) for column, names in columns[family].items()
                         if len(names) >= 0.5 * total), reverse=True)
        print(f"  {family} ({total} structures): {len(common)} columns in at least half of them")
        for count, column, names in common[:12]:
            top = Counter(names).most_common(3)
            role = f"  [{roles[column]}]" if column in roles else ""
            print(f"    column {column:<5} {100 * count // total:>3}%  "
                  + ", ".join(f"{name.capitalize()} {100 * n // count}%" for name, n in top) + role)
    print("\n  columns shared between families (at least half the structures of each):")
    shared = defaultdict(list)
    for family in families:
        total = members[family]
        for column, names in columns[family].items():
            if len(names) >= 0.5 * total:
                shared[column].append(family)
    for column, group in sorted(shared.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        if len(group) >= 3:
            names = Counter(n for family in group for n in columns[family][column])
            role = f"  [{roles[column]}]" if column in roles else ""
            print(f"    column {column:<5} {', '.join(group):<40} "
                  + ", ".join(f"{name.capitalize()} {count}" for name, count in names.most_common(4))
                  + role)


def section_aromatics(data, args):
    surfaces, proteins, sasa, braces, boundaries, _ = data
    print("\n6 AROMATICS  per protein, by substrate class\n")
    per_class = defaultdict(list)
    for struct_id, record in surfaces.items():
        entry = proteins.get(struct_id)
        label = _substrate_class(entry) if entry else None
        if not label:
            continue
        names = [residue["name"] for residue in record["residues"]]
        per_class[label].append((sum(1 for n in names if n in ("TYR", "TRP", "PHE")),
                                 sum(1 for n in names if n == "TYR"),
                                 sum(1 for n in names if n in ("ASN", "GLN")),
                                 len(names), entry["code"]))
    print(f"  {'substrate':<26}{'n':>4}{'aromatics':>12}{'Tyr':>8}{'Asn+Gln':>10}"
          f"{'none at all':>13}")
    for label, rows in sorted(per_class.items(), key=lambda kv: -len(kv[1])):
        aromatics = [r[0] for r in rows]
        print(f"  {label:<26}{len(rows):>4}{st.median(aromatics):>9.1f} med"
              f"{st.median([r[1] for r in rows]):>8.1f}{st.median([r[2] for r in rows]):>10.1f}"
              f"{sum(1 for a in aromatics if a == 0):>13}")


SECTIONS = {"composition": section_composition, "where": section_where,
            "substrate": section_substrate, "position": section_position,
            "columns": section_columns, "aromatics": section_aromatics}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--section", choices=sorted(SECTIONS), action="append")
    parser.add_argument("--family", help="only this family, where a section splits by family")
    args = parser.parse_args()
    data = load()
    surfaces = data[0]
    print(f"{len(surfaces)} structures with a binding surface, "
          f"{sum(r['n_residues'] for r in surfaces.values())} residues in total\n")
    for name in (args.section or ["composition", "where", "substrate", "position", "columns",
                                  "aromatics"]):
        SECTIONS[name](data, args)


if __name__ == "__main__":
    main()
