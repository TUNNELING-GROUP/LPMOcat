# Step 8 · Second sphere

Step 8 describes the surroundings of the copper site: the axial aromatic residue, the
contacts of the two brace histidines, the second-sphere hydrogen-bonding motif, the residue
shielding the copper from the solvent, and the sequence motifs that follow the brace
histidines.

| Name for `--steps` | `second_sphere` or `8` |
| --- | --- |
| Code | `brace_interactions.py`, `pipeline.py` `step_second_sphere` |
| Reads | [Step 4](04-alignment.md): alignment, aligned structures · [step 6](06-catalytic-domain.md): `histidine_brace.json`, `catalytic_domain_boundaries.json` (families, for the log) |
| Writes | `08_second_sphere/brace_interactions.json` |
| Run time | a few seconds |

## 8a · Axial aromatic

**With a metal**: the Phe or Tyr with any atom closest to the active-site metal.
**Without**: the Phe or Tyr at the alignment columns where the metal-bound structures have
theirs, the most frequent column first, then columns up to 3 away (`BRACE_COLUMN_WINDOW`),
nearest first. This is the same approach as for the brace in
[6f](06-catalytic-domain.md#6f-histidine-brace).

## 8b · Contacts of the brace histidines

For every structure with a complete brace:

| Contact | Criterion |
| --- | --- |
| Hydrogen bond | A ring nitrogen (ND1, NE2) of a brace histidine within **3.5 Å** of an acceptor atom of another residue (ND1, NE2, OD1, OD2, OE1, OE2, OG, OG1, OH) or of a donor atom (ND1, NE2, OG, OG1, OH, backbone N) |
| π–π stacking | Ring centroid of a brace histidine within **6.5 Å** of the ring centroid of another histidine, the ring planes less than 80° apart |
| Cation–π | Lys NZ, or the centre of Arg NE, NH1, NH2, within **6.0 Å** of a brace ring centroid and lying over the ring (the direction to it makes \|cos\| > 0.5 with the ring normal) |

The hydrogen-bond test is by distance between heavy atoms only; no hydrogen positions or
angles are used. Apart from the backbone N as a donor, partners are side-chain atoms.

## The copper position for 8c and 8d

Both use the [copper position](06-catalytic-domain.md#copper-position) of step 6f: the
metal when the rings of both brace histidines hold it, otherwise the pseudo-copper. The
output records which (`copper`).

## 8c · Second-sphere motif

The motif is read backwards from the axial aromatic (8a) and written with one-letter codes
and the spacing between its residues, e.g. `H-X8-Q-X-Y`.

**If the residue two before the aromatic is Gln or Glu**, the motif's first residue is the
one that holds it, searched up to 16 residues before it:

1. a histidine π-stacked on the Internal-His (8b);
2. else the residue required by a published sequence motif (`KNOWN_MOTIFS`), which is used
   before the contact search because that residue is not in contact in every structure;
3. else the residue with a side-chain N or O closest to the Gln/Glu side chain, within 3.5 Å.

Without any of these the motif is just `Q-X-Y` or `E-X-F`.

**Otherwise** (the chitin-active AA10s, which have Asn there), the motif is headed by the
Gln or Glu whose side-chain N or O is closest to the copper, within 8.0 Å. It lies elsewhere
in the sequence, so the motif is written with an ellipsis: `E…N-X-F`.

## 8d · Shielding residue

The residue covering the copper's axial position on the solvent side, opposite the axial
aromatic (the conserved Ala of AA10s). The reference direction is the normal of the plane
through the N-His N, N-His ND1 and Internal-His NE2, pointing away from the aromatic. The
shielding residue is the one with a heavy side-chain atom (for Gly, its CA) closest to the
copper among atoms within **4.5 Å** of the copper and within **65°** of that normal. If no
atom qualifies, the position is open (`{"name": null}`).

The cone is 65° wide because the Ala Cβ of some AA10s lies 53–61° off the normal, while atoms
in the equatorial plane lie near 90°.

## 8e · Motifs after the brace histidines

The residues following each brace histidine are compared with a list of motifs; the first
that fits is taken.

| Histidine | Motifs, in order | Sequence stored |
| --- | --- | --- |
| N-His | `H-G-X4-P-X2-R`, `H-G-X4-P`, `H-G` | 10 residues from the histidine |
| Internal-His | `H-X-G-P`, `H-X2-G-P`, `H-X-G` | 5 residues from the histidine |

A position after a chain break (peptide C–N distance above 2.0 Å) counts as missing. When no
motif fits, `motif` is null and the sequence is still stored, so the variants can be
compared.

## Parameters

Constants of `brace_interactions.py` unless noted.

| Constant | Value | Meaning |
| --- | --- | --- |
| `BRACE_COLUMN_WINDOW` (`pipeline.py`) | `3` columns | Columns searched next to the known aromatic columns (8a) |
| `CUTOFF` (in `detect_interactions`) | `3.5` Å | Hydrogen-bond distance |
| (in `detect_interactions`) | `6.5` Å, 80° | π–π stacking: centroid distance, largest angle between ring planes |
| (in `detect_interactions`) | `6.0` Å, `0.5` | Cation–π: distance, smallest \|cos\| to the ring normal |
| `MOTIF_WINDOW` | `16` residues | How far before the Gln/Glu the first motif residue is searched |
| `MOTIF_HBOND` | `3.5` Å | Side-chain N/O distance of a first residue found by contact |
| `MOTIF_GLN_GLU_TO_CU` | `8.0` Å | Largest copper distance of the Gln/Glu heading an `E…N-X-F` motif |
| `KNOWN_MOTIFS` | R-X4-E-X-F, R-X5-E-X-F, N-X-E-X-Y, Q-X2-Q-X-Y, Q-X8-Q-X-Y, H-X5…8-Q-X-Y | Published sequence motifs (Vu et al. 2019 and the works cited there) |
| `DISTAL_MAX_DISTANCE` | `4.5` Å | Largest copper distance of a shielding atom |
| `DISTAL_MAX_ANGLE` | `65`° | Largest angle from the brace-plane normal |
| `N_HIS_WINDOW`, `INTERNAL_HIS_WINDOW` | `10`, `5` residues | Sequence stored after each histidine |
| `PEPTIDE_BOND_MAX` | `2.0` Å | C–N distance above which the chain is broken |

## Output

### brace_interactions.json

--8<-- "_formats/brace_interactions.md"

## What the log tells you

How many aromatics came from the metal and from the alignment, with their columns; the
number of contacts of each kind; in how many structures the metal or the pseudo-copper was
used; and per family the counts of each second-sphere motif, shielding residue and brace
motif.

## Known limits

- π–π stacking is tested against other histidines only, not against Phe, Tyr or Trp.
- In a model that folds the active site wrongly, the motif and the shielding residue are
  wrong with it; the brace nitrogen distance of [6f](06-catalytic-domain.md) is the check.
- The output keeps the name `distal_axial` for the shielding residue.
