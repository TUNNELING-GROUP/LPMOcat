# Step 9 · Distance matrix

Step 9 computes the pairwise distances between the catalytic domains, from which the web app
draws its similarity tree.

| Name for `--steps` | `distance_matrix` or `9` |
| --- | --- |
| Code | `pipeline.py` `step_distance_matrix` |
| Reads | [Step 4](04-alignment.md): alignment, aligned structures · [step 6](06-catalytic-domain.md): `catalytic_domain_boundaries.json` |
| Writes | `09_distance_matrix/distance_matrix.csv` |
| External tools | Biopython |
| Run time | about 2 min |

## What it does

1. Each row of the structural alignment is **masked after the catalytic domain end**: the
   residues past it become gaps.
2. Biopython's `DistanceCalculator("blosum62")` gives the distance of every pair: one minus
   the pair's summed BLOSUM62 score divided by the larger of the two sequences' self-scores,
   over the columns where neither has a gap.

The residues compared are those the structural alignment pairs, so the distance is a
sequence distance over structurally equivalent positions.

## Why

- **Catalytic domains only.** A CBM or linker present in one protein and absent in another
  would otherwise pull two similar catalytic domains apart.
- **Gap columns skipped per pair.** The alignment is very gappy by construction (see
  [step 4](04-alignment.md)); counting gaps as differences would measure loop presence, not
  sequence similarity.
- **BLOSUM62 and not identity.** At the low identity between LPMO families a substitution
  matrix still separates conservative from non-conservative replacements.

The web app clusters the matrix by average linkage. The tree groups proteins by sequence
similarity of their catalytic domains; it is not a phylogeny.

## Output

### distance_matrix.csv

--8<-- "_formats/distance_matrix.md"

## Known limits

- A structure without a catalytic domain end is compared over its whole chain; the log names
  it.
- If a structure's residue count differs from its alignment row, its residues are numbered
  from 1 for the masking, with a `WARNING`.
