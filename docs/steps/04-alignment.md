# Step 4 · Alignment

Step 4 superposes all prepared structures at once. The result gives every residue a column in
one alignment and puts all structures in one coordinate frame.

| Name for `--steps` | `align` or `4` |
| --- | --- |
| Code | `align_structures.py` |
| Reads | [Step 3](03-preparation.md): `preprocessed/<id>.pdb` · [step 2](02-data.md): `pdb_selection.csv` |
| Writes | `04_alignment/alignment.fasta`, `matrices.json`, `aligned/<id>.pdb` |
| External tools | US-align 20241108, bundled in `lpmocat_pipeline/tools/USalign/` |
| Run time | about 25 min for 160 structures |

## What it does

```text
USalign -dir preprocessed/ structure_list.txt -mm 4 -mol prot -m matrices.txt -o result
```

`-mm 4` is US-align's multiple structure alignment of monomeric chains: one alignment for all
structures, not a set of pairwise ones. `-mol prot` restricts it to protein atoms. The
alignment is read from US-align's output, the superposed structures are saved as
`aligned/<id>.pdb`, and the transformation matrices are stored for the web app's 3D viewer.

## Why structure only

LPMO families share only 2–19 % sequence identity with each other, while the
β-sandwich fold and the copper site are conserved. A sequence aligner cannot place such
sequences reliably; a structural alignment pairs the residues that occupy the same place in
the fold. The alignment is therefore **not refined by sequence** afterwards. Its gaps are
structural differences: a loop present in one family and absent in another gets its own
columns. This is also why the alignment is wide (thousands of columns for a few hundred
residues per protein), and why later steps find equivalent residues by column, e.g. the
Internal-His of a model at the column where metal-bound structures have theirs.

## Consequences for reruns

The alignment is global: adding, removing or changing one structure changes the columns of
all. Every step from 5 on reads the aligned files or the columns, so after step 4 all later
steps must be rerun. Files in `aligned/` of structures no longer in the alignment are removed.

The step fails if US-align returns fewer structures than it was given.

## Outputs

--8<-- "_formats/alignment.md"

## Known limits

- The last one or two residues of a chain are sometimes placed in a column far from the rest
  of the domain. [Step 6c](06-catalytic-domain.md#known-limits) is affected by this.
- Run time grows faster than the number of structures; this is the slowest step after
  CABS-flex.
