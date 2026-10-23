# Conventions

## Structure identifiers

Every protein of the catalogue has one structure, named by its `struct_id`
(`selected_pdb` in `pdb_selection.csv`). It keys every output file.

| Form | Source | Example |
| --- | --- | --- |
| PDB code | RCSB PDB | `8RRY` |
| `AF-<UniProt accession>` | AlphaFold DB | `AF-G2RB72` |
| `C3D-<GenBank accession>` (the version dot written as `_`) | CAZyme3D | `C3D-BDB32680_1` |
| `AF3-<protein code>` | An AlphaFold 3 model of an in-house database entry | |

## Residue numbers, sequence positions and alignment columns { #numbering }

A residue is addressed in three ways. Which one a file uses is stated in its description;
mixing them up is the most common mistake when reading the files.

| | What it is | Used by |
| --- | --- | --- |
| **Residue number** | The number in the structure file: the author numbering of a PDB entry, the UniProt position in an AlphaFold DB model. It rarely starts at 1 and can have gaps | `sasa.json`, `dssp.json`, `disulfides.json`, `catalytic_domain_boundaries.json`, `histidine_brace.json`, `brace_interactions.json`, `flexibility_raw.json`, `hole_hopping.json`, PQR files |
| **Sequence position** | The 1-based index of the residue in the prepared structure, in file order. Position 1 is the first residue, normally the N-His | `beta_core.json`, `loop_regions.json`, `position` in `binding_surface.json` |
| **Alignment column** | The 1-based column of the residue in `alignment.fasta` | `his1_column`, `his2_column`, `column`, `columns` |

To convert:

- position → residue number: the *n*-th item of `residues` in `sasa.json`;
- position → column: the *n*-th non-gap letter of the structure's row in `alignment.fasta`.

Some PDB entries use insertion codes (a residue `79E` after `79`). The files address
residues by number alone, so the residues sharing a number are not distinguished except in
`binding_surface.json`, which stores the insertion code.

## N-His and Internal-His { #brace-histidines }

The copper of an LPMO is held by two histidines, the histidine brace. This documentation and
the web app call them the **N-His** (the N-terminal histidine, which binds through its amino
group and a ring nitrogen) and the **Internal-His**. They are not named by sequence
position, because most structures are not numbered from 1. In file keys, code and log lines
they appear as `his1` / `His-1` and `his2` / `His-2`.

## Copper and pseudo-copper

Step 6f decides the [copper position](steps/06-catalytic-domain.md#copper-position) once
for every structure, and all later steps use it: the active-site metal when the rings of
both brace histidines hold it, otherwise the
[pseudo-copper](steps/06-catalytic-domain.md#pseudo-copper), the centroid of the brace
nitrogens.

## Experimental structures and models

`source` in `pdb_selection.csv` is `PDB` for experimental structures and `AlphaFold`,
`CAZyme3D` or `Local model` for models. For models, the B-factor column holds pLDDT, the
N-terminus may be trimmed by a SignalP prediction ([3h](steps/03-preparation.md)), and there
is never a metal unless the model was built with one.

## Sub-step labels

Sub-steps are labelled with the step number and a letter (2a, 3h, 6e). The same labels
appear in the code comments and in the log.

## Loop and strand names

Core strands are `b1`–`b9` in files and β1–β9 in text. Loops are L2, L3, LS, L8 and LC; see
[step 7](steps/07-core-fold.md).
