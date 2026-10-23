# LPMOcat data pipeline

LPMOcat is a catalogue of lytic polysaccharide monooxygenases (LPMOs, CAZy families AA9–AA11
and AA13–AA17). This documentation describes the pipeline that builds its data: it collects
every characterized or structurally known LPMO, gives each one structure, aligns all
structures, and computes the same annotations and properties for each. The web app only
loads and displays the results; the one property it adds itself is hydrophobicity, a fixed
value per residue type on each of seven published scales (see
[Hydrophobicity scales](reference/hydrophobicity-scales.md)).

It is written for readers who want to understand how a value in the catalogue was obtained,
rerun the pipeline, or extend it.

## The 12 steps

The pipeline is one script, `lpmocat_pipeline/pipeline.py`. Each step writes into a folder
named after it and reads only files of earlier steps, so any step can be rerun on its own.

| | Step | Folder | What it produces |
| --- | --- | --- | --- |
| 1 | [Proteins](steps/01-proteins.md) | `01_proteins` | The list of LPMOs: CAZy, plus characterized LPMOs from Swiss-Prot, BRENDA and curation |
| 2 | [Data](steps/02-data.md) | `02_data` | One structure per protein (PDB, AlphaFold DB, CAZyme3D or an own model); UniProt accessions, domains, EC numbers, taxonomy, reported activities |
| 3 | [Preparation](steps/03-preparation.md) | `03_preparation` | One clean chain per structure, starting at the N-His, with its active-site metal only |
| 4 | [Alignment](steps/04-alignment.md) | `04_alignment` | A multiple structure alignment and the superposed structures |
| 5 | [DSSP and SASA](steps/05-dssp-sasa.md) | `05_dssp_sasa` | Secondary structure and solvent accessibility per residue |
| 6 | [Catalytic domain](steps/06-catalytic-domain.md) | `06_catalytic_domain` | Disulfide bonds, the catalytic domain end, the histidine brace and the pseudo-copper |
| 7 | [Core fold](steps/07-core-fold.md) | `07_core_fold` | β-sandwich core strands, loops, binding surface |
| 8 | [Second sphere](steps/08-second-sphere.md) | `08_second_sphere` | Axial aromatic, brace contacts, second-sphere motif, shielding residue, brace motifs |
| 9 | [Distance matrix](steps/09-distance-matrix.md) | `09_distance_matrix` | Pairwise distances between catalytic domains, for the similarity tree |
| 10 | [Flexibility](steps/10-flexibility.md) | `10_flexibility` | B-factors or pLDDT, GNM fluctuations, CABS-flex RMSF |
| 11 | [Electrostatics](steps/11-electrostatics.md) | `11_electrostatics` | Atomic partial charges and radii |
| 12 | [Hole hopping](steps/12-hole-hopping.md) | `12_hole_hopping` | Ranked charge-transfer paths between the surface and the copper |

Steps 1–4 collect, prepare and align the structures; steps 5–8 annotate them; steps 9–12
compute their properties. After the pipeline, the [Django import](django-import.md) loads the
files into the catalogue's database.

## What goes in and what comes out

**In**: public databases (CAZy, RCSB PDB, AlphaFold DB, CAZyme3D, UniProt and UniParc,
InterPro, NCBI Taxonomy, Europe PMC), a downloaded BRENDA release, and the
[curation files](reference/curation.md) kept with the code.

**Out**: per structure, a prepared and a superposed PDB file, and JSON or CSV files with
per-residue and per-protein results. All of them are described under
[Output formats](reference/outputs.md).

## Principles that hold for every step

- **One method for all proteins.** Every annotation and property is computed the same way
  for experimental structures and models, so values are comparable across the catalogue.
- **Structural alignment only.** Residues are compared through one multiple structure
  alignment; it is never refined by sequence.
- **Structures are trimmed only at the N-terminus**, to the N-His. The catalytic domain end
  is stored and applied where needed.
- **No waters or ligands** in any calculation; the active-site metal is the only non-protein
  atom kept.
- **Rules before exceptions.** Structures are handled by general rules; the few exceptions
  are explicit rows in the curation files, each with its reason.
- **Sources are kept apart.** Where databases report activities, every claim keeps its
  source and paper; nothing is merged silently.

## Where to start

- To run it: [Installation](installation.md), then [Running the pipeline](running.md).
- To read the output files: [Conventions](conventions.md), then
  [Output formats](reference/outputs.md).
- To understand one result: the page of the [step](steps/index.md) that produced it.
