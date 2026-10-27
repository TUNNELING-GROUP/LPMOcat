# Output formats

Every file the pipeline writes, in step order. Paths are relative to the data folder.
Read [Conventions](../conventions.md) first: the files use residue numbers, sequence
positions or alignment columns, and it matters which.

| File | Step | Content |
| --- | --- | --- |
| [`01_proteins/cazy_lpmos.csv`](#cazy_lpmoscsv) | 1a | LPMOs listed by CAZy |
| [`01_proteins/other_lpmos.csv`](#other_lpmoscsv) | 1b–1e | Characterized LPMOs from other lists |
| `02_data/<id>.pdb`, `.cif` | 2a–2d | Structure files as downloaded |
| [`02_data/pdb_selection.csv`](#pdb_selectioncsv) | 2 | The structure chosen for each protein |
| [`02_data/uniprot_ids.json`](#uniprot_idsjson) | 2f | UniProt accessions per protein |
| [`02_data/uniprot_references.json`](#uniprot_referencesjson) | 2f | Papers cited by those entries |
| [`02_data/domain_annotations.json`](#domain_annotationsjson) | 2g | Domains, signal peptides |
| [`02_data/taxonomy.json`](#taxonomyjson) | 2i | NCBI lineages |
| [`02_data/activities.json`](#activitiesjson) | 2j | Reported substrates and oxidised positions |
| [`03_preparation/filled/modelled_residues.json`](#modelled_residuesjson) | 3a | Rebuilt loops |
| [`03_preparation/signal_peptides.json`](#signal_peptidesjson) | 3h | SignalP predictions |
| [`03_preparation/modified_residues.json`](#modified_residuesjson) | 3f | Modified residues of the deposited chains |
| [`03_preparation/preprocessed/<id>.pdb`](#prepared-structures) | 3 | Prepared structures |
| [`04_alignment/…`](#alignment-files) | 4 | Alignment, matrices, superposed structures |
| [`05_dssp_sasa/sasa.json`](#sasajson) | 5a | Solvent accessibility |
| [`05_dssp_sasa/secondary_structure.json`](#secondary_structurejson) | 5b | Three-state secondary structure |
| [`05_dssp_sasa/dssp.json`](#dsspjson) | 5b | Full DSSP result |
| [`06_catalytic_domain/disulfides.json`](#disulfidesjson) | 6a | Disulfide bonds |
| [`06_catalytic_domain/catalytic_domain_boundaries.json`](#catalytic_domain_boundariesjson) | 6b–6e | Catalytic domain end |
| [`06_catalytic_domain/histidine_brace.json`](#histidine_bracejson) | 6f | Brace histidines, metal, pseudo-copper |
| [`07_core_fold/beta_core.json`](#beta_corejson) | 7a | Core strands |
| [`07_core_fold/loop_regions.json`](#loop_regionsjson) | 7b, 7c | Loops and binding surface as ranges |
| [`07_core_fold/binding_surface.json`](#binding_surfacejson) | 7c | Binding plane and its residues |
| [`08_second_sphere/brace_interactions.json`](#brace_interactionsjson) | 8 | Active-site surroundings |
| [`09_distance_matrix/distance_matrix.csv`](#distance_matrixcsv) | 9 | Pairwise distances |
| [`10_flexibility/flexibility_raw.json`](#flexibility_rawjson) | 10 | Flexibility measures |
| [`11_electrostatics/<id>.pqr`](#pqr-files) | 11 | Charges and radii |
| [`12_hole_hopping/hole_hopping.json`](#hole_hoppingjson) | 12 | Hole-hopping paths |

JSON files keyed by structure use the [`struct_id`](../conventions.md#structure-identifiers).

## Step 1

### cazy_lpmos.csv

--8<-- "_formats/cazy_lpmos.md"

### other_lpmos.csv

--8<-- "_formats/other_lpmos.md"

## Step 2

### pdb_selection.csv

--8<-- "_formats/pdb_selection.md"

### uniprot_ids.json

--8<-- "_formats/uniprot_ids.md"

### uniprot_references.json

--8<-- "_formats/uniprot_references.md"

### domain_annotations.json

--8<-- "_formats/domain_annotations.md"

### taxonomy.json

--8<-- "_formats/taxonomy.md"

### activities.json

--8<-- "_formats/activities.md"

## Step 3

### modelled_residues.json

--8<-- "_formats/modelled_residues.md"

### signal_peptides.json

--8<-- "_formats/signal_peptides.md"

### modified_residues.json

--8<-- "_formats/modified_residues.md"

### Prepared structures

--8<-- "_formats/preprocessed_pdb.md"

## Step 4

### Alignment files

--8<-- "_formats/alignment.md"

## Step 5

### sasa.json

--8<-- "_formats/sasa.md"

### secondary_structure.json

--8<-- "_formats/secondary_structure.md"

### dssp.json

--8<-- "_formats/dssp.md"

## Step 6

### disulfides.json

--8<-- "_formats/disulfides.md"

### catalytic_domain_boundaries.json

--8<-- "_formats/catalytic_domain_boundaries.md"

### histidine_brace.json

--8<-- "_formats/histidine_brace.md"

## Step 7

### beta_core.json

--8<-- "_formats/beta_core.md"

### loop_regions.json

--8<-- "_formats/loop_regions.md"

### binding_surface.json

--8<-- "_formats/binding_surface.md"

## Step 8

### brace_interactions.json

--8<-- "_formats/brace_interactions.md"

## Step 9

### distance_matrix.csv

--8<-- "_formats/distance_matrix.md"

## Step 10

### flexibility_raw.json

--8<-- "_formats/flexibility_raw.md"

## Step 11

### PQR files

--8<-- "_formats/pqr.md"

## Step 12

### hole_hopping.json

--8<-- "_formats/hole_hopping.md"
