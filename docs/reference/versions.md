# Versions and data releases

## Software

| Software | Version | Used in |
| --- | --- | --- |
| Python | 3.10 | Pipeline |
| Biopython | 1.84 | Steps 1–3, 5, 9, 12 |
| pandas / numpy / scipy | 2.2.3 / 2.1.3 / 1.14.1 | Throughout |
| ProDy | 2.6.1 | Step 10 |
| CABS-flex (CABSflex_standalone) | 3.0.12, commit `b847a35` | Step 10c |
| MODELLER | 10.8 | Step 3a |
| SignalP | 6.0h, fast mode | Step 3h |
| US-align | 20241108 | Step 4 |
| mkdssp | 4.4.11 | Step 5b |
| pdb2pqr | 2.1.1 | Step 11 |
| EHPath | Teo et al. 2019, modified (see [Installation](../installation.md#changes-to-ehpath)) | Step 12 |
| AlphaFold 3 | | Models of the in-house database entries (step 1e) |

The exact Python package versions are pinned in `lpmocat_pipeline/requirements.txt`. The
pipeline gives identical output files with numpy 1.26.4 and 2.1.3.

## Databases

Most sources are queried live, so a run reflects the databases on the day it was made. The
pipeline records what it can:

| Source | What is recorded | Where |
| --- | --- | --- |
| UniProt, UniParc | Date of the lookup | `checked` in `uniprot_ids.json` |
| BRENDA | Release and file name | `sources.brenda` in `activities.json` |
| Activities (all sources) | Date | `sources.checked` in `activities.json` |
| MODELLER | Version and settings per rebuilt structure | `method` in `modelled_residues.json` |
| AlphaFold DB | Model version | `TITLE` and `REMARK` records of the downloaded file |
| Every run | Date, time, steps, and tool versions where a step prints them | `logs/pipeline_<date>_<time>.log` |

For the current catalogue data:

| Source | Release or retrieval date |
| --- | --- |
| CAZy (characterized and structure pages) | Read 2026-10-21 |
| UniProtKB / UniParc | Release 2026_03 (2 September 2026), queried 2026-10-21 |
| RCSB PDB | Queried 2026-10-21 |
| AlphaFold DB | Model version 6 (`model_v6` files), queried 2026-10-21 |
| CAZyme3D | Queried 2026-10-21 |
| InterPro | Release 110.0 (3 September 2026), queried 2026-10-21 |
| NCBI Taxonomy | Queried 2026-10-21 |
| BRENDA | Release 2026.1 (`brenda_2026_1.json.tar.gz`) |
| Europe PMC | Queried 2026-10-21 (steps 1c and 2j) |
| wwPDB Chemical Component Dictionary | The 20 standard amino acids taken on 2026-09-16 (entries last modified 2024-09-27) |
| ExplorEnz (EC definitions) | Fixed table in the Django import, checked 2026-09-15 |

!!! note "Release numbers"
    The queries the pipeline makes do not return the release numbers of UniProt, InterPro and
    AlphaFold DB. They were read from the services afterwards (UniProt's `X-UniProt-Release`
    header, the InterPro API's database list, the AlphaFold DB file versions); after a new run,
    note them on the day of the run.

## Rebuilding the catalogue from scratch

1. Install the software ([Installation](../installation.md)) and download the BRENDA release.
2. Run the pipeline with the simulations:

    ```bash
    python lpmocat_pipeline/pipeline.py --cabsflex
    ```

    This takes about two hours for 160 structures. Read the log for `WARNING` lines, in
    particular N-termini cut at the first histidine (step 3h), braces with distant nitrogens
    (6f) and structures without a binding surface (7c).

3. If step 7 stops because the number of core-strand blocks changed, the loop definitions
   need revising ([step 7](../steps/07-core-fold.md#known-limits)).
4. Back up the database, then run the [Django import](../django-import.md) commands in
   order.
5. Note the database release numbers and dates in the table above.

Because the sources are queried live, a rebuild at a later date will not reproduce the
catalogue exactly: CAZy adds proteins, UniProt deletes and merges entries, new PDB entries
may be selected, and each added structure changes the alignment columns of all. The rules
are reproducible; the data are a snapshot. To reproduce a given release exactly, keep its
data folder: with `01_proteins/` and `02_data/` in place, `--steps 3-12` needs no database
except the RCSB Chemical Component Dictionary lookup in step 3f.
