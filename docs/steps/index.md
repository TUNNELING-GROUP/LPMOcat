# Pipeline steps

| | Name for `--steps` | Folder | Code | External tools | Run time¹ |
| --- | --- | --- | --- | --- | --- |
| [1](01-proteins.md) | `scrape` | `01_proteins` | `scrape_cazy.py`, `other_sources.py` | | 20 s |
| [2](02-data.md) | `fetch` | `02_data` | `fetch_structures.py`, `taxonomy.py`, `activities.py` | | 15–20 min |
| [3](03-preparation.md) | `prepare` | `03_preparation` | `fill_gaps.py`, `preprocess_structures.py` | MODELLER², SignalP 6² | 2–3 min |
| [4](04-alignment.md) | `align` | `04_alignment` | `align_structures.py` | US-align | 25 min |
| [5](05-dssp-sasa.md) | `dssp_sasa` | `05_dssp_sasa` | `sasa_dssp.py` | mkdssp | 1–2 min |
| [6](06-catalytic-domain.md) | `catalytic_domain` | `06_catalytic_domain` | `pipeline.py` | | 1 s |
| [7](07-core-fold.md) | `core_fold` | `07_core_fold` | `beta_core.py`, `binding_surface.py` | | 1.5 min |
| [8](08-second-sphere.md) | `second_sphere` | `08_second_sphere` | `brace_interactions.py` | | 2 s |
| [9](09-distance-matrix.md) | `distance_matrix` | `09_distance_matrix` | `pipeline.py` | | 2 min |
| [10](10-flexibility.md) | `flexibility` | `10_flexibility` | `pipeline.py` | ProDy, CABS-flex³ | 1 min; 55 min with all CABS-flex runs |
| [11](11-electrostatics.md) | `electrostatics` | `11_electrostatics` | `pipeline.py` | pdb2pqr | 1 min |
| [12](12-hole-hopping.md) | `hole_hopping` | `12_hole_hopping` | `hole_hopping_ehpath.py` | EHPath (bundled) | 6 min; seconds if unchanged |

¹ For 160 structures on a 20-core workstation. ² Optional; the step continues without them.
³ Only with `--cabsflex`.

Every step function lives in `pipeline.py`; the modules listed hold the methods.

## Which step reads which

| Step | Reads the output of |
| --- | --- |
| 2 | 1 |
| 3 | 2 |
| 4 | 2, 3 |
| 5 | 4 |
| 6 | 2, 4, 5 |
| 7 | 2, 4, 5, 6 |
| 8 | 4, 6 |
| 9 | 4, 6 |
| 10 | 2, 3 |
| 11 | 4, 6 |
| 12 | 4, 5, 6 |

The practical consequences:

- After **step 4** everything from step 5 on must be rerun: the alignment is global and the
  aligned files change.
- After **step 6**, rerun steps 7, 8, 9, 11 and 12.
- **Step 10** depends only on the prepared structures, not on the alignment.

`--steps 5-12` after a new alignment, and `--steps 6-9,11-12` after a change to step 6, are
the usual reruns.
