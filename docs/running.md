# Running the pipeline

```bash
python lpmocat_pipeline/pipeline.py                         # all 12 steps
python lpmocat_pipeline/pipeline.py --steps scrape,fetch    # by name
python lpmocat_pipeline/pipeline.py --steps 5-12            # by number, ranges allowed
python lpmocat_pipeline/pipeline.py --steps 2,6-9           # mixed
python lpmocat_pipeline/pipeline.py --steps 10 --cabsflex   # run missing CABS-flex simulations
python lpmocat_pipeline/pipeline.py --data-dir /path/to/run # another data folder
```

## Options

| Option | Default | Meaning |
| --- | --- | --- |
| `--steps` | all | Comma-separated step names or numbers; ranges such as `5-12`. Steps always run in pipeline order, whatever order they are given in |
| `--families` | all | Comma-separated families (`AA15,AA16`) for steps 1 and 2. Their output files then hold **only** these families, replacing the full lists |
| `--cabsflex` | off | In step 10, run CABS-flex for structures without a current run |
| `--data-dir` | `data/` | The folder holding `01_proteins` … `12_hole_hopping` and `logs` |

Step names: `scrape`, `fetch`, `prepare`, `align`, `dssp_sasa`, `catalytic_domain`,
`core_fold`, `second_sphere`, `distance_matrix`, `flexibility`, `electrostatics`,
`hole_hopping`.

## Environment variables

| Variable | Used by | Meaning |
| --- | --- | --- |
| `BRENDA_JSON` | steps 1c, 2j | Path of the BRENDA release; default: the newest `brenda_*.json*` in the data folder |
| `NCBI_API_KEY` | step 2i | Sent with NCBI requests when set |
| `MODELLER_PYTHON` | step 3a | Python interpreter of the MODELLER environment |
| `SIGNALP6` | step 3h | The `signalp6` executable |
| `MKDSSP` | step 5b | The `mkdssp` executable |

### The BRENDA file { #brenda }

Steps 1c and 2j read a downloaded BRENDA JSON release (`.json` or `.json.tar.gz`). Only the
four LPMO EC entries are parsed from it. If several files match, the last in alphabetical
order is used. The release identifier is recorded in `activities.json`.

## The data folder

```text
<data-dir>/
  01_proteins/ … 12_hole_hopping/     one folder per step
  logs/pipeline_<date>_<time>.log     one log per run
  brenda_<release>.json.tar.gz        the downloaded BRENDA release (steps 1c, 2j)
  django_outputs/                     written by the Django commands (see Loading into the catalogue)
```

Curation files and bundled tools are *not* in the data folder: they are read from
`lpmocat_pipeline/curation/` and `lpmocat_pipeline/tools/` next to the script, whatever
`--data-dir` is. A second data folder is therefore a full, independent run of the same code
and curation. To rerun later steps on a copy, copy or link the folders of the earlier steps
into it.

## Logs and failures

Everything printed is also written to the log file. Each step prints its start and its
duration. A step that cannot run or fails stops the pipeline with `ERROR in step N`, exit
code 1, and the later steps are not run. Typical causes are a missing input file ("Run the
'…' step first"), a missing required program, or a structure that could not be processed.

`WARNING` lines do not stop the run. They mark things to look at, such as an N-terminus cut
at the first histidine or a brace with distant nitrogens.

## Run times

For 160 structures on a 20-core workstation, a full run takes about an hour without
CABS-flex and about two hours with all simulations. The alignment (25 min) and the database
queries of step 2 (15–20 min) dominate; see the table under [Pipeline steps](steps/index.md).

## What a rerun reuses

Most steps recompute everything and overwrite their files. Three things are reused:

| Step | Reused | Until |
| --- | --- | --- |
| 3a | Rebuilt structures of an earlier run, **only when MODELLER is unavailable** | The gaps to fill differ |
| 10c | CABS-flex runs | The prepared structure's atoms change (fingerprint in `cabsflex_results/<id>/input.sha1`) |
| 12 | Hole-hopping results per structure | The aligned structure, brace, domain end, relative SASA or scripts change (fingerprint in `hole_hopping.json`) |

Several steps also remove files of structures that are no longer selected
(`preprocessed/`, `aligned/`, PQR files, hole-hopping files), so the folders always match
the current selection.

## Standalone scripts

Three parts of steps 1 and 2 can be rerun without the rest, on the existing files:

```bash
python lpmocat_pipeline/steps/other_sources.py [--data-dir DIR] [--brenda FILE]         # 1b–1e
python lpmocat_pipeline/steps/taxonomy.py      [--data-dir DIR]                         # 2i
python lpmocat_pipeline/steps/activities.py    [--data-dir DIR] [--brenda FILE] [--no-csv]   # 2j
```

`other_sources.py` rewrites `other_lpmos.csv` for the existing `cazy_lpmos.csv`, which keeps
the CAZy list fixed while the other lists are updated. `activities.py` also writes
`activities_review.csv`, one claim per row, for curation. `binding_surface.py` can be run on
its own to inspect planes in PyMOL ([step 7c](steps/07-core-fold.md#7c-binding-surface)).

## After the pipeline

The results reach the web app through the [Django import](django-import.md).
