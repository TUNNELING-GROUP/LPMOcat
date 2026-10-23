# Installation

The pipeline runs on Linux. It needs a Python environment, three external programs, and two
optional licensed programs. Some tools come with the repository.

## Python environment

Python 3.10 and the packages pinned in `lpmocat_pipeline/requirements.txt`:

```bash
python3.10 -m venv lpmo-pipeline-env
source lpmo-pipeline-env/bin/activate
pip install -r lpmocat_pipeline/requirements.txt
```

| Package | Version | Used for |
| --- | --- | --- |
| pandas | 2.2.3 | Tables throughout |
| numpy | 2.1.3 | Geometry (steps 6–8, 12) |
| scipy | 1.14.1 | Distances (step 12) |
| requests | 2.32.3 | Database queries (steps 1–3) |
| biopython | 1.84 | Structure parsing, SASA, alignments, BLOSUM62 distances |
| beautifulsoup4 | 4.14.3 | Reading the CAZy pages (step 1a) |
| networkx | 3.4.2 | EHPath (step 12) |
| ProDy | 2.6.1 | GNM and RMSF (step 10) |
| CABS | 3.0.12 (CABSflex_standalone, commit `b847a35`) | CABS-flex simulations (step 10c, only with `--cabsflex`) |

The web app has its own file, `lpmocat_web/requirements.txt` (also used by the Docker
image). The packages both need have the same versions in the two files, so one environment
can hold both:

```bash
pip install -r lpmocat_web/requirements.txt -r lpmocat_pipeline/requirements.txt
```

## External programs

| Program | Version | Step | How the pipeline finds it |
| --- | --- | --- | --- |
| mkdssp (DSSP) | 4.4.11 | 5b, required | `MKDSSP`, else `mkdssp` on `PATH`, else the conda environment `dssp` |
| pdb2pqr | 2.1.1 | 11, required | `pdb2pqr` on `PATH` |
| CABS-flex | see above | 10c | `python3 -m CABS` |

```bash
conda create -n dssp -c conda-forge dssp=4.4.11
sudo apt install pdb2pqr          # Ubuntu 22.04 provides 2.1.1
```

- **mkdssp**: use 4.4.11. The conda-forge build of 4.6.1 failed to start with a library
  error. The pipeline
  supplies the data files mkdssp needs itself (see [step 5b](steps/05-dssp-sasa.md)).
- **pdb2pqr**: the pipeline uses the 2.x command line (`--ff=amber`). Version 3 has other
  option names.

## Optional licensed tools

Both are free for academic use but must be obtained from their authors. The pipeline runs
without them, with the stated consequence.

| Program | Version | Step | How it is found | Without it |
| --- | --- | --- | --- | --- |
| MODELLER | 10.8 | 3a | `MODELLER_PYTHON` (the Python of its environment), else the conda environment `modeller` | Internal gaps of crystal structures are not rebuilt; rebuilt files of an earlier run are reused where they still apply |
| SignalP | 6.0 (fast mode) | 3h | `SIGNALP6`, else `signalp6` on `PATH` | Models without an annotated N-terminus are cut at their first histidine, with a warning |

```bash
conda create -n modeller -c salilab modeller      # then enter the licence key
```

MODELLER runs in its own environment as a subprocess (`modeller_worker.py`), so it need not
be installed into the pipeline's environment. SignalP 6 is installed from the package DTU
Health Tech provides.

## Bundled tools

In `lpmocat_pipeline/tools/`:

| Folder | Content | Licence |
| --- | --- | --- |
| `USalign/` | US-align version 20241108: source (`USalign.cpp`) and a Linux x86-64 binary, unmodified | As stated in the source file |
| `dssp/components.cif` | The 20 standard amino acids from the wwPDB Chemical Component Dictionary, which mkdssp needs to read PDB-format files | wwPDB data |
| `hole_hopping/` | The pipeline's own scripts `01_classify_resids.py`, `02_get_EH_results.py`, `03_analyse_paths.py` (step 12) and `04_create_pymol_sessions.py` (standalone) | |
| `hole_hopping/EHPath/` | EHPath (Teo et al. 2019), with `EHPath.py` modified from the upstream `original_EHPath.py` | GPLv3 |

On another platform, compile US-align from the bundled source:

```bash
g++ -O3 -ffast-math -lm -o USalign USalign.cpp
```

### Changes to EHPath

`EHPath.py` is EHPath (Teo et al. 2019, GPLv3) modified for the LPMO active site. The copper
node follows Sayler et al. 2025 ([doi:10.1073/pnas.2411229121](https://doi.org/10.1073/pnas.2411229121)),
whose code, provided by its authors, is `original_EHPath.py`. `EHPath.py` differs from it as
follows (the full list against the published EHPath is in `MODIFICATIONS.md` next to it):

- input files and the donor residue come from the command line instead of prompts;
- the settings are fixed: `cutoff_num` 5, `total_paths` 5, `type_da` `electron`,
  `alpha_reorg` 1;
- backbone and Cβ atoms of a histidine donor are left out, so the node is its ring;
- node parameters are defined for histidine and copper (λ 1.0 eV, E 1.0 V, radius 1 Å);
- the heme, iron–sulfur and other system-specific parameter blocks are commented out.

## Data to download

**BRENDA** (steps 1c and 2j): download the JSON release from
[brenda-enzymes.org/download.php](https://www.brenda-enzymes.org/download.php) (CC BY 4.0)
and save it as `data/brenda_<release>.json.tar.gz`, or point `BRENDA_JSON` to it.
Without it, the BRENDA source is skipped.

## Checking the installation

```bash
python lpmocat_pipeline/pipeline.py --help
mkdssp --version                 # or: conda run -n dssp mkdssp --version
pdb2pqr --version
signalp6 --version               # optional
python -m CABS flex --help       # optional
lpmocat_pipeline/tools/USalign/USalign | head -3
```

Step 5 names the mkdssp version it uses, and step 3 the MODELLER version, at the start of
their log output.
