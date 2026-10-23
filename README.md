# LPMO Catalogue (LPMOcat)

A web catalogue of lytic polysaccharide monooxygenases (LPMOs) for browsing the families and
comparing their structures, active sites and properties side by side:
**<https://lpmo-cat.polsl.pl>**. Its About page describes the content, how to use it and how the
catalogue is built.

This repository holds the data pipeline (`lpmocat_pipeline/`) that collects the proteins,
prepares and aligns their structures and computes their properties, and the Django web app
(`lpmocat_web/`) that loads and displays the results.

## How the data are made

`lpmocat_pipeline/pipeline.py` runs 12 steps; each writes into its own folder under `data/` and
reads only the output of earlier steps:

| Step | What it does |
|------|--------------|
| 1 scrape | CAZy's characterized and structure lists; characterized LPMOs CAZy misses (Swiss-Prot, BRENDA, curated lists, own models) |
| 2 fetch | best PDB entry → AlphaFold DB → CAZyme3D; UniProt, InterPro, EC numbers, NCBI Taxonomy, substrates and oxidised positions with their papers |
| 3 prepare | gap filling (MODELLER), protein + active-site metal only, N-terminal trim to the brace His, C-terminal tag removal |
| 4 align | structural alignment of all structures (US-align) |
| 5 dssp_sasa | DSSP secondary structure, per-residue SASA and relative SASA |
| 6 catalytic_domain | disulfide bonds, catalytic domain end, histidine brace |
| 7 core_fold | β-sandwich core strands, loop regions, binding surface |
| 8 second_sphere | second-sphere residues and motifs, shielding residue, brace histidine motifs |
| 9 distance_matrix | BLOSUM62 distances between catalytic domains (the trees) |
| 10 flexibility | B-factors / pLDDT, GNM, CABS-flex |
| 11 electrostatics | PDB2PQR (AMBER) |
| 12 hole_hopping | hole-hopping paths to the copper (EHPath) |

Every step, its rules, parameters, curation files and output formats are described in the
**pipeline documentation** (`docs/`), served by the app at
[/docs/](https://lpmo-cat.polsl.pl/docs/).

## Running the app

The repository holds the code. The data the app shows (the database and a few pipeline outputs)
are not in git: either use the **data bundle**, a 44 MB archive of the current release, or make
them yourself with the pipeline.

### 1. Locally, with the data bundle (no calculations)

```bash
git clone https://github.com/TUNNELING-GROUP/LPMOcat.git && cd LPMOcat
python3 -m venv venv && source venv/bin/activate
pip install -r lpmocat_web/requirements.txt
tar -xzf lpmocat-data-<release>.tar.gz     # -> lpmocat_web/db.sqlite3 and data/
cd lpmocat_web
python manage.py migrate
python manage.py runserver                 # http://127.0.0.1:8000/
```

The bundle (`deploy/make_data_bundle.sh` makes it) holds `lpmocat_web/db.sqlite3`,
`data/06_catalytic_domain/catalytic_domain_boundaries.json`,
`data/09_distance_matrix/distance_matrix.csv` and `data/django_outputs/`. Without a `.env`
file the app runs in development mode.

### 2. Locally, with your own data (running the pipeline)

1. Install the pipeline's packages and programs (Python 3.10, mkdssp, PDB2PQR, CABS-flex;
   optionally MODELLER and SignalP 6; a BRENDA release): see *Installation* in the pipeline
   documentation (`docs/installation.md`).
2. Run the pipeline, which takes hours, mostly in the database queries and CABS-flex:
   ```bash
   pip install -r lpmocat_web/requirements.txt -r lpmocat_pipeline/requirements.txt
   python lpmocat_pipeline/pipeline.py --cabsflex
   ```
3. Load the results into the database (*Loading into the catalogue*, `docs/django-import.md`):
   ```bash
   cd lpmocat_web
   python manage.py migrate
   python manage.py populate_db --force
   python manage.py update_protein_codes
   python manage.py rename_lpmo_to_aa --apply
   python manage.py clean_other_names --apply
   python manage.py backfill_references
   python manage.py backfill_cazy_references
   python manage.py map_domain_ranges
   python manage.py pregenerate_plots --force
   python manage.py build_downloads          # the Downloads page and the data release
   python manage.py runserver
   ```

The databases change over time, so a new run gives a newer catalogue than the bundle's.

### 3. On a server (Docker)

The site runs in Docker behind the server's web server (reverse proxy): the code comes by git,
the data by the bundle, mounted into the container, and the settings from `.env`.

First installation:

```bash
git clone --depth 1 https://github.com/TUNNELING-GROUP/LPMOcat.git && cd LPMOcat
cp .env.example .env                       # set the secret key and the host names
tar -xzf lpmocat-data-<release>.tar.gz
docker compose up -d --build               # listens on 127.0.0.1:${LPMOCAT_PORT:-8000}
docker compose exec web python manage.py warm_compare_cache   # whole-family comparisons, optional
docker compose exec web python manage.py createsuperuser      # an admin account, optional
```

New code:

```bash
git pull
docker compose up -d --build
docker builder prune -f                    # the previous build's cache
```

New data (a new bundle, after a pipeline run and the Django import on another machine):

```bash
docker compose down -v                     # -v also empties the plot cache of the old data
cp lpmocat_web/db.sqlite3 db-before-update.sqlite3
rm -r data lpmocat_web/db.sqlite3
tar -xzf lpmocat-data-<release>.tar.gz
docker compose up -d
docker compose exec web python manage.py warm_compare_cache
```

The container refuses to start, and names the missing file, when the bundle has not been
unpacked. The pipeline is not meant to run on the web server: it needs hours of CPU, its own
tools and licences, and a review of its results before they are published.

### Documentation and tests

```bash
pip install -r docs/requirements.txt
mkdocs build --strict                      # -> lpmocat_web/site_root/docs/, served at /docs/
mkdocs serve                               # or a live preview
cd lpmocat_web && python manage.py test catalogue
```

## Project structure

```
lpmo-catalogue/
├── lpmocat_pipeline/         # the data pipeline
│   ├── pipeline.py           # the entry point: runs steps 1-12
│   ├── steps/                # the steps' modules (scrape_cazy, fetch_structures, ...)
│   ├── curation/             # hand-curated inputs (extra proteins, own models, activities, trims)
│   ├── tools/                # bundled tools (US-align, DSSP dictionary, hole hopping)
│   ├── analysis/             # standalone analysis scripts
│   └── requirements.txt
├── lpmocat_web/              # the Django web app
│   ├── manage.py
│   ├── config/               # settings, urls, wsgi
│   ├── catalogue/            # the app: models, views, templates, static files, management
│   │                         # commands (populate_db, ...), tests, data/ (hydrophobicity scales)
│   └── requirements.txt
├── data/                     # what the pipeline writes (not in git): 01_proteins/ …
│                             # 12_hole_hopping/, django_outputs/, logs/
├── docs/                     # the pipeline documentation (MkDocs)
├── deploy/                   # make_data_bundle.sh
├── graphics/                 # source of the LPMO graphic on the About page
└── Dockerfile, docker-compose.yml, .env.example
```

## Data sources

CAZy, RCSB PDB, AlphaFold DB, CAZyme3D, UniProt (Swiss-Prot and TrEMBL), InterPro, NCBI
Taxonomy, SIFTS, Europe PMC, BRENDA (CC BY 4.0) and the IUBMB enzyme list (ExplorEnz). The
About page of the app lists them with references, together with the methods and their papers.

## License

- **Code**: MIT licence ([LICENSE](LICENSE)).
- **Data** (the database, the pipeline outputs, the data bundle and the website's downloads):
  CC BY 4.0 ([LICENSE-DATA.md](LICENSE-DATA.md)); the sources keep their own terms.
- **Bundled third-party code** keeps its licence: EHPath in
  `lpmocat_pipeline/tools/hole_hopping/EHPath/` is GPLv3 (its `LICENSE`; our changes are listed in
  its `MODIFICATIONS.md`), US-align in `lpmocat_pipeline/tools/USalign/` has the licence in the
  header of `USalign.cpp`.
- **Own AlphaFold 3 models** in `lpmocat_pipeline/curation/models/`: AlphaFold Server Output Terms
  of Use, non-commercial use only.

## Contact

lpmocat@tunnelinggroup.pl
