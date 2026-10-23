"""Where the pipeline's pieces live.

    lpmocat_pipeline/            the code: pipeline.py, steps/ (these modules), curation/, tools/
    data/                        what the pipeline writes (01_proteins ... 12_hole_hopping,
                                 django_outputs/, logs/) and the BRENDA download; git-ignored

pipeline.py and every standalone script take --data-dir to work on another data folder."""
import os

STEPS_DIR = os.path.dirname(os.path.abspath(__file__))
PIPELINE_DIR = os.path.dirname(STEPS_DIR)                  # lpmocat_pipeline/
CURATION_DIR = os.path.join(PIPELINE_DIR, "curation")     # hand-curated inputs, kept with the code
TOOLS_DIR = os.path.join(PIPELINE_DIR, "tools")           # bundled US-align, DSSP dictionary, hole hopping
DEFAULT_DATA_DIR = os.path.join(os.path.dirname(PIPELINE_DIR), "data")
