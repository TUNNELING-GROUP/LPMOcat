"""
pregenerate_plots.py — Pre-generate per-protein compare-page plots.

Each protein's plots are rendered ONCE against the full MSA coordinate system
(makePlots(..., no_trim=True)). At request time, catalogue.plot_subset.subset_figure
drops the selection's all-gap columns and renumbers — reproducing the per-combination
output exactly (verified), but without rebuilding any Plotly figure. This removes the
multi-minute per-combination build on the compare page.

Only the 4 combination-INDEPENDENT (subsettable) plot types are pre-generated:
Sequence, Charge, Hole Hopping, and Hydrophobicity (one per scale). The 3 line-gradient
types (Flexibility, Secondary Structure and Loops, Relative Surface Accessibility)
draw per-segment traces between column-adjacent residues, so their structure depends on
the selection — those are regenerated on demand instead.

Output: data/django_outputs/pregenerated_plots/<protein_code>/<type>.json.gz
Each file is a gzip-compressed Plotly figure JSON string.

Usage:
  python manage.py pregenerate_plots                 # all proteins, skip existing
  python manage.py pregenerate_plots --force         # regenerate all
  python manage.py pregenerate_plots --protein TrAA9 --protein CsAA9A
"""

import gzip
from pathlib import Path

import plotly.io as pio
from django.core.management.base import BaseCommand

from catalogue.models import Protein
from catalogue.makePlots import makePlots

PROJECT_ROOT = Path(__file__).parents[4]
OUTPUT_DIR = PROJECT_ROOT / "data" / "django_outputs" / "pregenerated_plots"

# Combination-independent plot types safe to pre-generate + subset.
MARKER_TYPES = ["Sequence", "Charge", "Hole Hopping"]


def type_slug(ptype, scale=None):
    base = ptype.lower().replace(" ", "_").replace("/", "_")
    if scale is None:
        return base
    return "hydrophobicity__" + scale.lower().replace(" ", "_").replace("/", "_")


class Command(BaseCommand):
    help = "Pre-generate full-alignment marker plots per protein for fast compare-page assembly."

    def add_arguments(self, parser):
        parser.add_argument("--protein", action="append", dest="proteins",
                            help="Protein code(s) to (re)generate; default is all proteins.")
        parser.add_argument("--force", action="store_true",
                            help="Overwrite existing pre-generated files.")
        parser.add_argument("--type", action="append", dest="types",
                            choices=MARKER_TYPES + ["Hydrophobicity"],
                            help="Plot type(s) to (re)generate, e.g. --type Sequence; default is all.")

    def handle(self, *args, **opts):
        qs = Protein.objects.all()
        if opts.get("proteins"):
            qs = qs.filter(protein_code__in=opts["proteins"])
        qs = qs.exclude(aligned_sequence__isnull=True).exclude(aligned_sequence="")

        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        total = qs.count()
        if total == 0:
            self.stdout.write(self.style.WARNING("No proteins with an aligned_sequence to process."))
            return

        done = 0
        for protein in qs.iterator():
            code = protein.protein_code
            pdir = OUTPUT_DIR / code
            pdir.mkdir(exist_ok=True)

            # Single-protein, full-alignment instance. totalProteinCount is forced >1 so the
            # hover text carries "Alignment position" exactly as in a multi-protein view.
            solo = makePlots([protein], {"name": [code], "seq": [protein.aligned_sequence]},
                             no_trim=True)
            solo.totalProteinCount = 2

            jobs = [(t, None) for t in MARKER_TYPES]
            jobs += [("Hydrophobicity", sc) for sc in solo.hydrophobicityScaleNames]
            if opts.get("types"):
                jobs = [(ptype, scale) for ptype, scale in jobs if ptype in opts["types"]]

            written = 0
            for ptype, scale in jobs:
                out = pdir / (type_slug(ptype, scale) + ".json.gz")
                if out.exists() and not opts["force"]:
                    continue
                fig_json = pio.to_json(
                    solo.make_plot_for_protein(ptype, 0, scale_name=scale),
                    False, False, True, "orjson")
                with gzip.open(out, "wt", encoding="utf-8") as fh:
                    fh.write(fig_json)
                written += 1

            done += 1
            self.stdout.write("[%d/%d] %s (%d files)" % (done, total, code, written))

        self.stdout.write(self.style.SUCCESS(
            "Pre-generated %d protein(s) -> %s" % (done, OUTPUT_DIR)))
