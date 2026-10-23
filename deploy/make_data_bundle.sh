#!/bin/sh
# Pack the data the deployed app needs but git does not carry: the database and the pipeline
# outputs the app reads at run time. Run from the repository root on the machine that holds the
# data; copy the archive to the server and unpack it in the repository root there (README, "Running the app").
#
#   sh deploy/make_data_bundle.sh            -> lpmocat-data-<date>-<commit>.tar.gz
set -e
cd "$(dirname "$0")/.."

for f in lpmocat_web/db.sqlite3 \
         data/06_catalytic_domain/catalytic_domain_boundaries.json \
         data/09_distance_matrix/distance_matrix.csv \
         data/django_outputs/pregenerated_plots; do
    [ -e "$f" ] || { echo "missing: $f" >&2; exit 1; }
done

stage=$(mktemp -d)
trap 'rm -rf "$stage"' EXIT
commit=$(git rev-parse --short HEAD)
name="lpmocat-data-$(date +%Y%m%d)-$commit.tar.gz"

# A consistent copy of the database even while a local server has it open
mkdir -p "$stage/lpmocat_web"
python3 -c "import sqlite3, sys; s = sqlite3.connect(sys.argv[1]); d = sqlite3.connect(sys.argv[2]); s.backup(d); d.close()" \
    lpmocat_web/db.sqlite3 "$stage/lpmocat_web/db.sqlite3"

mkdir -p "$stage/data/06_catalytic_domain" "$stage/data/09_distance_matrix"
cp data/06_catalytic_domain/catalytic_domain_boundaries.json "$stage/data/06_catalytic_domain/"
cp data/09_distance_matrix/distance_matrix.csv "$stage/data/09_distance_matrix/"
cp -r data/django_outputs "$stage/data/"

printf 'LPMOcat data bundle\nmade:   %s\ncommit: %s (%s)\nproteins: %s\n' \
    "$(date -Iseconds)" "$commit" "$(git rev-parse --abbrev-ref HEAD)" \
    "$(python3 -c "import sqlite3; print(sqlite3.connect('$stage/lpmocat_web/db.sqlite3').execute('select count(*) from catalogue_protein').fetchone()[0])")" \
    > "$stage/DATA_BUNDLE.txt"

tar -czf "$name" -C "$stage" lpmocat_web data DATA_BUNDLE.txt
echo "wrote $name ($(du -h "$name" | cut -f1))"
cat "$stage/DATA_BUNDLE.txt"
