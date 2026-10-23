#!/bin/sh

set -e

cd /app/lpmocat_web

# The data are mounted from the host (docker-compose.yml); they come with the data bundle
# (README, "Running the app"), not with git. A mount whose source is missing turns up as an empty directory,
# so check for every one before starting.
missing=""
[ -f /app/lpmocat_web/db.sqlite3 ] || missing="$missing lpmocat_web/db.sqlite3"
[ -f /app/data/06_catalytic_domain/catalytic_domain_boundaries.json ] || missing="$missing data/06_catalytic_domain/catalytic_domain_boundaries.json"
[ -f /app/data/09_distance_matrix/distance_matrix.csv ] || missing="$missing data/09_distance_matrix/distance_matrix.csv"
[ -d /app/data/django_outputs/pregenerated_plots ] || missing="$missing data/django_outputs/pregenerated_plots/"
if [ -n "$missing" ]; then
    echo "LPMOcat: data missing on the host (unpack the data bundle, see the README):" >&2
    for f in $missing; do echo "  $f" >&2; done
    exit 1
fi

python manage.py migrate --noinput
python manage.py collectstatic --noinput --clear

exec gunicorn config.wsgi:application --bind 0.0.0.0:8000 --workers 4 --timeout 400
