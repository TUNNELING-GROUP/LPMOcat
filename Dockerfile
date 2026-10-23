
#     # dockerfile for building container. this is not production ready

# FROM ubuntu:latest

# RUN mkdir app/
# COPY lpmocat_web /app/lpmocat_web
# WORKDIR app/lpmocat_web
# RUN apt update && apt -y install --no-install-recommends\
#                 python3\
#                 python3-pip
# RUN pip install --no-cache-dir --break-system-packages -r requirements.txt
# RUN python3 manage.py collectstatic --noinput
# EXPOSE 8000
# RUN chmod +x ./entrypoint.sh
# ENTRYPOINT /app/lpmocat_web/entrypoint.sh


# The pipeline documentation, built from docs/ with MkDocs; only the built pages go into the
# image below. --strict fails the build on a broken link.
FROM python:3.12-slim AS docs
WORKDIR /src
COPY docs/requirements.txt docs/requirements.txt
RUN pip install --no-cache-dir -r docs/requirements.txt
COPY mkdocs.yml .
COPY docs docs
RUN mkdocs build --strict --site-dir /site

FROM python:3.12-slim AS base

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc g++ make \
    && rm -rf /var/lib/apt/lists/*

COPY lpmocat_web/requirements.txt requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

COPY lpmocat_web /app/lpmocat_web
# served at /docs/ (WHITENOISE_ROOT in settings.py); pre-compressed like the static files
COPY --from=docs /site /app/lpmocat_web/site_root/docs
RUN python -m whitenoise.compress /app/lpmocat_web/site_root/docs
# so that `docker compose exec web python manage.py ...` runs where manage.py is
WORKDIR /app/lpmocat_web

COPY lpmocat_web/entrypoint.sh /entrypoint.sh
RUN sed -i 's/\r$//' /entrypoint.sh && chmod +x /entrypoint.sh

EXPOSE 8000

ENTRYPOINT ["/entrypoint.sh"]

