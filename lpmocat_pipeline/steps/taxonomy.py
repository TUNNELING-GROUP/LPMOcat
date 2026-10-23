"""
2i: NCBI Taxonomy lineage of the source organisms (E-utilities), for grouping the catalogue
by taxon (fungi, oomycetes, insects, bacterial phyla ...) beyond CAZy's kingdom.

An organism name as CAZy gives it (often with a strain) is looked up as a scientific name,
then as any name (synonyms), then as its species (first two words), then as its genus; the
first query with exactly one hit wins.

Output (02_data/taxonomy.json):
    {organism name: {"taxid", "name", "rank", "matched": "name" | "species" | "genus",
                     "lineage": [[rank, name], ...] from the root down to the taxon itself,
                     "species": {"taxid", "name"} of the taxon's species (the taxon itself
                                when it is one), or null above species level}}
    with null for a name NCBI does not know.

Runs as part of pipeline step 2 (fetch); standalone it fills taxonomy.json for the current
pdb_selection.csv:
    python lpmocat_pipeline/steps/taxonomy.py [--data-dir DIR]

NCBI asks for at most 3 requests per second without an API key (NCBI_API_KEY is used when
set).
"""

import argparse
import json
import os
import time
import xml.etree.ElementTree as ET

import requests

from paths import DEFAULT_DATA_DIR

EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
REQUEST_INTERVAL = 0.34     # s between requests (NCBI: <= 3 per second without a key)
FETCH_BATCH = 100           # taxids per efetch request

_last_request = 0.0


def _get(endpoint, params):
    """GET an E-utilities endpoint, paced and retried (NCBI answers bursts with 429)."""
    global _last_request
    params = {"tool": "lpmocat", **params}
    if os.environ.get("NCBI_API_KEY"):
        params["api_key"] = os.environ["NCBI_API_KEY"]
    for attempt in range(4):
        wait = REQUEST_INTERVAL - (time.monotonic() - _last_request)
        if wait > 0:
            time.sleep(wait)
        _last_request = time.monotonic()
        try:
            r = requests.get(f"{EUTILS_URL}/{endpoint}", params=params, timeout=30)
            if r.status_code == 200:
                return r
        except requests.RequestException:
            pass
        time.sleep(2 ** attempt)
    raise RuntimeError(f"NCBI {endpoint} failed for {params}")


def _search(term):
    """Taxids matching an E-utilities taxonomy query."""
    r = _get("esearch.fcgi", {"db": "taxonomy", "term": term, "retmode": "json"})
    return r.json().get("esearchresult", {}).get("idlist", [])


def _queries(name):
    """(query, how it matched) in the order tried."""
    words = name.replace('"', "").split()
    queries = [(f'"{" ".join(words)}"[Scientific Name]', "name"),
               (f'"{" ".join(words)}"[All Names]', "name")]
    if len(words) > 2:
        species = " ".join(words[:2])
        queries += [(f'"{species}"[Scientific Name]', "species"),
                    (f'"{species}"[All Names]', "species")]
    if len(words) > 1:
        queries.append((f'"{words[0]}"[Scientific Name]', "genus"))
    return queries


def find_taxid(name):
    """(taxid, how it matched) of an organism name, or (None, None)."""
    for term, matched in _queries(name):
        ids = _search(term)
        if len(ids) == 1:
            return int(ids[0]), matched
    return None, None


def fetch_taxa(taxids):
    """{taxid: {"taxid", "name", "rank", "lineage", "species"}} from NCBI Taxonomy records."""
    taxa = {}
    taxids = sorted(set(taxids))
    for i in range(0, len(taxids), FETCH_BATCH):
        batch = taxids[i:i + FETCH_BATCH]
        r = _get("efetch.fcgi", {"db": "taxonomy", "id": ",".join(map(str, batch)),
                                 "retmode": "xml"})
        for taxon in ET.fromstring(r.content).findall("Taxon"):
            taxid = int(taxon.findtext("TaxId"))
            ancestors = taxon.find("LineageEx")
            nodes = (list(ancestors) if ancestors is not None else []) + [taxon]
            lineage = [[node.findtext("Rank"), node.findtext("ScientificName")] for node in nodes]
            species = next(({"taxid": int(node.findtext("TaxId")), "name": node.findtext("ScientificName")}
                            for node in nodes if node.findtext("Rank") == "species"), None)
            taxa[taxid] = {"taxid": taxid, "name": taxon.findtext("ScientificName"),
                           "rank": taxon.findtext("Rank"), "lineage": lineage, "species": species}
    return taxa


def fetch_lineages(organism_names):
    """{organism name: taxonomy entry or None}; see the module docstring."""
    found = {}
    for name in sorted({n.strip() for n in organism_names if isinstance(n, str) and n.strip()}):
        found[name] = find_taxid(name)
    taxa = fetch_taxa(taxid for taxid, _ in found.values() if taxid)
    result = {}
    for name, (taxid, matched) in found.items():
        entry = taxa.get(taxid)
        result[name] = dict(entry, matched=matched) if entry else None
    return result


def rank_name(entry, rank):
    """The lineage's name at `rank` (e.g. 'phylum'), or None."""
    return next((n for r, n in (entry or {}).get("lineage", []) if r == rank), None)


def report(taxonomy):
    """Print how the names were matched and the domain/kingdom/phylum counts."""
    missing = sorted(n for n, e in taxonomy.items() if e is None)
    by_match = {}
    for name, entry in taxonomy.items():
        if entry and entry["matched"] != "name":
            by_match.setdefault(entry["matched"], []).append(f"{name} -> {entry['name']}")
    species = {e["species"]["taxid"] for e in taxonomy.values() if e and e.get("species")}
    print(f"  {len(taxonomy) - len(missing)} of {len(taxonomy)} organisms found in NCBI Taxonomy"
          f" ({len(species)} species)")
    for matched, names in sorted(by_match.items()):
        print(f"  matched by {matched} ({len(names)}): {'; '.join(sorted(names))}")
    if missing:
        print(f"  WARNING: not found ({len(missing)}): {'; '.join(missing)}")
    groups = {}
    for entry in taxonomy.values():
        if entry:
            key = " / ".join(filter(None, (rank_name(entry, "domain") or rank_name(entry, "superkingdom")
                                           or rank_name(entry, "acellular root") or entry["lineage"][0][1],
                                           rank_name(entry, "kingdom"), rank_name(entry, "phylum"))))
            groups[key] = groups.get(key, 0) + 1
    for key, n in sorted(groups.items()):
        print(f"    {n:3d}  {key}")


def main():
    import pandas as pd

    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--data-dir", default=DEFAULT_DATA_DIR)
    args = parser.parse_args()
    structures = os.path.join(args.data_dir, "02_data")
    selection = pd.read_csv(os.path.join(structures, "pdb_selection.csv"))
    names = selection.loc[selection.source != "NONE", "organism"].dropna()
    print(f"2i: NCBI Taxonomy lineage of {names.nunique()} organisms...")
    taxonomy = fetch_lineages(names)
    report(taxonomy)
    out = os.path.join(structures, "taxonomy.json")
    with open(out, "w") as f:
        json.dump(taxonomy, f, indent=1)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
