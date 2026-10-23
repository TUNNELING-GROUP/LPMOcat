"""
Scraper for CAZy characterized and structure pages for LPMO families (AA9-AA17).

Scrapes both the "characterized" and "structure" tabs for each family,
then merges them into a single deduplicated list with two boolean columns:
  - in_cazy_characterized: protein appears on the characterized page
  - in_cazy_structure: protein appears on the structure page

A protein can appear on both pages (characterized enzyme that also has
a crystal structure) or on just one.

Output: data/01_proteins/cazy_lpmos.csv

Usage:
    python lpmocat_pipeline/steps/scrape_cazy.py
"""

import re
import requests
import pandas as pd
from bs4 import BeautifulSoup
from pathlib import Path

import unicodedata

LPMO_FAMILIES = ["AA9", "AA10", "AA11", "AA13", "AA14", "AA15", "AA16", "AA17"]
CAZY_BASE_URL = "https://www.cazy.org"
from paths import DEFAULT_DATA_DIR

OUTPUT_DIR = Path(DEFAULT_DATA_DIR) / "01_proteins"   # where step 1 keeps it


def _clean_text(text):
    """Clean scraped text: fix broken UTF-8 encoding, replace non-ASCII
    characters with ASCII equivalents, and normalize whitespace."""
    # Fix common mojibake patterns (broken UTF-8 decoded as latin-1)
    replacements = {
        "â€\x93": "-", "â€\x94": "-",  # en-dash, em-dash
        "â€œ": '"', "â€\x9d": '"',     # smart quotes
        "â€˜": "'", "â€™": "'",         # smart apostrophes
        "â€¦": "...",                     # ellipsis
        "Î±": "alpha", "Î²": "beta", "Î³": "gamma", "Î´": "delta",
    }
    for broken, fixed in replacements.items():
        text = text.replace(broken, fixed)
    # Remove any remaining non-printable / control characters
    text = "".join(c for c in text if unicodedata.category(c)[0] != "C"
                   or c in ("\n", "\t"))
    # Replace non-ASCII with closest ASCII or drop
    text = unicodedata.normalize("NFKD", text)
    text = text.encode("ascii", "ignore").decode("ascii")
    # Normalize whitespace
    text = " ".join(text.split())
    return text


def _extract_pdb_codes(cell):
    """Extract unique PDB codes from a table cell containing PDB links."""
    pdb_entries = []
    for a in cell.find_all("a"):
        pdb_text = a.get_text(strip=True)
        match = re.match(r"(\w{4})", pdb_text)
        if match:
            pdb_entries.append(match.group(1))
    return list(dict.fromkeys(pdb_entries))  # unique, preserve order


def _parse_kingdom_row(row):
    """Extract kingdom name from a separator row, or None."""
    if not (row.get("class") and "royaume" in row.get("class", [])):
        return None
    span = row.find("span", class_="tabulation")
    if span:
        text = span.get_text(strip=True)
        if text and text != "Top":
            return text
    return None


# ---------------------------------------------------------------------------
# Characterized page
# ---------------------------------------------------------------------------

def scrape_characterized_page(family):
    """Scrape a single CAZy characterized page.

    Columns: Protein Name | EC# | Reference | Organism | GenBank | Uniprot | PDB/3D
    """
    url = f"{CAZY_BASE_URL}/{family}_characterized.html"
    response = requests.get(url, timeout=30)
    response.raise_for_status()

    soup = BeautifulSoup(response.text, "html.parser")
    tables = soup.find_all("table")
    if len(tables) < 2:
        return []

    proteins = []
    current_kingdom = None

    for row in tables[1].find_all("tr"):
        kingdom = _parse_kingdom_row(row)
        if kingdom is not None:
            current_kingdom = kingdom
            continue

        if row.get("id") == "line_titre":
            continue

        cells = row.find_all("td")
        if len(cells) < 7:
            continue

        ec_links = cells[1].find_all("a")
        ref_links = cells[2].find_all("a")
        genbank_links = cells[4].find_all("a")
        uniprot_links = cells[5].find_all("a")

        proteins.append({
            "family": family,
            "kingdom": current_kingdom,
            "protein_name": _clean_text(cells[0].get_text(strip=True)),
            "ec_numbers": "; ".join(a.get_text(strip=True) for a in ec_links),
            "organism": _clean_text(cells[3].get_text(strip=True)),
            "genbank_ids": "; ".join(a.get_text(strip=True) for a in genbank_links),
            "uniprot_ids": "; ".join(a.get_text(strip=True) for a in uniprot_links),
            "pdb_codes": ", ".join(_extract_pdb_codes(cells[6])),
            "references": "; ".join(a.get("href", "") for a in ref_links),
        })

    return proteins


# ---------------------------------------------------------------------------
# Structure page
# ---------------------------------------------------------------------------

def scrape_structure_page(family):
    """Scrape a single CAZy structure page.

    Columns: Protein Name | EC# | Organism | GenBank | Uniprot | PDB/3D (nested) | (empty)
    Note: no References column; data rows identified by onmouseover attribute.
    """
    url = f"{CAZY_BASE_URL}/{family}_structure.html"
    response = requests.get(url, timeout=30)
    response.raise_for_status()

    soup = BeautifulSoup(response.text, "html.parser")

    proteins = []
    current_kingdom = None

    for row in soup.find_all("tr"):
        kingdom = _parse_kingdom_row(row)
        if kingdom is not None:
            current_kingdom = kingdom
            continue

        if not row.get("onmouseover"):
            continue

        cells = row.find_all("td", recursive=False)
        if len(cells) < 6:
            continue

        ec_links = cells[1].find_all("a")
        genbank_links = cells[3].find_all("a")
        uniprot_links = cells[4].find_all("a")

        proteins.append({
            "family": family,
            "kingdom": current_kingdom,
            "protein_name": _clean_text(cells[0].get_text(strip=True)),
            "ec_numbers": "; ".join(a.get_text(strip=True) for a in ec_links),
            "organism": _clean_text(cells[2].get_text(strip=True)),
            "genbank_ids": "; ".join(a.get_text(strip=True) for a in genbank_links),
            "uniprot_ids": "; ".join(a.get_text(strip=True) for a in uniprot_links),
            "pdb_codes": ", ".join(_extract_pdb_codes(cells[5])),
            "references": "",
        })

    return proteins


# ---------------------------------------------------------------------------
# Merge characterized + structure into one deduplicated list
# ---------------------------------------------------------------------------

def _genbank_set(entry):
    """Return the set of GenBank IDs for an entry."""
    return {g.strip() for g in entry["genbank_ids"].split("; ") if g.strip()}


def _pdb_set(entry):
    """Return the set of uppercase PDB codes for an entry."""
    return {p.strip().upper() for p in entry["pdb_codes"].split(", ") if p.strip()}


def _uniprot_set(entry):
    """Return the set of UniProt IDs for an entry."""
    return {u.strip() for u in entry["uniprot_ids"].split("; ") if u.strip()}


def scrape_all_lpmo_families():
    """Scrape both pages for all families, merge with in_cazy_characterized
    and in_cazy_structure columns."""

    all_characterized = []
    all_structures = []

    for family in LPMO_FAMILIES:
        print(f"Scraping {family} ...")

        try:
            char = scrape_characterized_page(family)
            print(f"  characterized: {len(char)} proteins")
        except Exception as e:
            print(f"  characterized ERROR: {e}")
            char = []

        try:
            struct = scrape_structure_page(family)
            print(f"  structure:     {len(struct)} proteins")
        except Exception as e:
            print(f"  structure ERROR: {e}")
            struct = []

        all_characterized.extend(char)
        all_structures.extend(struct)

    # Build lookup indices for structure entries
    struct_gb_index = {}   # genbank_id -> entry
    struct_pdb_index = {}  # pdb_code -> entry
    struct_uni_index = {}  # uniprot_id -> entry
    for e in all_structures:
        for g in _genbank_set(e):
            struct_gb_index[g] = e
        for p in _pdb_set(e):
            struct_pdb_index[p] = e
        for u in _uniprot_set(e):
            struct_uni_index[u] = e

    # Same for characterized (needed to check structure-only entries)
    char_gb_index = {}
    char_pdb_index = {}
    char_uni_index = {}
    for e in all_characterized:
        for g in _genbank_set(e):
            char_gb_index[g] = e
        for p in _pdb_set(e):
            char_pdb_index[p] = e
        for u in _uniprot_set(e):
            char_uni_index[u] = e

    def _is_in_struct(entry):
        for g in _genbank_set(entry):
            if g in struct_gb_index:
                return True
        for p in _pdb_set(entry):
            if p in struct_pdb_index:
                return True
        for u in _uniprot_set(entry):
            if u in struct_uni_index:
                return True
        return False

    def _is_in_char(entry):
        for g in _genbank_set(entry):
            if g in char_gb_index:
                return True
        for p in _pdb_set(entry):
            if p in char_pdb_index:
                return True
        for u in _uniprot_set(entry):
            if u in char_uni_index:
                return True
        return False

    # Step 1: Add all characterized entries
    merged = []
    seen_genbanks = set()
    seen_pdbs = set()

    for entry in all_characterized:
        gbs = _genbank_set(entry)
        pdbs = _pdb_set(entry)
        if gbs & seen_genbanks or pdbs & seen_pdbs:
            continue
        seen_genbanks.update(gbs)
        seen_pdbs.update(pdbs)

        merged.append({
            **entry,
            "in_cazy_characterized": "yes",
            "in_cazy_structure": "yes" if _is_in_struct(entry) else "no",
        })

    # Step 2: Add structure-only entries (not already added from characterized)
    for entry in all_structures:
        gbs = _genbank_set(entry)
        pdbs = _pdb_set(entry)
        if gbs & seen_genbanks or pdbs & seen_pdbs:
            continue
        seen_genbanks.update(gbs)
        seen_pdbs.update(pdbs)

        merged.append({
            **entry,
            "in_cazy_characterized": "no",
            "in_cazy_structure": "yes",
        })

    return merged


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    proteins = scrape_all_lpmo_families()

    df = pd.DataFrame(proteins)
    df = df.sort_values(["family", "organism", "protein_name"]).reset_index(drop=True)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_path = OUTPUT_DIR / "cazy_lpmos.csv"
    df.to_csv(output_path, index=False)
    print(f"\nSaved {len(proteins)} proteins to {output_path}")

    # Summary
    char_yes = df.in_cazy_characterized == "yes"
    struct_yes = df.in_cazy_structure == "yes"
    with_pdb = df.pdb_codes.notna() & (df.pdb_codes != "")

    print("\n=== SUMMARY ===")
    print(f"Total unique LPMOs:     {len(df)}")
    print(f"  In characterized:     {char_yes.sum()}")
    print(f"  In structure:         {struct_yes.sum()}")
    print(f"  In both:              {(char_yes & struct_yes).sum()}")
    print(f"  Characterized only:   {(char_yes & ~struct_yes).sum()}")
    print(f"  Structure only:       {(~char_yes & struct_yes).sum()}")
    print(f"  With PDB codes:       {with_pdb.sum()}")

    print("\n=== PER FAMILY ===")
    for fam in LPMO_FAMILIES:
        f = df[df.family == fam]
        fc = f.in_cazy_characterized == "yes"
        fs = f.in_cazy_structure == "yes"
        print(
            f"  {fam:5s}: {len(f):3d} total | "
            f"char={fc.sum():2d} struct={fs.sum():2d} both={(fc & fs).sum():2d}"
        )
