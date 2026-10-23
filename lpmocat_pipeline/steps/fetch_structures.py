"""
Step 2 of the LPMOcat pipeline: find and download one structure per CAZy entry.

pipeline.step_fetch takes, per entry of cazy_lpmos.csv, the first found of
  2a  the best PDB entry CAZy lists: RCSB metadata ranks X-ray first, then resolution,
      then chain completeness (modeled residues / sample sequence length); downloaded as
      PDB, or mmCIF when there is no PDB-format file
  2b  the AlphaFold DB model of a UniProt ID CAZy lists
  2c  the AlphaFold DB model of a UniProt entry with the same sequence as CAZy's GenBank
      protein, found through UniParc
  2d  the CAZyme3D model of the GenBank protein
then drops identical-sequence duplicates (2e), collects every UniProt accession of the
selected proteins with its status, EC numbers and papers (2f) and fetches their
InterPro/UniProt domain annotations (2g).
"""

import json
import os
import re

import requests

# Pfam IDs for LPMO catalytic domains
LPMO_PFAM_IDS = {
    "PF03443",  # AA9 (formerly GH61)
    "PF03067",  # AA10, AA11, AA13, AA15, AA16 (LPMO_10)
    "PF22810",  # AA14
}


# ---------------------------------------------------------------------------
# 2a: Query RCSB for PDB metadata and select the best entry
# ---------------------------------------------------------------------------

# 2a: the metal ligands that can sit in the histidine brace -- the set preprocessing (3g)
# keeps, plus CU1, the Cu(I) ligand code. Copper is the native metal and ranks first; the
# others are substitutes or soaks.
COPPER_LIGANDS = {"CU", "CU1"}
METAL_LIGANDS = COPPER_LIGANDS | {"ZN", "NI", "CO", "MN", "FE", "MO", "MG"}
# How much resolution a structure may give up to be preferred for its metal (user decision
# 2026-09-23; being wild type counts whatever the resolution since 2026-09-30). It keeps a much worse entry from winning on the metal alone:
# 8C5N 0.75 A stays over the apo 7SQX at 3.0 A.
MAX_RESOLUTION_PENALTY = 1.0

# 2a: entries chosen by hand (curation/structures.csv): one of them wins whenever a protein
# lists it, before the ranking below -- for the cases the rule cannot see, such as the
# high-dose end of an X-ray dose series winning on resolution.
PINNED_STRUCTURES_CSV = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                     "curation", "structures.csv")


def pinned_structures(path=PINNED_STRUCTURES_CSV):
    """{PDB code: reason} of curation/structures.csv ({} without the file)."""
    import csv
    if not os.path.exists(path):
        return {}
    with open(path, newline="") as f:
        return {row["pdb_entry"].strip().upper(): row.get("reason", "")
                for row in csv.DictReader(f) if (row.get("pdb_entry") or "").strip()}


def fetch_pdb_metadata(pdb_ids):
    """Fetch resolution, method, chain info and completeness for a list of PDB IDs.

    Returns dict: pdb_id -> {resolution, method, metals, mutations, chains: [{chain_id,
    seq_length, modeled_count, coverage}]}, where metals are the metal ligands of the entry
    (BRACE_METALS, the set preprocessing keeps) and mutations the engineered mutations its
    polymer entities declare -- both used by select_best_pdb.
    """
    metadata = {}

    # GraphQL for resolution, method, chains
    for i in range(0, len(pdb_ids), 80):
        batch = pdb_ids[i:i + 80]
        query = """
        {
          entries(entry_ids: %s) {
            rcsb_id
            rcsb_entry_info {
              resolution_combined
              experimental_method
            }
            nonpolymer_entities {
              nonpolymer_comp {
                chem_comp {
                  id
                }
              }
            }
            polymer_entities {
              rcsb_polymer_entity {
                pdbx_description
                pdbx_mutation
              }
              rcsb_polymer_entity_container_identifiers {
                asym_ids
                auth_asym_ids
              }
              entity_poly {
                pdbx_strand_id
                rcsb_sample_sequence_length
                rcsb_entity_polymer_type
              }
            }
          }
        }
        """ % json.dumps(batch)

        r = requests.post("https://data.rcsb.org/graphql", json={"query": query},
                          timeout=30)
        r.raise_for_status()
        data = r.json()

        for entry in data["data"]["entries"]:
            pdb_id = entry["rcsb_id"]
            res_list = entry["rcsb_entry_info"]["resolution_combined"]
            resolution = res_list[0] if res_list else None
            method = entry["rcsb_entry_info"]["experimental_method"]

            ligands = {(c.get("nonpolymer_comp") or {}).get("chem_comp", {}).get("id", "").upper()
                       for c in (entry.get("nonpolymer_entities") or [])}
            metals = sorted(ligands & METAL_LIGANDS)
            mutations = [(pe.get("rcsb_polymer_entity") or {}).get("pdbx_mutation")
                         for pe in (entry.get("polymer_entities") or [])]
            mutations = [m for m in mutations if m]

            chains = []
            for pe in (entry.get("polymer_entities") or []):
                ep = pe.get("entity_poly") or {}
                if ep.get("rcsb_entity_polymer_type") != "Protein":
                    continue
                seq_len = ep.get("rcsb_sample_sequence_length", 0)
                desc = (pe.get("rcsb_polymer_entity") or {}).get(
                    "pdbx_description", "")
                # Get auth_asym_ids and label asym_ids mapping
                cids = pe.get(
                    "rcsb_polymer_entity_container_identifiers") or {}
                asym_ids = cids.get("asym_ids") or []
                auth_asym_ids = cids.get("auth_asym_ids") or []
                for idx, auth_id in enumerate(auth_asym_ids):
                    label_id = asym_ids[idx] if idx < len(asym_ids) else auth_id
                    chains.append({
                        "chain_id": auth_id,
                        "asym_id": label_id,
                        "seq_length": seq_len,
                        "description": desc,
                    })

            metadata[pdb_id] = {
                "resolution": resolution,
                "method": method,
                "metals": metals,
                "mutations": mutations,
                "chains": chains,
            }

        print(f"  Fetched metadata batch {i // 80 + 1}"
              f" ({len(batch)} entries)")

    # REST calls for modeled_residue_count per chain
    # Note: REST API uses label asym_id, not auth asym_id
    for pdb_id, meta in metadata.items():
        for chain in meta["chains"]:
            asym_id = chain.get("asym_id", chain["chain_id"])
            try:
                r = requests.get(
                    f"https://data.rcsb.org/rest/v1/core/"
                    f"polymer_entity_instance/{pdb_id}/{asym_id}",
                    timeout=10)
                if r.status_code == 200:
                    d = r.json()
                    modeled = (d.get("rcsb_polymer_instance_info") or {}).get(
                        "modeled_residue_count", 0)
                    chain["modeled_count"] = modeled
                    chain["coverage"] = (
                        modeled / chain["seq_length"]
                        if chain["seq_length"] > 0 else 0
                    )
                else:
                    chain["modeled_count"] = 0
                    chain["coverage"] = 0
            except Exception:
                chain["modeled_count"] = 0
                chain["coverage"] = 0

    return metadata


def select_best_pdb(pdb_codes_str, all_metadata, pinned=None):
    """Select the best PDB code and chain for alignment.

    Criteria (in order):
      1. Wild type over an entry with engineered mutations
      2. Method: X-ray > Neutron > Electron microscopy > NMR > other
      3. The active-site metal: copper > another brace metal > none
      4. Resolution: lower is better
      5. Chain completeness: higher coverage is better (fewer missing residues)
      6. Chain length: longer protein chain preferred

    1 comes first, whatever the resolution (user decision 2026-09-30): an entry with engineered
    mutations is not the protein, and CAZy and UniProt list a protein's mutant structures
    under it. A mutant is taken only when the protein has no wild-type entry at all.
    3 comes before resolution because the catalogue is about the active site: the brace
    geometry, the second-sphere distances and the shielding residue's angle are measured to
    the copper. It only decides among entries within MAX_RESOLUTION_PENALTY of the best
    resolution available, so a poor structure never wins on its metal alone.

    `pinned` ({PDB code: reason}, curation/structures.csv): when the protein lists one of these
    entries (and RCSB has it), only the pinned entries compete, so a curated choice wins.

    Returns (pdb_code, chain_id, metadata_dict) or (None, None, None).
    """
    codes = [c.strip().upper() for c in pdb_codes_str.split(",") if c.strip()]
    if pinned:
        chosen = [c for c in codes if c in pinned and all_metadata.get(c) and all_metadata[c]["chains"]]
        if chosen:
            codes = chosen

    METHOD_RANK = {
        "X-ray": 0, "Neutron": 1, "Electron Microscopy": 2, "NMR": 3
    }

    def resolution_of(meta):
        return meta["resolution"] if meta["resolution"] is not None else 99

    def metal_rank(meta):
        metals = set(meta.get("metals") or ())
        if metals & COPPER_LIGANDS:
            return 0
        return 1 if metals else 2

    usable = [(code, all_metadata[code]) for code in codes
              if all_metadata.get(code) and all_metadata[code]["chains"]]
    if not usable:
        return None, None, None
    # the metal window is measured among the entries that can win: the wild-type ones if any
    wild_type = [meta for _, meta in usable if not meta.get("mutations")]
    best_resolution = min(resolution_of(meta) for meta in (wild_type or [m for _, m in usable]))

    best = None
    best_score = None

    for code, meta in usable:
        mutant = bool(meta.get("mutations"))
        method_score = METHOD_RANK.get(meta["method"], 4)
        resolution = resolution_of(meta)
        # only a structure of comparable quality may win on its metal
        close = resolution <= best_resolution + MAX_RESOLUTION_PENALTY
        metal = metal_rank(meta) if close else 3

        for chain in meta["chains"]:
            coverage = chain.get("coverage", 0)
            seq_len = chain.get("seq_length", 0)
            score = (mutant, method_score, metal, resolution, -coverage, -seq_len)

            if best_score is None or score < best_score:
                best_score = score
                best = (code, chain["chain_id"], meta, chain)

    if best:
        code, chain_id, meta, chain = best
        return code, chain_id, {
            "pdb_code": code,
            "chain_id": chain_id,
            "resolution": meta["resolution"],
            "method": meta["method"],
            "seq_length": chain["seq_length"],
            "modeled_count": chain.get("modeled_count", ""),
            "coverage": chain.get("coverage", ""),
        }
    return None, None, None


# ---------------------------------------------------------------------------
# 2a-2d: Download structures (RCSB, AlphaFold DB, UniParc, CAZyme3D)
# ---------------------------------------------------------------------------

def download_pdb(pdb_code, outdir):
    """Download a PDB structure from RCSB. Tries PDB format first, falls back
    to mmCIF."""
    # Try PDB format
    url = f"https://files.rcsb.org/download/{pdb_code}.pdb"
    r = requests.get(url, timeout=30)
    if r.status_code == 200 and r.text.startswith(("HEADER", "ATOM", "REMARK")):
        filepath = os.path.join(outdir, f"{pdb_code}.pdb")
        with open(filepath, "w") as f:
            f.write(r.text)
        return filepath

    # Fall back to mmCIF
    url = f"https://files.rcsb.org/download/{pdb_code}.cif"
    r = requests.get(url, timeout=30)
    if r.status_code == 200:
        filepath = os.path.join(outdir, f"{pdb_code}.cif")
        with open(filepath, "w") as f:
            f.write(r.text)
        return filepath

    return None


def check_alphafold_available(uniprot_id):
    """Check if an AlphaFold model exists for the given UniProt ID.
    Returns True/False without downloading."""
    try:
        api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
        r = requests.get(api_url, timeout=10)
        if r.status_code == 200:
            data = r.json()
            if data and data[0].get("pdbUrl"):
                return True
    except Exception:
        pass
    return False


def download_alphafold(uniprot_id, outdir):
    """Download an AlphaFold model from the AlphaFold EBI database."""
    try:
        api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
        r = requests.get(api_url, timeout=15)
        if r.status_code != 200:
            return None

        data = r.json()
        if not data:
            return None

        pdb_url = data[0].get("pdbUrl", "")
        if not pdb_url:
            return None

        r2 = requests.get(pdb_url, timeout=30)
        if r2.status_code != 200:
            return None

        filepath = os.path.join(outdir, f"AF-{uniprot_id}.pdb")
        with open(filepath, "w") as f:
            f.write(r2.text)
        return filepath
    except Exception as e:
        print(f"  Warning: AlphaFold download failed for {uniprot_id}: {e}")
        return None


UNIPARC_URL = "https://rest.uniprot.org/uniparc"


_uniparc_records = {}


def _uniparc_uniprot_refs(db_id):
    """UniParc record of a sequence database ID (GenBank protein or UniProt accession, also a
    deleted one): (UniParc ID, [UniProtKB cross-references: dicts with 'id', 'active',
    'database', 'organism']), or (None, []) when UniParc has no record or the lookup fails.
    A cross-reference is active while the UniProt entry has exactly this sequence."""
    if db_id in _uniparc_records:
        return _uniparc_records[db_id]
    try:
        r = requests.get(f"{UNIPARC_URL}/search", timeout=30, params={
            "query": f"dbid:{db_id}", "fields": "upi", "format": "json", "size": 1})
        r.raise_for_status()
        hits = r.json().get("results") or []
        if not hits:
            _uniparc_records[db_id] = (None, [])
            return None, []
        upi = hits[0]["uniParcId"]
        r = requests.get(f"{UNIPARC_URL}/{upi}/databases", timeout=30, params={
            "dbTypes": "UniProtKB/Swiss-Prot,UniProtKB/TrEMBL", "format": "json", "size": 500})
        r.raise_for_status()
    except requests.RequestException as e:
        print(f"  Warning: UniParc lookup failed for {db_id}: {e}")
        return None, []
    _uniparc_records[db_id] = (upi, r.json().get("results") or [])
    return _uniparc_records[db_id]


def uniparc_uniprot_ids(genbank_id, organism=""):
    """UniProt accessions with a sequence identical to a GenBank protein, via UniParc.

    CAZy often lists no UniProt ID although UniProt holds the same sequence (under another
    strain, or in a newer entry after the original was deleted), and AlphaFold DB models are
    keyed by UniProt accession. Order: active entries first -- Swiss-Prot, then the CAZy
    organism (genus + species), then the rest -- followed by inactive ones, for which
    AlphaFold DB can still serve a model (MaLPMO10C: A0A1C6SBM0).

    Returns (UniParc ID, [accessions]), or (None, []) when UniParc has no record.
    """
    upi, refs = _uniparc_uniprot_refs(genbank_id)
    species = " ".join(organism.split()[:2]).lower()

    def rank(ref):
        name = (ref.get("organism") or {}).get("scientificName", "")
        return (not ref.get("active"),
                ref.get("database") != "UniProtKB/Swiss-Prot",
                " ".join(name.split()[:2]).lower() != species)

    accessions = []
    for ref in sorted(refs, key=rank):
        if ref["id"] not in accessions:
            accessions.append(ref["id"])
    return upi, accessions


def uniparc_active_accessions(accession):
    """Active UniProt accessions with the same sequence as `accession`, Swiss-Prot first.

    UniProt deletes entries (18 of the catalogue's accessions by 2026-09) while AlphaFold DB
    keeps serving their models. A deleted entry has no features and InterPro no longer lists
    it, but an active entry with the identical sequence has the same annotations in the same
    residue numbering."""
    _, refs = _uniparc_uniprot_refs(accession)
    active = sorted((ref for ref in refs if ref.get("active") and ref["id"] != accession),
                    key=lambda ref: ref.get("database") != "UniProtKB/Swiss-Prot")
    return list(dict.fromkeys(ref["id"] for ref in active))


# ---------------------------------------------------------------------------
# 2f: UniProt accessions of a catalogue protein
# ---------------------------------------------------------------------------

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb"

_uniprot_session = requests.Session()
_uniprot_session.mount("https://", requests.adapters.HTTPAdapter(max_retries=requests.adapters.Retry(
    total=3, backoff_factor=2, status_forcelist=(429, 500, 502, 503, 504),
    allowed_methods=("GET",))))
_uniprot_entries = {}

# Reference scopes UniProt curators give to papers about the protein itself. Genome
# sequencing papers ("NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA]") have none of them.
CURATED_REFERENCE_SCOPES = ("FUNCTION", "CATALYTIC ACTIVITY", "X-RAY CRYSTALLOGRAPHY", "STRUCTURE",
                            "BIOPHYSICOCHEMICAL", "COFACTOR", "METHYLATION", "MUTAGENESIS",
                            "SUBSTRATE")


def uniprot_entry(accession):
    """The UniProtKB JSON entry of an accession, also of a deleted one (cached), or None when
    the lookup fails."""
    if accession not in _uniprot_entries:
        try:
            r = _uniprot_session.get(f"{UNIPROT_URL}/{accession}.json", timeout=30)
            r.raise_for_status()
            _uniprot_entries[accession] = r.json()
        except (requests.RequestException, ValueError) as e:
            print(f"  Warning: UniProt lookup failed for {accession}: {e}")
            return None
    return _uniprot_entries[accession]


def _first_author_surname(name):
    """'van den Berg B.' -> 'van den Berg' (UniProt writes initials last)."""
    words = [w for w in name.split() if not w.endswith(".")]
    return " ".join(words) or name


def uniprot_summary(accession):
    """Status, organism, EC numbers and papers of a UniProt entry:
    {"status": "reviewed" | "unreviewed" | "deleted", "organism", "ec_numbers": [...],
     "references": [{"pmid", "doi", "author", "n_authors", "year", "title", "journal",
                     "scopes"}]}, or None when the lookup fails.

    Papers are those UniProt cites for a curated scope (CURATED_REFERENCE_SCOPES), with a
    PubMed ID or DOI, and not marked as describing an incorrect function."""
    entry = uniprot_entry(accession)
    if entry is None:
        return None
    entry_type = entry.get("entryType", "")
    if entry_type == "Inactive":
        return {"status": "deleted", "organism": None, "ec_numbers": [], "references": [],
                "deleted_reason": (entry.get("inactiveReason") or {}).get("deletedReason")}

    description = entry.get("proteinDescription", {})
    names = ([description.get("recommendedName") or {}] + description.get("submissionNames", [])
             + description.get("alternativeNames", []))
    ec_numbers = sorted({ec["value"] for name in names for ec in name.get("ecNumbers", [])})

    references = []
    for ref in entry.get("references", []):
        scopes = ref.get("referencePositions", [])
        text = " ".join(scopes).upper()
        if "INCORRECT" in text or not any(scope in text for scope in CURATED_REFERENCE_SCOPES):
            continue
        citation = ref.get("citation", {})
        xrefs = {x["database"]: x["id"] for x in citation.get("citationCrossReferences", [])}
        if not xrefs.get("PubMed") and not xrefs.get("DOI"):
            continue
        authors = citation.get("authors") or []
        references.append({
            "pmid": xrefs.get("PubMed"),
            "doi": xrefs.get("DOI"),
            "author": _first_author_surname(authors[0]) if authors else None,
            "n_authors": len(authors),
            "year": (citation.get("publicationDate") or "")[:4] or None,
            "title": citation.get("title"),
            "journal": citation.get("journal"),
            "scopes": scopes,
        })
    return {"status": "reviewed" if "Swiss-Prot" in entry_type else "unreviewed",
            "organism": (entry.get("organism") or {}).get("scientificName"),
            "ec_numbers": ec_numbers, "references": references}


def uniprot_identity(listed, genbank_ids, variants=()):
    """2f: every UniProt accession of a catalogue protein, with its status.

    listed: [(accession, source)] from CAZy, the AlphaFold DB or CAZyme3D model and UniParc
    (2c), in order of preference; genbank_ids: the GenBank proteins of the CAZy entry. Every
    UniProt entry that currently has the sequence of a listed accession or GenBank protein
    (UniParc) is added as "identical sequence" -- deleted entries of that kind are not, as
    they come from other genomes and carry no data. variants: [(accession, source)] of the same
    protein from another strain (a strain variant of another list, step 2), added as they are.

    Returns {"uniparc": [UniParc IDs], "accessions": {accession: {"source", "status",
    "organism", "ec_numbers", "references"}}}, listed accessions first."""
    upis, same_sequence = [], []
    for db_id in list(genbank_ids) + [acc for acc, _ in listed]:
        upi, refs = _uniparc_uniprot_refs(db_id)
        if upi and upi not in upis:
            upis.append(upi)
            same_sequence += [ref["id"] for ref in refs if ref.get("active")]

    accessions = {}
    for acc, source in listed:
        accessions.setdefault(acc, {"source": source})
    for acc in same_sequence:
        accessions.setdefault(acc, {"source": "identical sequence"})
    for acc, source in variants:
        accessions.setdefault(acc, {"source": source})
    for acc, info in accessions.items():
        info.update(uniprot_summary(acc)
                    or {"status": "unknown", "organism": None, "ec_numbers": [], "references": []})
    accessions = {acc: info for acc, info in accessions.items()
                  if info["source"] != "identical sequence" or info["status"] != "deleted"}
    return {"uniparc": upis, "accessions": accessions}


# CAZyme3D (Yin lab, UNL; doi:10.1016/j.jmb.2025.169001) holds an AlphaFold2 model for CAZy
# entries, keyed by CAZy's own GenBank accession: a copy of the AlphaFold DB model when the
# accession maps to UniProt, otherwise its own prediction. There is no API; the model link
# is read from the entry page.
CAZYME3D_URL = "https://pro.unl.edu/CAZyme3D"

# The server now and then answers a burst of requests with 403; retry with backoff
_cazyme3d_session = requests.Session()
_cazyme3d_session.mount("https://", requests.adapters.HTTPAdapter(max_retries=requests.adapters.Retry(
    total=3, backoff_factor=5, status_forcelist=(403, 429, 500, 502, 503, 504),
    allowed_methods=("GET",))))


def cazyme3d_struct_id(genbank_id):
    """Structure ID of a CAZyme3D model: C3D-<GenBank accession>, with the version dot as an
    underscore -- a dot would break the '<id>.pdb:A' headers of the alignment."""
    return "C3D-" + genbank_id.strip().replace(".", "_")


def download_cazyme3d(genbank_id, outdir):
    """Download the CAZyme3D model for a GenBank accession.

    Last fallback, after AlphaFold DB via CAZy's and UniParc's UniProt IDs: for entries with
    no UniProt entry at all, e.g. ChiS, whose only model is CAZyme3D's own AlphaFold2
    prediction.

    Returns (filepath, info) with info = {"uniprot": mapped UniProt ID or "",
    "plddt": average pLDDT or None}, or (None, None) when CAZyme3D has no model.
    """
    try:
        page = _cazyme3d_session.get(f"{CAZYME3D_URL}/get_entry.php",
                                     params={"cazyid": genbank_id}, timeout=30)
        page.raise_for_status()
        link = re.search(r"mapped_uniprot_cazy/([^'\"/]+\.pdb)", page.text)
        if not link:  # "No records found for CAZyID"
            return None, None
        r = _cazyme3d_session.get(f"{CAZYME3D_URL}/mapped_uniprot_cazy/{link.group(1)}", timeout=60)
        r.raise_for_status()
    except requests.RequestException as e:
        print(f"  Warning: CAZyme3D download failed for {genbank_id}: {e}")
        return None, None
    if not re.search(r"^ATOM  ", r.text, flags=re.M):
        return None, None

    filepath = os.path.join(outdir, f"{cazyme3d_struct_id(genbank_id)}.pdb")
    with open(filepath, "w") as f:
        f.write(r.text)
    uniprot = re.match(r"AF-([A-Z0-9]+)-F\d+-model", link.group(1))
    plddt = re.search(r"pLDDT = ([\d.]+)", page.text)
    return filepath, {"uniprot": uniprot.group(1) if uniprot else "",
                      "plddt": float(plddt.group(1)) if plddt else None}


# ---------------------------------------------------------------------------
# 2g: Fetch domain annotations from InterPro + UniProt
# ---------------------------------------------------------------------------

def fetch_domain_annotations(uniprot_id):
    """Fetch domain annotations for a UniProt ID from InterPro and UniProt.

    Returns dict with:
      - domains: list of {source, accession, name, start, end}
      - signal_peptide: {start, end} or None
      - chain: {start, end, description} or None
      - catalytic_domain: {start, end, pfam_id, name} or None (LPMO domain)
    """
    result = {
        "uniprot_id": uniprot_id,
        "domains": [],
        "signal_peptide": None,
        "chain": None,
        "catalytic_domain": None,
    }

    # 1. InterPro API — all domain annotations with boundaries
    try:
        r = requests.get(
            f"https://www.ebi.ac.uk/interpro/api/entry/all/protein/"
            f"uniprot/{uniprot_id}",
            timeout=15)
        if r.status_code == 200:
            data = r.json()
            for entry in data.get("results", []):
                meta = entry.get("metadata", {})
                source = meta.get("source_database", "")
                accession = meta.get("accession", "")
                name = meta.get("name", "")
                for prot in entry.get("proteins", []):
                    for loc in prot.get("entry_protein_locations", []):
                        for frag in loc.get("fragments", []):
                            domain = {
                                "source": source,
                                "accession": accession,
                                "name": name,
                                "start": frag["start"],
                                "end": frag["end"],
                            }
                            result["domains"].append(domain)
                            # Check if this is the LPMO catalytic domain
                            if (source == "pfam"
                                    and accession in LPMO_PFAM_IDS):
                                result["catalytic_domain"] = {
                                    "start": frag["start"],
                                    "end": frag["end"],
                                    "pfam_id": accession,
                                    "name": name,
                                }
        elif r.status_code != 204:   # 204: InterPro does not list the protein
            print(f"  Warning: InterPro returned HTTP {r.status_code} for {uniprot_id}")
    except Exception as e:
        print(f"  Warning: InterPro lookup failed for {uniprot_id}: {e}")

    # 2. UniProt API — signal peptide, chain boundaries
    try:
        data = uniprot_entry(uniprot_id)
        if data is not None:
            for feat in data.get("features", []):
                ftype = feat["type"]
                loc = feat.get("location", {})
                start = loc.get("start", {}).get("value")
                end = loc.get("end", {}).get("value")
                desc = feat.get("description", "")

                if ftype == "Signal" and start and end:
                    result["signal_peptide"] = {
                        "start": int(start), "end": int(end)}
                elif ftype == "Chain" and start and end:
                    result["chain"] = {
                        "start": int(start), "end": int(end),
                        "description": desc}
                elif ftype == "Domain" and start and end:
                    # Add UniProt domain annotations (may overlap with
                    # InterPro but provides extra info like domain names)
                    result["domains"].append({
                        "source": "uniprot",
                        "accession": "",
                        "name": desc,
                        "start": int(start),
                        "end": int(end),
                    })
    except Exception as e:
        print(f"  Warning: UniProt lookup failed for {uniprot_id}: {e}")

    return result


def _has_annotation(annotation):
    return bool(annotation["catalytic_domain"] or annotation["chain"] or annotation["domains"])


def fetch_all_domain_annotations(uniprot_ids):
    """Fetch domain annotations for multiple UniProt IDs.

    An accession without any annotation (a deleted UniProt entry) gets those of the first
    active entry with the same sequence that has some (uniparc_active_accessions), recorded
    as "annotation_from".

    Returns dict: uniprot_id -> annotation dict
    """
    annotations = {}
    for uid in uniprot_ids:
        if not uid:
            continue
        annotation = fetch_domain_annotations(uid)
        if not _has_annotation(annotation):
            for other in uniparc_active_accessions(uid)[:10]:
                other_annotation = fetch_domain_annotations(other)
                if _has_annotation(other_annotation):
                    annotation = {**other_annotation, "uniprot_id": uid, "annotation_from": other}
                    break
        annotations[uid] = annotation
    return annotations


# Local models (step 1e, curation/local_models.csv): characterized LPMOs no public database holds
# (JGI-only), modelled for the catalogue. AlphaFold Server writes the Cu as a chain of its own; the
# preprocessing keeps the selected chain, so the metal moves into the protein chain.
LOCAL_MODEL_METALS = {"CU", "ZN", "NI", "CO", "MN", "FE"}


def import_local_model(path, struct_id, outdir):
    """Write a local model (mmCIF or PDB) as <outdir>/<struct_id>.pdb: the first model, the protein
    as chain A with its B-factors (pLDDT) and the metal ions of other chains moved into chain A.
    Returns the written path."""
    from Bio.PDB import MMCIFParser, PDBParser, PDBIO

    parser = MMCIFParser(QUIET=True) if path.lower().endswith(".cif") else PDBParser(QUIET=True)
    model = parser.get_structure(struct_id, path)[0]
    chains = list(model)
    protein = chains[0]
    for chain in chains[1:]:
        for residue in list(chain):
            if residue.get_resname().strip().upper() in LOCAL_MODEL_METALS:
                chain.detach_child(residue.id)
                last = max(r.id[1] for r in protein)
                residue.id = (f"H_{residue.get_resname().strip()}", last + 1, " ")
                protein.add(residue)
        if not len(chain):
            model.detach_child(chain.id)
    if protein.id != "A":
        protein.id = "A"
    io = PDBIO()
    io.set_structure(model)
    out = os.path.join(outdir, f"{struct_id}.pdb")
    io.save(out)
    return out
