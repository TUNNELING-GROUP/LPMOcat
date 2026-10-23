"""
1b-1d: characterized LPMOs that CAZy's lists miss (pipeline step 1, after the CAZy scrape).

  1b Swiss-Prot  reviewed UniProt entries with an LPMO EC number (1.14.99.53-56) whose catalytic
                 activity is experimentally shown (ECO:0000269) and cites a paper; the family is
                 the one the entry's curated similarity comment names ("Belongs to the
                 polysaccharide monooxygenase AA9 family"), else that of the most similar
                 CAZy-listed protein
  1c BRENDA      proteins of EC 1.14.99.53-56 in the downloaded BRENDA release that have a
                 UniProt accession and a substrate row citing a paper that names the protein's
                 genus in its title or abstract (Europe PMC; BRENDA attaches the Fusarium
                 fujikuroi entry S0DX13 to a paper on the Bacillus amyloliquefaciens AA10); the
                 family is that of the most similar CAZy-listed protein
  1d curated     lpmocat_pipeline/curation/extra_proteins.csv (accession, protein, family, pmid,
                 reason): unreviewed UniProt entries with published experimental results, as an
                 exception to 1b
  1e local model lpmocat_pipeline/curation/local_models.csv (struct_id, protein, family, organism,
                 kingdom, file, pmid, reason): characterized LPMOs in no public sequence database,
                 with a model file made for the catalogue (BaAA9A and PchAA9E, JGI-only: AlphaFold 3
                 models in curation/models/); step 2 takes the file instead of fetching

Dropped: entries CAZy's lists already hold (a UniProt accession, PDB code or GenBank protein
they list, or the sequence of a UniProt entry or GenBank protein they list - UniParc gives it
for deleted accessions), deleted UniProt entries, entries outside AA9-AA17, and entries with the
sequence of an entry already added or >= 95 % of it in the same species (Swiss-Prot has
T. thermophilus LPMO9B and LPMO9C under two strains; the entry citing more papers comes first,
and step 2f collects the other as an identical-sequence entry). A strain variant (>= 95 % identical to a CAZy-listed
protein of the same species) is kept with that protein's ID in variant_of: step 2 drops it when
the CAZy entry gets a structure, and uses it when not (Swiss-Prot's Q4X066 of A. fumigatus Af293
for AfAA9A, which CAZy lists by a GenBank protein without an AlphaFold DB model). "The most similar CAZy-listed protein" is the best local alignment
(BLOSUM62) to the sequences of the UniProt entries CAZy lists: at least 40 % identity over at
least 100 aligned residues and 5 points above the best alignment to another family. An entry
without such a family is reported, not added (trsA, an Egh16-like protein with an LPMO EC
number, aligns at 36 % to an AA10): add it to curation/extra_proteins.csv with its family if it
is an LPMO. UniProt's own domain annotations cannot tell the families apart (the LPMO_10 domain
PF03067 is in AA10, AA11, AA13, AA15 and AA16, and unreviewed AA15 and AA16 entries are labelled
AA13).

Output (01_proteins/other_lpmos.csv): the columns of cazy_lpmos.csv (in_cazy_characterized and
in_cazy_structure "no"; references as PubMed links, like CAZy's) plus list_source (Swiss-Prot |
BRENDA | curated | local model), evidence (why the entry is in), variant_of, and struct_id and
local_file (a path under lpmocat_pipeline/) for 1e. Step 2 reads it after cazy_lpmos.csv.

Standalone, it writes other_lpmos.csv for the current cazy_lpmos.csv:
    python lpmocat_pipeline/steps/other_sources.py [--data-dir DIR] [--brenda FILE]
"""

import argparse
import os
import re
from collections import defaultdict

import pandas as pd
import requests
from Bio import Align
from Bio.Align import substitution_matrices

from paths import CURATION_DIR, DEFAULT_DATA_DIR, PIPELINE_DIR

CURATED_CSV = os.path.join(CURATION_DIR, "extra_proteins.csv")
LOCAL_MODELS_CSV = os.path.join(CURATION_DIR, "local_models.csv")
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb"
LPMO_EC_NUMBERS = ("1.14.99.53", "1.14.99.54", "1.14.99.55", "1.14.99.56")
FAMILIES = ("AA9", "AA10", "AA11", "AA13", "AA14", "AA15", "AA16", "AA17")
EXPERIMENTAL = "ECO:0000269"
MIN_IDENTITY = 0.40        # family by similarity: identity of the best local alignment ...
MIN_ALIGNED = 100          # ... over at least this many aligned residues ...
MIN_MARGIN = 0.05          # ... and this much above the best alignment to another family
SAME_PROTEIN = 0.95        # identity to a CAZy-listed protein of the same species: a strain variant
COLUMNS = ["family", "kingdom", "protein_name", "ec_numbers", "organism", "genbank_ids",
           "uniprot_ids", "pdb_codes", "references", "in_cazy_characterized", "in_cazy_structure",
           "list_source", "evidence", "variant_of", "struct_id", "local_file"]

_session = requests.Session()
_cazy_species = {}         # species (first two words, lower case) -> as CAZy's organisms write it
_session.mount("https://", requests.adapters.HTTPAdapter(max_retries=requests.adapters.Retry(
    total=5, backoff_factor=1, status_forcelist=(429, 500, 502, 503, 504))))


# ---------------------------------------------------------------------------
# UniProt entries
# ---------------------------------------------------------------------------

def uniprot_entries(accessions, batch=100):
    """{accession: UniProtKB JSON entry} of the active entries among `accessions`."""
    entries = {}
    accessions = sorted(set(accessions))
    for i in range(0, len(accessions), batch):
        r = _session.get(f"{UNIPROT_URL}/accessions", timeout=120,
                         params={"accessions": ",".join(accessions[i:i + batch]), "format": "json"})
        r.raise_for_status()
        for entry in r.json().get("results", []):
            entries[entry["primaryAccession"]] = entry
    return entries


def uniparc_sequence(db_id):
    """The sequence UniParc holds for a UniProt accession (also a deleted one) or a GenBank
    protein, or ''."""
    try:
        r = _session.get("https://rest.uniprot.org/uniparc/search", timeout=60,
                         params={"query": f"dbid:{db_id}", "format": "fasta", "size": 1})
        r.raise_for_status()
    except requests.RequestException:
        return ""
    return "".join(line.strip() for line in r.text.splitlines() if line and not line.startswith(">"))


def swissprot_lpmo_entries():
    """The reviewed entries with an LPMO EC number."""
    query = "reviewed:true AND (" + " OR ".join(f"ec:{ec}" for ec in LPMO_EC_NUMBERS) + ")"
    r = _session.get(f"{UNIPROT_URL}/stream", params={"query": query, "format": "json"}, timeout=300)
    r.raise_for_status()
    return {entry["primaryAccession"]: entry for entry in r.json().get("results", [])}


def sequence(entry):
    return (entry.get("sequence") or {}).get("value", "")


def _species(name):
    return " ".join(name.split()[:2]).lower()


def organism(entry):
    """UniProt's organism without the strain ('Aspergillus fumigatus (strain ATCC MYA-4609 / ...
    / Af293)' -> 'Aspergillus fumigatus'); CAZy's species name instead when UniProt's name, one
    of its synonyms or its "common name" (which holds other scientific names) is a species of
    CAZy's organisms, so the catalogue has one organism per species (UniProt's Emericella
    nidulans -> CAZy's Aspergillus nidulans)."""
    info = entry.get("organism") or {}
    name = re.sub(r"\s*\(.*\)\s*$", "", info.get("scientificName", "")).strip()
    for other in [name] + (info.get("synonyms") or []) + [info.get("commonName") or ""]:
        if _species(other) in _cazy_species:
            return _cazy_species[_species(other)]
    return name


def kingdom(entry):
    """The top of the lineage, as CAZy's kingdom column has it (Bacteria, Eukaryota, ...)."""
    lineage = (entry.get("organism") or {}).get("lineage") or []
    return lineage[0] if lineage else ""


def protein_name(entry, extra=""):
    """The names as UniProt lists them, 'Full name (Short) (Other full) (Other short)', so that
    update_protein_codes finds a short name such as AfAA9C; `extra` adds a literature name."""
    description = entry.get("proteinDescription", {})
    blocks = ([description.get("recommendedName") or {}] + description.get("alternativeNames", [])
              + description.get("submissionNames", []))
    names = []
    for block in blocks:
        names.append((block.get("fullName") or {}).get("value", ""))
        names += [s.get("value", "") for s in block.get("shortNames", [])]
    names = [n for n in names if n]
    if extra and extra not in names:
        names.append(extra)
    return (names[0] + "".join(f" ({n})" for n in names[1:])) if names else entry["primaryAccession"]


def genbank_ids(entry):
    return sorted({p["value"] for x in entry.get("uniProtKBCrossReferences", []) if x["database"] == "EMBL"
                   for p in x.get("properties", []) if p["key"] == "ProteinId" and p["value"] != "-"})


def pdb_codes(entry):
    return sorted({x["id"] for x in entry.get("uniProtKBCrossReferences", []) if x["database"] == "PDB"})


def curated_family(entry):
    """The family a reviewed entry's similarity comment names, or None."""
    if "Swiss-Prot" not in entry.get("entryType", ""):
        return None
    for comment in entry.get("comments", []):
        if comment["commentType"] == "SIMILARITY":
            for text in comment.get("texts", []):
                m = re.search(r"polysaccharide monooxygenase (AA\d+) family", text["value"])
                if m:
                    return m.group(1)
    return None


def experimental_activities(entry):
    """[(EC number, [PubMed IDs])] of the entry's experimentally shown catalytic activities."""
    found = []
    for comment in entry.get("comments", []):
        if comment["commentType"] != "CATALYTIC ACTIVITY":
            continue
        reaction = comment["reaction"]
        pmids = sorted({e["id"] for e in reaction.get("evidences", [])
                        if e.get("evidenceCode") == EXPERIMENTAL and e.get("source") == "PubMed"})
        if reaction.get("ecNumber") in LPMO_EC_NUMBERS and pmids:
            found.append((reaction["ecNumber"], pmids))
    return found


def pubmed_links(pmids):
    return "; ".join(f"http://pubmed.ncbi.nlm.nih.gov/{p}" for p in pmids)


# ---------------------------------------------------------------------------
# Family by similarity to the CAZy-listed proteins
# ---------------------------------------------------------------------------

class FamilyFinder:
    """The family of the most similar CAZy-listed protein (see the module docstring)."""

    def __init__(self, references):
        self.references = [(family, seq) for family, seq in references if seq]
        self.aligner = Align.PairwiseAligner(mode="local", open_gap_score=-10, extend_gap_score=-0.5,
                                             substitution_matrix=substitution_matrices.load("BLOSUM62"))

    def identity(self, seq, ref):
        """Identity and aligned length of the best local alignment."""
        counts = self.aligner.align(seq, ref)[0].counts()
        aligned = counts.identities + counts.mismatches
        return (counts.identities / aligned if aligned else 0.0), aligned

    def family(self, seq):
        """(family, identity, note): the family, or None with the reason in the note."""
        if not seq or not self.references:
            return None, 0.0, "no sequence"
        # the three best-scoring references of each family, then their identities
        by_family = defaultdict(list)
        for family, ref in self.references:
            by_family[family].append((self.aligner.score(seq, ref), ref))
        best = {}
        for family, scored in by_family.items():
            for _, ref in sorted(scored, key=lambda x: -x[0])[:3]:
                identity, aligned = self.identity(seq, ref)
                if aligned >= MIN_ALIGNED and identity > best.get(family, 0.0):
                    best[family] = identity
        ranked = sorted(best.items(), key=lambda x: -x[1])
        if not ranked:
            return None, 0.0, "no alignment of 100 residues to any family"
        (family, identity), runner = ranked[0], (ranked[1] if len(ranked) > 1 else (None, 0.0))
        note = f"{family} {identity:.0%}" + (f", {runner[0]} {runner[1]:.0%}" if runner[0] else "")
        if identity < MIN_IDENTITY:
            return None, identity, f"best {note}: below {MIN_IDENTITY:.0%}"
        if identity - runner[1] < MIN_MARGIN:
            return None, identity, f"best {note}: families too close"
        return family, identity, note

    def same_protein(self, seq, candidates):
        """The first (accession, identity) of `candidates` [(accession, sequence)] that `seq` matches
        at >= SAME_PROTEIN over most of its length, or None."""
        for acc, ref in candidates:
            identity, aligned = self.identity(seq, ref)
            if identity >= SAME_PROTEIN and aligned >= 0.8 * min(len(seq), len(ref)):
                return acc, identity
        return None


# ---------------------------------------------------------------------------
# The three lists
# ---------------------------------------------------------------------------

def _row(entry, family, list_source, ec_numbers, pmids, evidence, extra_name="", variant_of=""):
    return {
        "family": family, "kingdom": kingdom(entry), "protein_name": protein_name(entry, extra_name),
        "ec_numbers": "; ".join(ec_numbers), "organism": organism(entry),
        "genbank_ids": "; ".join(genbank_ids(entry)), "uniprot_ids": entry["primaryAccession"],
        "pdb_codes": ", ".join(pdb_codes(entry)), "references": pubmed_links(pmids),
        "in_cazy_characterized": "no", "in_cazy_structure": "no",
        "list_source": list_source, "evidence": evidence, "variant_of": variant_of,
    }


_paper_texts = {}


def paper_text(pmid):
    """Title and abstract of a PubMed paper (Europe PMC); None when the lookup fails.

    None is not the same as "": "" means Europe PMC holds no text for the paper, None that the
    request did not get through. The caller must not read None as "the paper does not name the
    organism" -- on 2026-09-23 one request timed out and TgAA11, whose paper does name
    Trichoderma guizhouense in its title, was dropped from the catalogue without a word. Only
    successful lookups are cached, so a later paper of the same entry can still be checked."""
    if pmid not in _paper_texts:
        try:
            r = _session.get("https://www.ebi.ac.uk/europepmc/webservices/rest/search", timeout=60,
                             params={"query": f"EXT_ID:{pmid} AND SRC:MED", "resultType": "core",
                                     "format": "json"})
            r.raise_for_status()
            hits = r.json().get("resultList", {}).get("result", [])
            _paper_texts[pmid] = " ".join((hits[0].get("title") or "", hits[0].get("abstractText") or "")) if hits else ""
        except (requests.RequestException, ValueError) as e:
            print(f"  WARNING: Europe PMC lookup of PubMed {pmid} failed ({e.__class__.__name__})")
            return None
    return _paper_texts[pmid]


def genera(entry):
    """The genus names UniProt gives the organism (scientific name, synonyms, "common name")."""
    info = entry.get("organism") or {}
    names = [info.get("scientificName") or ""] + (info.get("synonyms") or []) + [info.get("commonName") or ""]
    return {n.split()[0] for n in names if n and n[0].isupper()}


def brenda_candidates(path):
    """{accession: {"ecs": set, "pmids": set}} of BRENDA's LPMO proteins with a UniProt accession
    and a substrate row citing a paper."""
    from activities import load_brenda, NOT_SUBSTRATE

    release, entries = load_brenda(path)
    found = defaultdict(lambda: {"ecs": set(), "pmids": set()})
    for ec, entry in entries.items():
        refs = entry.get("reference", {})
        for row in entry.get("substrates_products", []):
            substrate = row.get("value", "").split(" = ")[0].split(" + ")[0].strip()
            if NOT_SUBSTRATE.match(substrate):
                continue
            pmids = {str(refs[r]["pmid"]) for r in row.get("references", [])
                     if r in refs and refs[r].get("pmid")}
            if not pmids:
                continue
            for pid in row.get("proteins", []):
                for acc in entry["protein"].get(pid, {}).get("accessions") or []:
                    found[acc]["ecs"].add(ec)
                    found[acc]["pmids"] |= pmids
    return release, dict(found)


def _local_model_rows():
    """1e: curation/local_models.csv as dicts, with the optional `accession` column (the
    UniProt entry of the protein the model stands in for) defaulted to ""."""
    if not os.path.exists(LOCAL_MODELS_CSV):
        return []
    local = pd.read_csv(LOCAL_MODELS_CSV, dtype=str, keep_default_na=False, comment="#")
    if "accession" not in local.columns:
        local["accession"] = ""
    return [{k: str(v).strip() for k, v in row.items()} for row in local.to_dict("records")]


def collect(cazy, brenda_path=None, curated_path=CURATED_CSV, log=print):
    """The other_lpmos.csv rows (see the module docstring) for the CAZy table `cazy`."""
    _cazy_species.clear()
    for name in cazy.organism.dropna():
        words = name.split()
        if len(words) >= 2 and words[1] not in ("sp.", "cf."):
            _cazy_species.setdefault(_species(name), " ".join(words[:2]))
    cazy_accessions = {a.strip() for v in cazy.uniprot_ids.fillna("") for a in re.split(r"[;,\s]+", v)
                       if a.strip()}
    family_of_accession = {}
    for r in cazy.itertuples():
        for a in re.split(r"[;,\s]+", str(r.uniprot_ids) if pd.notna(r.uniprot_ids) else ""):
            if a.strip():
                family_of_accession.setdefault(a.strip(), r.family)
    cazy_pdbs = {c.strip().upper() for v in cazy.pdb_codes.fillna("") for c in v.split(",") if c.strip()}
    cazy_genbank = {g.strip() for v in cazy.genbank_ids.fillna("") for g in v.split(";") if g.strip()}
    cazy_entries = uniprot_entries(cazy_accessions)
    cazy_sequences = {sequence(e): acc for acc, e in cazy_entries.items() if sequence(e)}
    by_species = defaultdict(list)       # species -> [(accession, sequence)] of CAZy-listed entries
    for acc, e in cazy_entries.items():
        by_species[" ".join(organism(e).split()[:2])].append((acc, sequence(e)))
    finder = FamilyFinder([(family_of_accession[acc], sequence(e)) for acc, e in cazy_entries.items()])
    # entries CAZy lists by a deleted UniProt accession or by a GenBank protein alone: their
    # sequence from UniParc, so an identical reviewed entry (G2QZK6 for the deleted A0A446B9M7)
    # counts as listed
    unsequenced = 0
    for r in cazy.itertuples():
        ids = [a.strip() for a in re.split(r"[;,\s]+", str(r.uniprot_ids) if pd.notna(r.uniprot_ids) else "")
               if a.strip()]
        if any(a in cazy_entries for a in ids):
            continue
        genbank = [g.strip() for g in str(r.genbank_ids).split(";") if g.strip()] if pd.notna(r.genbank_ids) else []
        for db_id in ids + genbank:
            seq = uniparc_sequence(db_id)
            if seq:
                cazy_sequences.setdefault(seq, db_id)
                by_species[" ".join(str(r.organism).split()[:2])].append((db_id, seq))
                break
        else:
            unsequenced += 1
    log(f"  {len(cazy_entries)} active UniProt entries of CAZy's lists as the family references; "
        f"{len(cazy_sequences) - len(cazy_entries)} more CAZy sequences from UniParc, {unsequenced} entries without one")

    rows, taken, dropped, variants = [], set(), defaultdict(list), []
    # 1e: a local model can name the UniProt accession of the protein it models, so that the
    # entry keeps its UniProt annotations (EC numbers, papers, activities). The accession is
    # taken here, before 1b-1d, so the Swiss-Prot/BRENDA route does not add the protein a
    # second time -- the model replaces the database structure, not the entry.
    local_models = _local_model_rows()
    local_accessions = {r["accession"] for r in local_models if r["accession"]}
    taken |= local_accessions
    added_sequences, added_by_species = {}, defaultdict(list)

    def admit(entry, family_hint, list_source, ec_numbers, pmids, evidence, extra_name="", paper_filter=None):
        """Add the entry unless it is there already; paper_filter(entry, pmids) keeps the papers
        that count, and {pmids} in the evidence is replaced with them."""
        acc = entry["primaryAccession"]
        seq = sequence(entry)
        if acc in taken:
            return
        if acc in cazy_accessions:
            dropped["listed by CAZy"].append(acc)
            return
        shared_pdb = sorted(set(pdb_codes(entry)) & cazy_pdbs)
        if shared_pdb:
            dropped["structure CAZy lists"].append(f"{acc} ({', '.join(shared_pdb)})")
            return
        shared_genbank = sorted(set(genbank_ids(entry)) & cazy_genbank)
        if shared_genbank:
            dropped["GenBank protein CAZy lists"].append(f"{acc} ({', '.join(shared_genbank)})")
            return
        if seq in cazy_sequences:
            dropped["sequence of a CAZy-listed entry"].append(f"{acc} (= {cazy_sequences[seq]})")
            return
        if seq in added_sequences:
            dropped["sequence of an entry already added"].append(f"{acc} (= {added_sequences[seq]})")
            return
        species = " ".join(organism(entry).split()[:2])
        same = finder.same_protein(seq, added_by_species[species])
        if same:
            dropped["strain variant of an entry already added"].append(f"{acc} ({same[1]:.0%} = {same[0]})")
            return
        if paper_filter:
            kept = paper_filter(entry, pmids)
            if not kept:
                dropped[f"{list_source} papers do not name the organism"].append(
                    f"{acc} {organism(entry)} (PubMed {', '.join(pmids)})")
                return
            pmids, evidence = kept, evidence.replace("{pmids}", ", ".join(kept))
        variant = finder.same_protein(seq, by_species.get(species, []))
        if variant:
            variants.append(f"{acc} ({variant[1]:.0%} = {variant[0]})")
        family = family_hint
        how = "curated" if family_hint else None
        if not family:
            family = curated_family(entry)
            how = "Swiss-Prot similarity comment" if family else None
        if not family:
            family, identity, note = finder.family(seq)
            how = f"most similar CAZy-listed protein ({note})" if family else None
            if not family:
                dropped["family not clear - curate it if it is an LPMO"].append(
                    f"{acc} {organism(entry)} ({note})")
                return
        if family not in FAMILIES:
            dropped["family outside AA9-AA17"].append(f"{acc} ({family})")
            return
        taken.add(acc)
        added_sequences[seq] = acc
        added_by_species[species].append((acc, seq))
        rows.append(_row(entry, family, list_source, ec_numbers, pmids, f"{evidence}; family: {how}", extra_name,
                         variant[0] if variant else ""))

    # 1b Swiss-Prot
    swissprot = swissprot_lpmo_entries()
    shown = []
    for acc, entry in swissprot.items():
        activities = experimental_activities(entry)
        if not activities:
            dropped["Swiss-Prot activity not experimentally shown"].append(acc)
            continue
        ecs = sorted({ec for ec, _ in activities})
        pmids = sorted({p for _, ps in activities for p in ps})
        shown.append((acc, entry, ecs, pmids))
    # the entry citing more papers first, of two with the same protein
    for acc, entry, ecs, pmids in sorted(shown, key=lambda x: (-len(x[3]), x[0])):
        admit(entry, None, "Swiss-Prot", ecs, pmids,
              f"Swiss-Prot: EC {', '.join(ecs)} experimentally shown (PubMed {', '.join(pmids)})")
    log(f"  1b Swiss-Prot: {len(swissprot)} reviewed entries with an LPMO EC number, "
        f"{sum(r['list_source'] == 'Swiss-Prot' for r in rows)} added")

    # 1c BRENDA
    if brenda_path and os.path.exists(brenda_path):
        release, candidates = brenda_candidates(brenda_path)
        entries = uniprot_entries(candidates)
        missing = sorted(set(candidates) - set(entries))
        if missing:
            dropped["deleted UniProt entry"] += missing
        unchecked = []

        def naming_papers(entry, pmids):
            """The papers whose title or abstract names a genus of the entry's organism. A paper
            Europe PMC did not answer for counts as naming it: a failed request is not evidence
            against the protein, and dropping it would lose it from the catalogue silently."""
            kept = []
            for p in pmids:
                text = paper_text(p)
                if text is None:
                    unchecked.append(f"{p} ({entry['primaryAccession']})")
                    kept.append(p)
                elif any(re.search(rf"\b{re.escape(g)}\b", text) for g in genera(entry)):
                    kept.append(p)
            return kept

        for acc in sorted(entries):
            ecs, pmids = sorted(candidates[acc]["ecs"]), sorted(candidates[acc]["pmids"])
            admit(entries[acc], None, "BRENDA", ecs, pmids,
                  f"BRENDA {release}: EC {', '.join(ecs)}, substrate rows citing PubMed {{pmids}}",
                  paper_filter=naming_papers)
        log(f"  1c BRENDA {release}: {len(candidates)} proteins with a UniProt accession and a cited "
            f"substrate row, {sum(r['list_source'] == 'BRENDA' for r in rows)} added")
        if unchecked:
            log(f"  1c WARNING: Europe PMC did not answer for {len(unchecked)} paper(s), kept as "
                f"naming the organism - check them: {', '.join(dict.fromkeys(unchecked))}")
    else:
        log("  1c BRENDA: no release file (data/brenda_*.json.tar.gz or BRENDA_JSON) - skipped")

    # 1d curated exceptions
    if curated_path and os.path.exists(curated_path):
        curated = pd.read_csv(curated_path, dtype=str, keep_default_na=False, comment="#")
        entries = uniprot_entries(curated.accession)
        for r in curated.itertuples():
            if r.accession not in entries:
                dropped["deleted UniProt entry"].append(r.accession)
                continue
            admit(entries[r.accession], r.family, "curated", [], [r.pmid] if r.pmid else [],
                  f"curated: {r.reason}", extra_name=r.protein)
        log(f"  1d curated: {len(curated)} entries, {sum(r['list_source'] == 'curated' for r in rows)} added")

    # 1e local models
    for r in local_models:
        if not os.path.exists(os.path.join(PIPELINE_DIR, r["file"])):
            dropped["local model file missing"].append(f'{r["struct_id"]} ({r["file"]})')
            continue
        rows.append({"family": r["family"], "kingdom": r["kingdom"], "protein_name": r["protein"],
                     "ec_numbers": "", "organism": r["organism"], "genbank_ids": "",
                     "uniprot_ids": r["accession"], "pdb_codes": "",
                     "references": pubmed_links([r["pmid"]] if r["pmid"] else []),
                     "in_cazy_characterized": "no", "in_cazy_structure": "no",
                     "list_source": "local model", "evidence": f'local model: {r["reason"]}',
                     "variant_of": "", "struct_id": r["struct_id"], "local_file": r["file"]})
    if local_models:
        log(f"  1e local models: {len(local_models)} entries, "
            f"{sum(r['list_source'] == 'local model' for r in rows)} added"
            + (f" (UniProt entries they replace: {', '.join(sorted(local_accessions))})"
               if local_accessions else ""))

    for reason, accs in sorted(dropped.items()):
        accs = list(dict.fromkeys(accs))     # an accession Swiss-Prot and BRENDA both give, once
        log(f"  dropped, {reason} ({len(accs)}): {', '.join(accs)}")
    if variants:
        log(f"  strain variants of CAZy-listed proteins, used in step 2 only if the CAZy entry gets no "
            f"structure ({len(variants)}): {', '.join(variants)}")
    return pd.DataFrame(rows, columns=COLUMNS)


def main():
    from activities import find_brenda_file

    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--data-dir", default=DEFAULT_DATA_DIR)
    parser.add_argument("--brenda", help="BRENDA JSON release (.json or .json.tar.gz)")
    args = parser.parse_args()
    cazy = pd.read_csv(os.path.join(args.data_dir, "01_proteins", "cazy_lpmos.csv"), dtype=str)
    print("1b-1d: characterized LPMOs outside CAZy's lists (Swiss-Prot, BRENDA, curated)...")
    other = collect(cazy, args.brenda or find_brenda_file(args.data_dir))
    out = os.path.join(args.data_dir, "01_proteins", "other_lpmos.csv")
    other.to_csv(out, index=False)
    print(f"Saved {len(other)} entries to {out}")


if __name__ == "__main__":
    main()
