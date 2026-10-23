"""
2j: substrates and oxidised positions of the selected proteins, each with the papers that show
it, from three curated sources:

  CAZy        the EC numbers of the protein's CAZy entries; CAZy lists its papers per entry,
              so they are not tied to one EC number
  Swiss-Prot  the reviewed UniProt entries holding the protein's sequence: catalytic
              activities with experimental evidence (ECO:0000269) and the sentences of the
              function text, each with its PubMed IDs (substrates, substrates without
              activity, C1/C4); other experimental comments and features (pH/temperature,
              kinetics, regulation, PTMs, mutants) as notes
  BRENDA      the proteins of EC 1.14.99.53-56 in a downloaded BRENDA JSON release (CC BY
              4.0; brenda-enzymes.org/download.php, saved as data/brenda_*.json.tar.gz
              or given with --brenda / BRENDA_JSON): substrate -> product rows, "no substrate"
              remarks and notes, each with its papers. A BRENDA protein is the catalogue
              protein holding its UniProt accession; a BRENDA protein without an accession
              (BRENDA pools the enzymes of an organism) is matched per row, to the only
              catalogue protein of its species citing the row's paper (among several, the one
              whose name holds a BRENDA synonym of that paper)

  Literature  hand-curated claims for what none of the three holds, one row per claim with its
              primary paper and where the claim was checked (lpmocat_pipeline/curation/
              activities.csv: struct_id, protein, substrate, form, result, oxidation, pmid, doi,
              evidence); reviews are never cited

The collector does not decide between sources: EC numbers are kept as EC numbers (their
C1/C4 reading is left to populate_db), and every claim keeps its source, text and papers.

Output (02_data/activities.json):
    {"sources": {"brenda": {"file", "release"}, "checked"},
     "proteins": {struct_id: {
         "ec": {ec: {source: [ref, ...]}},
         "observations": [{"source", "substrate", "form", "result": "active" | "inactive",
                           "oxidation": ["C1", "C4"], "only": true when the text excludes the
                           other position, "refs", "text", "accession" | "brenda_protein",
                           "match"}],
         "substrates": {substrate: {"result": "active" | "inactive" | "conflict",
                                    "forms": [...], "active": [ref, ...], "inactive": [ref, ...],
                                    "oxidation": {"C1": [ref, ...], "C4": [ref, ...]},
                                    "sources": [...]}},
         "notes": [{"source", "topic", "text", "refs"}]}},
     "references": {ref: {"title", "author", "n_authors", "year", "journal", "doi", "pmid"}},
     "unmatched": {"brenda": [{"ec", "protein", "organism", "accessions", "reason"}]},
     "unparsed": [substrate texts no category matched]}
    A ref is "PMID:<id>", "DOI:<doi>" or "BRENDA:<ec>:<reference id>" (a BRENDA paper
    without a PubMed ID). A substrate observation without a substrate (Swiss-Prot's
    "yielding C1 oxidation products") has substrate null.

activities_review.csv holds the same claims as one row each (protein, substrate, form, result,
oxidised positions, source, papers, the sentence or reaction it comes from, and the flag of a
substrate whose sources disagree), for reading through and curating.

Runs as part of pipeline step 2 (fetch); standalone it fills activities.json and
activities_review.csv for the current pdb_selection.csv:
    python lpmocat_pipeline/steps/activities.py [--data-dir DIR] [--brenda FILE] [--no-csv]
"""

import argparse
import glob
import io
import json
import os
import re
import tarfile
import time
from collections import Counter, defaultdict
from datetime import datetime

import pandas as pd
import requests

from paths import CURATION_DIR, DEFAULT_DATA_DIR, PIPELINE_DIR
LPMO_EC_NUMBERS = ("1.14.99.53", "1.14.99.54", "1.14.99.55", "1.14.99.56")
EC_SUBSTRATE = {"1.14.99.53": "chitin", "1.14.99.54": "cellulose", "1.14.99.55": "starch",
                "1.14.99.56": "cellulose"}
EXPERIMENTAL = "ECO:0000269"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb"
EUROPEPMC_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

# Substrate names -> (category, form). The first pattern that matches wins, so oligosaccharides
# and derivatives come before their polymers. The form is the matched name as written.
SUBSTRATE_PATTERNS = [
    (r"chito-?oligo\w*|chito(?:bi|tri|tetra|penta|hexa)ose|N-acetyl-beta-D-glucosaminyl\]\d|"
     r"\(GlcNAc\)\d|GlcNAc\d", "chito-oligosaccharides"),
    (r"chitosan", "chitosan"),
    (r"N-acetyl-beta-D-glucosaminyl\]n|\[\(1->4\)-N-acetyl|chitin", "chitin"),
    (r"cello-?oligo\w*|cello-?(?:bi|tri|tetra|penta|hexa)os\w*|gluco-oligosacch\w*", "cello-oligosaccharides"),
    (r"xylo-?glucan[ -]oligo\w*|xyloglucan oligomer\w*|XXXG\w*|xylosyl\]\d", "xyloglucan oligosaccharides"),
    (r"xyloglucan", "xyloglucan"),
    (r"glucomannan", "glucomannan"),
    (r"(?:galacto)?mannan|mannohexaose|manno-?oligo\w*", "mannan"),
    (r"xylo-?oligo\w*|xylo(?:bi|tri|tetra|penta|hexa)ose", "xylo-oligosaccharides"),
    (r"(?:arabino|glucurono|arabinoglucurono)?xylan", "xylan"),
    (r"mixed[- ]linkage[ -](?:beta-)?glucan|beta-\(1->3, ?1->4\)-glucan|\(1,3;1,4\)-beta-glucan|"
     r"lichenan|barley beta-glucan|soluble beta-glucan|beta-glucan", "mixed-linkage glucan"),
    (r"malto-?oligo\w*|malto(?:tri|tetra|penta|hexa|hepta)ose", "malto-oligosaccharides"),
    (r"starch|amylose|amylopectin|glycogen", "starch"),
    (r"pectin|homogalacturonan|polygalacturon\w*", "pectin"),
    (r"laminarin|beta-1,3-glucan|curdlan", "beta-1,3-glucan"),
    (r"carboxy ?methyl ?c?ellulose|CMC|cellulose acetate|hydroxyethyl cellulose|methyl cellulose",
     "cellulose derivatives"),
    (r"phosphoric acid[- ]swollen cellul\w*|PASC|avicel\w*|microcrystalline cellulose\w*|"
     r"amorphous cellulose|regenerated (?:amorphous )?cellulose|bacterial (?:microcrystalline )?cellulose|"
     r"cellulose nanofibrils?|nanocellulose|filter paper|cotton|kraft pulp|pulp fib\w*|"
     r"crystalline (?:forms of )?cellulose|\[\(1->4\)-beta-D-glucosyl\]n|cellulose", "cellulose"),
    (r"spruce|straw|soy spent flakes|lignocellulos\w*|biomass|bagasse|corn stover", "biomass"),
    (r"galactan", "galactan"),
    (r"peptidoglycan", "peptidoglycan"),
    (r"alginate", "alginate"),
]
# a name starts at a word start and is not the enzyme ("chitinase", "xylanolytic") or a
# qualifier ("cellulose-associated xylan")
# qualifiers kept with the form ("alpha-chitin", "konjac glucomannan")
QUALIFIER = (r"(?:(?:alpha|beta|colloidal|squid[- ]pen|crystalline|amorphous|soluble|insoluble|reduced|"
             r"konjac|tamarind|birch ?wood|beech ?wood|oat spelt|wheat|barley|ivory nut|corn|potato|carob|"
             r"guar|spruce)[- ])*")
_SUBSTRATE_RES = [(re.compile(rf"(?<![A-Za-z0-9]){QUALIFIER}(?:{pattern})"
                              rf"(?!ases?\b|olytic|-(?:associated|binding|bound))", re.I), category)
                  for pattern, category in SUBSTRATE_PATTERNS]
# names that are not polysaccharide substrates (assay compounds, cosubstrates)
NOT_SUBSTRATE = re.compile(r"^(?:more|2,6-dimethoxyphenol|H2O2|O2|ascorbate|ascorbic acid|"
                           r"hydrogen peroxide|Amplex red)$", re.I)


# ---------------------------------------------------------------------------
# Text helpers
# ---------------------------------------------------------------------------

def substrates_in(text):
    """[(category, form)] of the substrate names in a text, in order, one per category."""
    found, seen = [], set()
    spans = []
    for pattern, category in _SUBSTRATE_RES:
        for m in pattern.finditer(text):
            if any(a <= m.start() < b for a, b in spans):
                continue
            spans.append((m.start(), m.end()))
            if category not in seen:
                seen.add(category)
                found.append((m.start(), category, m.group(0)))
    return [(category, form) for _, category, form in sorted(found)]


def oxidation_in(text):
    """(positions, only): the oxidised carbons a text names, and whether it says the enzyme
    oxidises only there."""
    t = text.replace("‐", "-")
    # "C1 or C4" names no position: it is UniProt's family template ("... yielding C1 or C4
    # oxidation products"), used even where the cited paper found only C1 (AoAA11, 4MAI)
    t = re.sub(r"\bC-?1\s+or\s+C-?4\b|\bC-?4\s+or\s+C-?1\b", " ", t)
    positions = set()
    if re.search(r"C1\s*(?:/|and|,)\s*C4|C4\s*(?:/|and)\s*C1|C1- and C4-|C1 and C4|"
                 r"double \(?C4/C1\)?|both C1 and C4|mixed C1/C4", t):
        positions |= {"C1", "C4"}
    if re.search(r"\bC-?1[- ](?:oxidi|hydroxylat|oxidation)|\bC-?1\b(?! and| or|/)|aldonic acid|"
                 r"glucono-1,5-lactone|reducing[- ]end (?:oxidation|lactone)", t):
        positions.add("C1")
    if re.search(r"\bC-?4[- ](?:oxidi|dehydrogenat|oxidation)|\bC-?4\b|4-keto|4-dehydro|"
                 r"non-?reducing[- ]end (?:oxidation)|gem-?diol", t):
        positions.add("C4")
    if re.search(r"\bC-?6[- ](?:oxidi|aldehyde)", t):
        positions.add("C6")
    only = bool(re.search(r"\b(?:only|exclusively|strictly|solely)\s+(?:C-?[14]|yields? C-?[14])", t))
    return sorted(positions), only


NEGATION = re.compile(r"^(?:not |no |does not|do not|is not|shows? no|displays? no|lacks?|"
                      r"inactive|unable|without activity|cannot)|\bnot active|\bno (?:activity|"
                      r"detectable|oxidative|substrates?)\b|\binactive (?:on|towards|against)|"
                      r"\bdoes(?: not|n't) (?:cleave|oxidi|degrade|act|show)|\bunable to|"
                      r"\bno products?\b|\bnot (?:cleave|oxidi|degrade)", re.I)


def clauses(sentence):
    """Split a sentence where a positive and a negative statement meet ("... but not ...")."""
    parts = re.split(r",? but (?=not |no |shows no |is not |does not )|; whereas |, whereas |; ",
                     sentence)
    return [p.strip() for p in parts if p.strip()]


def sentences(text):
    """The sentences of a UniProt comment, split after the closing '(PubMed:...).'"""
    parts = re.split(r"(?<=[.)])\s+(?=[A-Z(])", text)
    return [p.strip() for p in parts if p.strip()]


def pubmed_refs(text):
    return [f"PMID:{pmid}" for pmid in re.findall(r"PubMed:(\d+)", text)]


def valid_pmid(pmid):
    """True for a plausible PubMed ID: 1-8 digits, no leading zero. PubMed assigns them in
    sequence (~40 million in 2026), so a longer number is a typo, not a future ID. BRENDA
    2026.1 carries one (315404984 for PMID 31540984, EC 1.14.99.54 ref 49 and .56 ref 25),
    which split that paper into two references, one of them a dead PubMed link."""
    pmid = str(pmid or "").strip()
    return bool(re.fullmatch(r"[1-9]\d{0,7}", pmid))


def _ref_from_url(url):
    """CAZy reference links -> 'PMID:...' / 'DOI:...'."""
    url = (url or "").strip()
    m = re.search(r"pubmed[^/]*/(?:\?term=)?(\d+)", url) or re.search(r"list_uids=(\d+)", url)
    if m:
        return f"PMID:{m.group(1)}"
    m = re.search(r"doi\.org/(10\..+)$", url)
    if m:
        return f"DOI:{m.group(1).lower()}"
    return None


def _ids(value):
    return [v.strip() for v in re.split(r"[;,|]", str(value or "")) if v.strip() and v.strip() != "nan"]


def _species(name):
    words = (name or "").replace("[", "").replace("]", "").split()
    if len(words) >= 3 and words[1] in ("sp.", "cf."):
        return " ".join(words[:3]).lower()
    return " ".join(words[:2]).lower()


def _add(mapping, key, refs):
    target = mapping.setdefault(key, [])
    for ref in refs:
        if ref not in target:
            target.append(ref)


# ---------------------------------------------------------------------------
# CAZy
# ---------------------------------------------------------------------------

def cazy_rows_of(selection, cazy):
    """{struct_id: [CAZy rows]}: the entry of the selection row, the entries 2c marked as its
    duplicates, and entries sharing a UniProt accession, GenBank protein or PDB code with it
    (e.g. the structure entry of a protein 2e dropped as a sequence duplicate)."""
    by_key = {(r.family, r.protein_name, r.organism): r for r in cazy.itertuples()}
    rows = defaultdict(list)
    kept = selection[selection.source != "NONE"]
    owner_of = {}
    for r in kept.itertuples():
        for acc in _ids(r.uniprot_ids):
            owner_of.setdefault(("uniprot", acc), r.selected_pdb)
        for gb in _ids(getattr(r, "genbank_id", "")):
            owner_of.setdefault(("genbank", gb), r.selected_pdb)
        for code in _ids(r.all_pdb_codes) + [r.selected_pdb]:
            owner_of.setdefault(("pdb", code.upper()), r.selected_pdb)
    for r in kept.itertuples():
        entry = by_key.get((r.family, r.protein_name, r.organism))
        if entry is not None:
            rows[r.selected_pdb].append(entry)
    names = {(r.family, r.protein_name, r.organism): r.selected_pdb for r in kept.itertuples()}
    for r in selection[selection.source == "NONE"].itertuples():
        if getattr(r, "duplicate_of", ""):
            # duplicate_of is "<family> <protein name>" of the entry that holds the structure
            owner = next((s for (fam, name, _), s in names.items()
                          if fam == r.family and r.duplicate_of in (name, f"{fam} {name}")), None)
            entry = by_key.get((r.family, r.protein_name, r.organism))
            if owner and entry is not None:
                rows[owner].append(entry)
    listed = {id(e) for entries in rows.values() for e in entries}
    for entry in cazy.itertuples():
        if id(by_key.get((entry.family, entry.protein_name, entry.organism))) in listed:
            continue
        keys = ([("uniprot", a) for a in _ids(entry.uniprot_ids)]
                + [("genbank", g) for g in _ids(entry.genbank_ids)]
                + [("pdb", p.upper()) for p in _ids(entry.pdb_codes)])
        owners = {owner_of[k] for k in keys if k in owner_of}
        if len(owners) == 1:
            rows[owners.pop()].append(entry)
    return rows


def collect_cazy(protein, entries):
    for entry in entries:
        refs = [ref for ref in map(_ref_from_url, _ids(entry.references)) if ref]
        for ec in _ids(entry.ec_numbers):
            _add(protein["ec"].setdefault(ec, {}), "CAZy", refs)
            if ec in EC_SUBSTRATE:
                protein["observations"].append({
                    "source": "CAZy", "substrate": EC_SUBSTRATE[ec], "form": None,
                    "result": "active", "oxidation": [], "only": False, "refs": refs,
                    "text": f"EC {ec} ({entry.protein_name})", "ec": ec})


# ---------------------------------------------------------------------------
# Swiss-Prot
# ---------------------------------------------------------------------------

_uniprot_session = requests.Session()
_uniprot_session.mount("https://", requests.adapters.HTTPAdapter(max_retries=requests.adapters.Retry(
    total=5, backoff_factor=1, status_forcelist=(429, 500, 502, 503, 504))))


def fetch_uniprot_entries(accessions, batch=50):
    """{accession: UniProtKB JSON entry} of active entries."""
    entries = {}
    accessions = sorted(set(accessions))
    for i in range(0, len(accessions), batch):
        chunk = accessions[i:i + batch]
        r = _uniprot_session.get(f"{UNIPROT_URL}/accessions",
                                 params={"accessions": ",".join(chunk), "format": "json"}, timeout=120)
        r.raise_for_status()
        for entry in r.json().get("results", []):
            entries[entry["primaryAccession"]] = entry
    return entries


def _evidence_refs(evidences, experimental_only=True):
    return [f"PMID:{e['id']}" for e in (evidences or [])
            if e.get("source") == "PubMed" and (not experimental_only or e.get("evidenceCode") == EXPERIMENTAL)]


def _is_experimental(evidences):
    return any(e.get("evidenceCode") == EXPERIMENTAL for e in (evidences or []))


NOTE_COMMENTS = ("ACTIVITY REGULATION", "BIOPHYSICOCHEMICAL PROPERTIES", "PTM", "INDUCTION",
                 "DISRUPTION PHENOTYPE", "BIOTECHNOLOGY", "COFACTOR", "DOMAIN", "SUBUNIT",
                 "PATHWAY")
NOTE_FEATURES = ("Mutagenesis", "Modified residue", "Glycosylation", "Disulfide bond", "Binding site",
                 "Active site", "Site")


def collect_swissprot(protein, accession, entry, references):
    """Catalytic activities, function sentences and notes of one reviewed entry."""
    for ref in entry.get("references", []):
        citation = ref.get("citation", {})
        xrefs = {x["database"]: x["id"] for x in citation.get("citationCrossReferences", [])}
        if xrefs.get("PubMed"):
            authors = citation.get("authors") or []
            references.setdefault(f"PMID:{xrefs['PubMed']}", {
                "pmid": xrefs["PubMed"], "doi": xrefs.get("DOI"), "title": citation.get("title"),
                "author": authors[0] if authors else None, "n_authors": len(authors),
                "year": (citation.get("publicationDate") or "")[:4] or None,
                "journal": citation.get("journal")})

    for comment in entry.get("comments", []):
        kind = comment["commentType"]
        if kind == "CATALYTIC ACTIVITY":
            reaction = comment["reaction"]
            ec = reaction.get("ecNumber")
            if ec and _is_experimental(reaction.get("evidences")):
                _add(protein["ec"].setdefault(ec, {}), "Swiss-Prot", _evidence_refs(reaction.get("evidences")))
        elif kind == "FUNCTION":
            for text in comment.get("texts", []):
                if not _is_experimental(text.get("evidences")):
                    continue
                for sentence in sentences(text["value"]):
                    refs = pubmed_refs(sentence)
                    if not refs:        # "(By similarity)" and uncited sentences
                        continue
                    _function_sentence(protein, accession, sentence, refs)
        elif kind in NOTE_COMMENTS:
            _swissprot_note(protein, accession, kind, comment)

    for feature in entry.get("features", []):
        if feature["type"] in NOTE_FEATURES and _is_experimental(feature.get("evidences")):
            loc = feature.get("location", {})
            start, end = loc.get("start", {}).get("value"), loc.get("end", {}).get("value")
            where = f"{start}" if start == end else f"{start}-{end}"
            text = f"{feature['type']} {where}"
            if feature.get("alternativeSequence"):
                alt = feature["alternativeSequence"]
                text += f" {alt.get('originalSequence', '')}->{'/'.join(alt.get('alternativeSequences', []))}"
            if feature.get("description"):
                text += f": {feature['description']}"
            protein["notes"].append({"source": "Swiss-Prot", "accession": accession,
                                     "topic": feature["type"].lower(), "text": text,
                                     "refs": _evidence_refs(feature.get("evidences"))})


def _function_sentence(protein, accession, sentence, refs):
    plain = re.sub(r"\s*\((?:PubMed:\d+(?:, )?)+\)", "", sentence)
    positions_all, only_all = oxidation_in(plain)
    found_any = False
    for clause in clauses(plain):
        found = substrates_in(clause)
        if not found:
            continue
        found_any = True
        negative = bool(NEGATION.search(clause))
        positions, only = oxidation_in(clause)
        for category, form in found:
            protein["observations"].append({
                "source": "Swiss-Prot", "accession": accession, "substrate": category, "form": form,
                "result": "inactive" if negative else "active",
                "oxidation": [] if negative else positions, "only": only and not negative,
                "refs": refs, "text": sentence})
    if not found_any and positions_all and re.search(r"oxidi|product|cleav|hydroxylat", plain, re.I):
        protein["observations"].append({
            "source": "Swiss-Prot", "accession": accession, "substrate": None, "form": None,
            "result": "active", "oxidation": positions_all, "only": only_all, "refs": refs,
            "text": sentence})
    elif not found_any:
        protein["notes"].append({"source": "Swiss-Prot", "accession": accession, "topic": "function",
                                 "text": sentence, "refs": refs})


def _swissprot_note(protein, accession, kind, comment):
    texts = []
    if kind == "COFACTOR":
        for cofactor in comment.get("cofactors", []):
            if _is_experimental(cofactor.get("evidences")):
                texts.append((cofactor["name"], _evidence_refs(cofactor.get("evidences"))))
        note = comment.get("note") or {}
        texts += [(t["value"], _evidence_refs(t.get("evidences"))) for t in note.get("texts", [])
                  if _is_experimental(t.get("evidences"))]
    elif kind == "BIOPHYSICOCHEMICAL PROPERTIES":
        kinetics = comment.get("kineticParameters") or {}
        for km in kinetics.get("michaelisConstants", []):
            if _is_experimental(km.get("evidences")):
                texts.append((f"KM={km.get('constant')} {km.get('unit')} for {km.get('substrate')}",
                              _evidence_refs(km.get("evidences"))))
        for key in ("phDependence", "temperatureDependence", "redoxPotential", "absorption"):
            block = comment.get(key) or {}
            texts += [(f"{key}: {t['value']}", _evidence_refs(t.get("evidences")))
                      for t in block.get("texts", []) if _is_experimental(t.get("evidences"))]
        for t in (kinetics.get("note") or {}).get("texts", []):
            if _is_experimental(t.get("evidences")):
                texts.append((t["value"], _evidence_refs(t.get("evidences"))))
    else:
        texts = [(t["value"], _evidence_refs(t.get("evidences"))) for t in comment.get("texts", [])
                 if _is_experimental(t.get("evidences"))]
    for text, refs in texts:
        refs = refs or pubmed_refs(text)
        protein["notes"].append({"source": "Swiss-Prot", "accession": accession, "topic": kind.lower(),
                                 "text": text, "refs": refs})


# ---------------------------------------------------------------------------
# BRENDA
# ---------------------------------------------------------------------------

def find_brenda_file(data_dir):
    path = os.environ.get("BRENDA_JSON")
    if path:
        return path
    files = sorted(glob.glob(os.path.join(data_dir, "brenda_*.json*")))
    return files[-1] if files else None


def load_brenda(path):
    """(release, {ec: BRENDA entry}) of the LPMO EC numbers. The release file is ~700 MB of
    indented JSON with one line per EC key, so only the LPMO entries are parsed."""
    if path.endswith((".tar.gz", ".tgz")):
        archive = tarfile.open(path)
        member = next(m for m in archive.getmembers() if m.name.endswith(".json"))
        stream = io.TextIOWrapper(archive.extractfile(member), encoding="utf-8")
    else:
        stream = open(path, encoding="utf-8")
    release, entries, current, lines = None, {}, None, []
    key_line = re.compile(r'^\s{4}"([^"]+)": \{$')
    with stream:
        for line in stream:
            if release is None:
                m = re.match(r'^\s{2}"release": "([^"]+)"', line)
                if m:
                    release = m.group(1)
            if current is None:
                m = key_line.match(line.rstrip("\n"))
                if m and m.group(1) in LPMO_EC_NUMBERS:
                    current, lines = m.group(1), ["{"]
            elif line.rstrip("\n") in ("    },", "    }"):
                entries[current] = json.loads("".join(lines + ["}"]))
                current = None
                if len(entries) == len(LPMO_EC_NUMBERS):
                    break
            else:
                lines.append(line)
    return release, entries


def _brenda_ref(ec, ref_id, entry, references, bad_pmids=None):
    ref = entry.get("reference", {}).get(ref_id)
    if not ref:
        return f"BRENDA:{ec}:{ref_id}"
    pmid = ref.get("pmid")
    if pmid and not valid_pmid(pmid):
        # A malformed ID would become a dead PubMed link and would not merge with the same
        # paper cited by CAZy or Swiss-Prot; keep the claim under BRENDA's own reference.
        if bad_pmids is not None:
            bad_pmids.setdefault(str(pmid), set()).add(f"{ec}:{ref_id}")
        pmid = None
    key = f"PMID:{pmid}" if pmid else f"BRENDA:{ec}:{ref_id}"
    authors = ref.get("authors") or []
    references.setdefault(key, {"pmid": str(pmid) if pmid else None, "doi": None,
                                "title": ref.get("title"), "author": authors[0] if authors else None,
                                "n_authors": len(authors),
                                "year": str(ref.get("year") or "") or None, "journal": ref.get("journal")})
    return key


BRENDA_NOTE_FIELDS = ("ph_optimum", "ph_range", "temperature_optimum", "temperature_range", "km_value",
                      "turnover_number", "kcat_km_value", "protein_variants", "activating_compound",
                      "inhibitor", "general_information", "posttranslational_modification",
                      "temperature_stability")
NO_SUBSTRATE = re.compile(r"\b(?:no substrates?|no activity (?:on|with|towards)|not active (?:on|with|towards)|"
                          r"inactive (?:on|towards|against))\s*:?\s*([^.]+)", re.I)


GENERIC_NAME = re.compile(r"^(?:lpmos?|pmos?|aa|gh|cbm|cel|egl?)[ -]?\d*$|monooxygenase|polysaccharide|"
                          r"cellulase|glucanase|chitinase|protein|enzyme|\s", re.I)


def _specific_name(name):
    """A BRENDA synonym naming one enzyme (NcLPMO9C, PsLPMOA, Cel61A), not a family or activity."""
    return len(name) >= 4 and not GENERIC_NAME.search(name) and re.search(r"\d|[A-Z].*[A-Z]", name[1:])


def _names_protein(name, names):
    """Whether a catalogue protein's names hold a synonym, alone or after an organism prefix
    (LPMO9C in NcLPMO9C)."""
    n = re.escape(name.lower())
    return re.search(rf"(?<![a-z0-9])(?:[a-z]{{1,4}})?{n}(?![a-z0-9])", names) is not None


def collect_brenda(proteins, catalogue, entries, references):
    """Assign the BRENDA rows to catalogue proteins (see the module docstring). catalogue:
    {struct_id: {"accessions", "species", "refs", "names"}}. Returns the unmatched proteins."""
    by_accession = {acc: sid for sid, c in catalogue.items() for acc in c["accessions"]}
    by_species = defaultdict(list)
    for sid, c in catalogue.items():
        by_species[c["species"]].append(sid)
    unmatched = []
    bad_pmids = {}                          # malformed BRENDA PubMed ID -> {"<ec>:<ref id>"}

    def bref(ec, ref_id, entry):
        return _brenda_ref(ec, ref_id, entry, references, bad_pmids)

    for ec, entry in entries.items():
        synonyms = defaultdict(list)            # BRENDA protein id -> [(name, refs)]
        for syn in entry.get("synonyms", []):
            for pid in syn.get("proteins", []):
                synonyms[pid].append((syn["value"], set(syn.get("references", []))))

        def targets(pid, row_refs):
            """[(struct_id, how)] for a BRENDA protein in a row citing row_refs: by UniProt
            accession; else among the catalogue proteins of its species, the one named by a
            specific BRENDA synonym of the row's paper, or without such a synonym the only one
            citing the paper."""
            record = entry["protein"].get(pid, {})
            accessions = record.get("accessions") or []
            sids = {by_accession[a] for a in accessions if a in by_accession}
            if sids:
                return [(sid, "uniprot") for sid in sorted(sids)]
            pmids = {key for key in (bref(ec, r, entry) for r in row_refs)
                     if key.startswith("PMID:")}
            candidates = by_species.get(_species(record.get("organism")), [])
            names = {name.lower() for name, refs in synonyms.get(pid, [])
                     if _specific_name(name) and (refs & set(row_refs) or not refs)}
            named = [sid for sid in candidates if any(_names_protein(n, catalogue[sid]["names"]) for n in names)]
            citing = [sid for sid in candidates if pmids & catalogue[sid]["refs"]]
            suffix = " (BRENDA's accession differs)" if accessions else ""
            if len(named) == 1:
                how = "species + name" + (" + paper" if named[0] in citing else "")
                used = sorted(n for n in names if _names_protein(n, catalogue[named[0]]["names"]))
                return [(named[0], f"{how} ({', '.join(used)}){suffix}")]
            if not named and not names and len(citing) == 1:
                return [(citing[0], "species + paper" + suffix)]
            return []

        seen_proteins = {}                      # (BRENDA protein, struct_id) -> how
        for row in entry.get("substrates_products", []) + entry.get("natural_substrates_products", []):
            value, comment = row.get("value", ""), row.get("comment", "")
            substrate_text = value.split(" = ")[0].split(" + ")[0].strip()
            product_text = value.split(" = ")[1] if " = " in value else ""
            for pid in row.get("proteins", []):
                matches = targets(pid, row.get("references", []))
                record = entry["protein"].get(pid, {})
                if not matches:
                    if NOT_SUBSTRATE.match(substrate_text) and not comment:
                        continue
                    unmatched.append({"ec": ec, "protein": pid, "organism": record.get("organism"),
                                      "accessions": record.get("accessions") or [],
                                      "row": value, "reason": "accession not in the catalogue"
                                      if record.get("accessions") else "no single catalogue protein of the species"})
                    continue
                refs = [bref(ec, r, entry) for r in row.get("references", [])]
                for sid, how in matches:
                    # one match per protein pair, whichever row found it: rows of other papers
                    # give the same name match without "+ paper", which would repeat every note
                    if "+ paper" in how or (pid, sid) not in seen_proteins:
                        seen_proteins[(pid, sid)] = how
                    protein = proteins[sid]
                    if how != "uniprot":
                        _add(protein["ec"].setdefault(ec, {}), "BRENDA", refs)
                    obs = []
                    if not NOT_SUBSTRATE.match(substrate_text):
                        cats = substrates_in(substrate_text)
                        if not cats:
                            protein.setdefault("_unparsed", []).append(substrate_text)
                        positions, only = oxidation_in(product_text)
                        for category, form in cats[:1]:
                            obs.append({"substrate": category, "form": substrate_text, "result": "active",
                                        "oxidation": positions, "only": only})
                    for m in NO_SUBSTRATE.finditer(comment):
                        for category, form in substrates_in(m.group(1)):
                            obs.append({"substrate": category, "form": form, "result": "inactive",
                                        "oxidation": [], "only": False})
                    text = value + (f" ({comment})" if comment else "")
                    for o in obs:
                        protein["observations"].append({"source": "BRENDA", **o, "refs": refs, "text": text,
                                                        "ec": ec, "brenda_protein": pid, "match": how})
                    if not obs and comment:
                        protein["notes"].append({"source": "BRENDA", "topic": "substrates", "text": text,
                                                 "refs": refs, "brenda_protein": pid, "match": how})

        # EC membership and notes of the matched BRENDA proteins
        for (pid, sid), how in sorted(seen_proteins.items()):
            record = entry["protein"][pid]
            if how == "uniprot":
                refs = [bref(ec, r, entry) for r in record.get("references", [])]
                _add(proteins[sid]["ec"].setdefault(ec, {}), "BRENDA", refs)
            for field in BRENDA_NOTE_FIELDS:
                for item in entry.get(field, []):
                    if pid not in item.get("proteins", []):
                        continue
                    refs = [bref(ec, r, entry) for r in item.get("references", [])]
                    if how != "uniprot" and not set(refs) & catalogue[sid]["refs"]:
                        continue      # a pooled organism record: keep only the protein's papers
                    text = str(item.get("value", ""))
                    if item.get("comment"):
                        text += f" ({item['comment']})"
                    proteins[sid]["notes"].append({"source": "BRENDA", "topic": field, "text": text,
                                                   "refs": refs, "ec": ec, "brenda_protein": pid, "match": how})
    for pmid, where in sorted(bad_pmids.items()):
        print(f"  WARNING: BRENDA reference with a malformed PubMed ID {pmid} "
              f"({', '.join(sorted(where))}); cited as the BRENDA reference instead")
    unique = {(u["ec"], u["protein"], u["reason"]): u for u in unmatched}
    return list(unique.values())


# ---------------------------------------------------------------------------
# Summary and references
# ---------------------------------------------------------------------------

DP_WORDS = {"bi": 2, "tri": 3, "tetra": 4, "penta": 5, "hexa": 6, "hepta": 7}


def form_class(category, form):
    """The kind of substrate a form names, to tell a real disagreement (activity and none on
    the same kind) from activity on one form and none on another (beta- vs alpha-chitin)."""
    f = (form or "").lower()
    if category == "chitin":
        for kind in ("alpha", "beta", "colloidal"):
            if kind in f:
                return kind
    elif category == "cellulose":
        if re.search(r"phosphoric|pasc|amorphous|regenerated|swollen", f):
            return "amorphous"
        if re.search(r"avicel|crystalline|filter paper|cotton|bacterial|kraft|pulp", f):
            return "crystalline"
        if re.search(r"nanofibril|nanocellulose", f):
            return "nanocellulose"
    elif category in ("cello-oligosaccharides", "chito-oligosaccharides", "malto-oligosaccharides"):
        m = re.search(r"(bi|tri|tetra|penta|hexa|hepta)os|glc(\d)|\](\d)", f)
        if m:
            return f"DP{DP_WORDS.get(m.group(1)) or m.group(2) or m.group(3)}"
    return "any"


def summarise(protein):
    """Per substrate: the papers showing activity and no activity, per form too, and the
    papers for each oxidised position. The result is "conflict" when some paper finds no
    activity on a form (or on the substrate in general) where another finds activity; activity
    on one form and none on another is "active" (see the forms)."""
    summary = {}
    for obs in protein["observations"]:
        if not obs["substrate"]:
            continue
        s = summary.setdefault(obs["substrate"], {"result": None, "forms": {}, "active": [], "inactive": [],
                                                  "oxidation": {}, "sources": []})
        form = s["forms"].setdefault(obs.get("form") or obs["substrate"],
                                     {"kind": form_class(obs["substrate"], obs.get("form")),
                                      "active": [], "inactive": []})
        _add(form, obs["result"], obs["refs"])
        _add(s, obs["result"], obs["refs"])
        for position in obs["oxidation"]:
            _add(s["oxidation"], position, obs["refs"])
        if obs["source"] not in s["sources"]:
            s["sources"].append(obs["source"])
    for s in summary.values():
        active_kinds = {f["kind"] for f in s["forms"].values() if f["active"]}
        inactive_kinds = {f["kind"] for f in s["forms"].values() if f["inactive"]}
        clash = active_kinds & inactive_kinds or (active_kinds and "any" in inactive_kinds)
        s["result"] = ("conflict" if clash else "active" if s["active"] else "inactive")
        # activity reported on the substrate in general only, none on a named form of it
        s["check"] = bool(not clash and active_kinds == {"any"} and inactive_kinds - {"any"})
    return summary


def uniprot_names(entry):
    """Protein and gene names of a UniProtKB entry, lower case."""
    description = entry.get("proteinDescription", {})
    names = []
    for block in [description.get("recommendedName") or {}] + description.get("alternativeNames", []) \
            + description.get("submissionNames", []):
        names.append((block.get("fullName") or {}).get("value", ""))
        names += [s.get("value", "") for s in block.get("shortNames", [])]
    for gene in entry.get("genes", []):
        names.append((gene.get("geneName") or {}).get("value", ""))
        for key in ("synonyms", "orderedLocusNames", "orfNames"):
            names += [g.get("value", "") for g in gene.get(key, [])]
    return " ".join(n for n in names if n).lower()


def rcsb_citations(pdb_codes):
    """{PDB code: 'PMID:...' or 'DOI:...'} of the entries' primary citations."""
    codes = sorted({c.upper() for c in pdb_codes if re.fullmatch(r"[0-9][A-Za-z0-9]{3}", c)})
    citations = {}
    for i in range(0, len(codes), 100):
        chunk = codes[i:i + 100]
        query = ("{ entries(entry_ids: %s) { rcsb_id rcsb_primary_citation "
                 "{ pdbx_database_id_PubMed pdbx_database_id_DOI } } }" % json.dumps(chunk))
        try:
            r = requests.post("https://data.rcsb.org/graphql", json={"query": query}, timeout=60)
            r.raise_for_status()
        except requests.RequestException as e:
            print(f"  Warning: RCSB citation lookup failed: {e}")
            continue
        for entry in (r.json().get("data") or {}).get("entries") or []:
            cite = entry.get("rcsb_primary_citation") or {}
            if cite.get("pdbx_database_id_PubMed"):
                citations[entry["rcsb_id"]] = f"PMID:{cite['pdbx_database_id_PubMed']}"
            elif cite.get("pdbx_database_id_DOI"):
                citations[entry["rcsb_id"]] = f"DOI:{cite['pdbx_database_id_DOI'].lower()}"
    return citations


CURATED_CSV = os.path.join(CURATION_DIR, "activities.csv")
CATEGORIES = {category for _, category in SUBSTRATE_PATTERNS}


def collect_curated(proteins, path=CURATED_CSV):
    """The hand-curated claims (source "Literature"); rows for proteins not in the selection or
    with values the collector does not know are reported and skipped."""
    if not os.path.exists(path):
        return 0
    curated = pd.read_csv(path, dtype=str, keep_default_na=False, comment="#")
    used = 0
    for row in curated.itertuples():
        oxidation = sorted({p for p in re.split(r"[/+, ]+", row.oxidation.upper()) if p})
        problems = [what for what, bad in (
            ("protein not in the selection", row.struct_id not in proteins),
            (f"substrate '{row.substrate}'", row.substrate not in CATEGORIES),
            (f"result '{row.result}'", row.result not in ("active", "inactive")),
            (f"oxidation '{row.oxidation}'", not set(oxidation) <= {"C1", "C4", "C6"}),
            ("no PubMed ID or DOI", not (row.pmid or row.doi)),
            (f"malformed PubMed ID '{row.pmid}'", bool(row.pmid) and not valid_pmid(row.pmid))) if bad]
        if problems:
            print(f"  WARNING: curated claim {row.struct_id} {row.substrate} skipped: {'; '.join(problems)}")
            continue
        ref = f"PMID:{row.pmid}" if row.pmid else f"DOI:{row.doi.lower()}"
        proteins[row.struct_id]["observations"].append({
            "source": "Literature", "substrate": row.substrate, "form": row.form or None,
            "result": row.result, "oxidation": oxidation if row.result == "active" else [],
            "only": False, "refs": [ref], "text": row.evidence})
        used += 1
    return used


def _title_key(title):
    return re.sub(r"[^a-z0-9]+", " ", (title or "").lower()).strip()


def merge_references_by_title(proteins, references):
    """Fold a "BRENDA:<ec>:<id>" reference onto the "PMID:..." one with the same title, so a
    paper BRENDA holds without a usable PubMed ID (none, or a malformed one) is not cited a
    second time beside the CAZy or Swiss-Prot entry for it. Titles are compared without case
    or punctuation, and only a BRENDA key is ever replaced."""
    by_title = {}
    for key, ref in references.items():
        if key.startswith("PMID:") and _title_key(ref.get("title")):
            by_title.setdefault(_title_key(ref["title"]), key)
    rename = {key: by_title[_title_key(ref.get("title"))]
              for key, ref in references.items()
              if key.startswith("BRENDA:") and _title_key(ref.get("title")) in by_title}
    if not rename:
        return {}

    def swap(node):
        if isinstance(node, dict):
            return {rename.get(k, k): swap(v) for k, v in node.items()}
        if isinstance(node, list):
            items = [swap(item) for item in node]
            if any(isinstance(item, str) and item in rename for item in node):
                items = list(dict.fromkeys(items))   # only the duplicate this merge created
            return items
        return rename.get(node, node) if isinstance(node, str) else node

    for sid, protein in proteins.items():
        proteins[sid] = swap(protein)
    for key in rename:
        references.pop(key, None)
    for old, new in sorted(rename.items()):
        print(f"  {old} merged into {new} (same title)")
    return rename


def fill_reference_metadata(references, refs):
    """Title, first author, year and journal of PubMed IDs not seen in the sources (Europe PMC)."""
    missing = sorted({r[5:] for r in refs if r.startswith("PMID:") and not references.get(r, {}).get("title")})
    for i in range(0, len(missing), 100):
        chunk = missing[i:i + 100]
        query = " OR ".join(f"EXT_ID:{p}" for p in chunk)
        try:
            r = requests.get(EUROPEPMC_URL, params={"query": f"({query}) AND SRC:MED", "format": "json",
                                                    "pageSize": 1000, "resultType": "lite"}, timeout=60)
            r.raise_for_status()
        except requests.RequestException as e:
            print(f"  Warning: Europe PMC lookup failed: {e}")
            continue
        for hit in r.json().get("resultList", {}).get("result", []):
            references[f"PMID:{hit['pmid']}"] = {
                "pmid": hit["pmid"], "doi": hit.get("doi"), "title": hit.get("title"),
                "author": (hit.get("authorString") or "").split(",")[0] or None,
                "n_authors": len([a for a in (hit.get("authorString") or "").split(",") if a.strip()]),
                "year": hit.get("pubYear"), "journal": hit.get("journalTitle")}
        time.sleep(0.2)
    for ref in refs:
        if ref.startswith("DOI:"):
            references.setdefault(ref, {"pmid": None, "doi": ref[4:], "title": None, "author": None,
                                        "year": None, "journal": None})


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def collect_activities(selection, cazy, identity, uniprot_references, taxonomy, brenda_path):
    """The activities.json content (see the module docstring)."""
    kept = selection[selection.source != "NONE"]
    proteins = {sid: {"ec": {}, "observations": [], "substrates": {}, "notes": []}
                for sid in kept.selected_pdb}
    references = {}

    # CAZy
    cazy_rows = cazy_rows_of(selection, cazy)
    for sid, entries in cazy_rows.items():
        collect_cazy(proteins[sid], entries)

    # Swiss-Prot: the reviewed entries of each protein (2f); the others only for their names
    active = {sid: [acc for acc, info in identity.get(sid, {}).get("accessions", {}).items()
                    if info.get("status") in ("reviewed", "unreviewed")] for sid in proteins}
    reviewed = {sid: [acc for acc in accs if identity[sid]["accessions"][acc]["status"] == "reviewed"]
                for sid, accs in active.items()}
    entries = fetch_uniprot_entries({acc for accs in active.values() for acc in accs})
    for sid, accs in reviewed.items():
        for acc in accs:
            if acc in entries:
                collect_swissprot(proteins[sid], acc, entries[acc], references)

    # hand-curated claims (curation/activities.csv)
    curated = collect_curated(proteins)
    if curated:
        print(f"  {curated} curated claims from {os.path.relpath(CURATED_CSV, PIPELINE_DIR)}")

    # what the catalogue knows of each protein, for matching BRENDA
    pdb_codes = {sid: {c.upper() for c in _ids(r.all_pdb_codes) + ([sid] if r.source == "PDB" else [])}
                 for sid, r in zip(kept.selected_pdb, kept.itertuples())}
    structure_papers = rcsb_citations({c for codes in pdb_codes.values() for c in codes})
    catalogue = {}
    for r in kept.itertuples():
        refs = {ref for entry in cazy_rows.get(r.selected_pdb, [])
                for ref in map(_ref_from_url, _ids(entry.references)) if ref}
        refs |= {structure_papers[c] for c in pdb_codes[r.selected_pdb] if c in structure_papers}
        refs |= {f"PMID:{p['pmid']}" for p in uniprot_references.get(r.selected_pdb, []) if p.get("pmid")}
        refs |= {ref for obs in proteins[r.selected_pdb]["observations"] for ref in obs["refs"]
                 if ref.startswith("PMID:")}
        species = ((taxonomy.get(r.organism) or {}).get("species") or {}).get("name") or r.organism
        names = " ".join([r.protein_name] + [e.protein_name for e in cazy_rows.get(r.selected_pdb, [])]
                         + [uniprot_names(entries[acc]) for acc in active[r.selected_pdb] if acc in entries]).lower()
        catalogue[r.selected_pdb] = {
            "accessions": set(identity.get(r.selected_pdb, {}).get("accessions", {})) | set(_ids(r.uniprot_ids)),
            "species": _species(species), "refs": refs, "names": names}

    release, unmatched = None, []
    if brenda_path and os.path.exists(brenda_path):
        release, brenda = load_brenda(brenda_path)
        unmatched = collect_brenda(proteins, catalogue, brenda, references)
    else:
        print("  WARNING: no BRENDA JSON (data/brenda_*.json.tar.gz or BRENDA_JSON) - BRENDA skipped")

    merge_references_by_title(proteins, references)

    unparsed = Counter()
    for protein in proteins.values():
        unparsed.update(protein.pop("_unparsed", []))
        unique = {}
        for obs in protein["observations"]:
            key = (obs["source"], obs["substrate"], obs.get("form"), obs["result"], tuple(obs["oxidation"]),
                   tuple(obs["refs"]), obs["text"], obs.get("accession"), obs.get("brenda_protein"))
            unique.setdefault(key, obs)
        protein["observations"] = list(unique.values())
        protein["substrates"] = summarise(protein)
    all_refs = {ref for p in proteins.values() for part in ("observations", "notes")
                for item in p[part] for ref in item["refs"]}
    all_refs |= {ref for p in proteins.values() for sources in p["ec"].values()
                 for refs in sources.values() for ref in refs}
    fill_reference_metadata(references, all_refs)

    return {"sources": {"checked": datetime.now().strftime("%Y-%m-%d"),
                        "brenda": {"file": os.path.basename(brenda_path) if brenda_path else None,
                                   "release": release}},
            "proteins": proteins,
            "references": {ref: references[ref] for ref in sorted(all_refs) if ref in references},
            "unmatched": {"brenda": unmatched},
            "unparsed": sorted(unparsed)}


CSV_COLUMNS = ("struct_id", "family", "protein", "organism", "substrate", "form", "form_kind", "result",
               "oxidation", "only", "source", "match", "papers", "refs", "substrate_result", "flag", "text")


def review_rows(result, selection):
    """One row per claim for reviewing: the observation with its papers, and the summary result
    of its substrate. EC numbers come first per protein, as rows with substrate "(EC number)"."""
    info = {r.selected_pdb: r for r in selection.itertuples()}

    def paper(ref):
        """'Karkehabadi S. 2008' -> 'Karkehabadi 2008' (the sources write the first author
        differently: 'Simmons T.J.', 'Simmons, T.J.', 'Simmons TJ')."""
        meta = result["references"].get(ref, {})
        words = (meta.get("author") or "").replace(",", " ").split()
        while words and (words[-1].endswith(".") or (len(words[-1]) <= 3 and words[-1].isupper())):
            words.pop()
        return f"{' '.join(words)} {meta.get('year') or ''}".strip() or ref

    rows = []
    for sid, protein in sorted(result["proteins"].items()):
        row = info.get(sid)
        base = {"struct_id": sid, "family": getattr(row, "family", ""),
                "protein": getattr(row, "protein_name", ""), "organism": getattr(row, "organism", "")}
        for ec, sources in sorted(protein["ec"].items()):
            for source, refs in sorted(sources.items()):
                rows.append({**base, "substrate": "(EC number)", "form": ec, "result": "",
                             "source": source, "papers": "; ".join(paper(r) for r in refs),
                             "refs": "; ".join(refs),
                             "flag": "EC numbers differ between sources"
                                     if len({frozenset(s) for s in
                                             [{e for e, ss in protein["ec"].items() if src in ss}
                                              for src in {k for ss in protein["ec"].values() for k in ss}]}) > 1
                                     else ""})
        for obs in sorted(protein["observations"], key=lambda o: (o["substrate"] or "", o["source"])):
            summary = protein["substrates"].get(obs["substrate"], {})
            flag = ("activity and no activity on the same substrate" if summary.get("result") == "conflict"
                    else "activity on the substrate, none on this form" if summary.get("check") else "")
            rows.append({**base, "substrate": obs["substrate"] or "(any substrate)",
                         "form": obs.get("form") or "", "form_kind": form_class(obs["substrate"], obs.get("form")),
                         "result": obs["result"], "oxidation": "/".join(obs["oxidation"]),
                         "only": "only" if obs["only"] else "", "source": obs["source"],
                         "match": obs.get("match") or obs.get("accession") or "",
                         "papers": "; ".join(paper(r) for r in obs["refs"]), "refs": "; ".join(obs["refs"]),
                         "substrate_result": summary.get("result", ""), "flag": flag,
                         "text": " ".join(obs["text"].split())})
    return rows


def write_review_csv(result, selection, path):
    rows = review_rows(result, selection)
    pd.DataFrame(rows, columns=CSV_COLUMNS).fillna("").to_csv(path, index=False)
    return len(rows)


def report(result):
    """Print coverage, disagreements between the sources and what could not be assigned."""
    proteins = result["proteins"]
    by_source = Counter()
    for p in proteins.values():
        for source in {o["source"] for o in p["observations"]}:
            by_source[source] += 1
    with_data = sum(bool(p["substrates"]) for p in proteins.values())
    print(f"  {with_data} of {len(proteins)} proteins with substrate data; per source: "
          + ", ".join(f"{s} {n}" for s, n in by_source.most_common()))
    oxidation = sum(any(s["oxidation"] for s in p["substrates"].values())
                    or any(o["oxidation"] for o in p["observations"]) for p in proteins.values())
    print(f"  {oxidation} proteins with a C1/C4 statement from a paper (Swiss-Prot text or BRENDA product)")
    categories = Counter(sub for p in proteins.values() for sub, s in p["substrates"].items()
                         if s["result"] != "inactive")
    print("  substrates with activity: " + ", ".join(f"{c} {n}" for c, n in categories.most_common()))
    inactive = Counter(sub for p in proteins.values() for sub, s in p["substrates"].items() if s["inactive"])
    print("  substrates reported without activity: " + ", ".join(f"{c} {n}" for c, n in inactive.most_common()))

    ec_disagree = []
    for sid, p in proteins.items():
        by = defaultdict(set)
        for ec, sources in p["ec"].items():
            for source in sources:
                by[source].add(ec)
        if len(by) > 1 and len({frozenset(v) for v in by.values()}) > 1:
            ec_disagree.append(f"{sid} " + " | ".join(f"{s} {'; '.join(sorted(v))}" for s, v in sorted(by.items())))
    print(f"  EC numbers differing between sources ({len(ec_disagree)}):")
    for line in ec_disagree:
        print(f"    {line}")
    conflicts = [f"{sid} {sub}" for sid, p in proteins.items() for sub, s in p["substrates"].items()
                 if s["result"] == "conflict"]
    print(f"  activity and no activity on the same substrate ({len(conflicts)}): {', '.join(conflicts)}")
    checks = [f"{sid} {sub}" for sid, p in proteins.items() for sub, s in p["substrates"].items() if s["check"]]
    print(f"  activity on the substrate, none on a named form of it ({len(checks)}): {', '.join(checks)}")
    names = Counter(o["match"] for p in proteins.values() for o in p["observations"]
                    if o["source"] == "BRENDA" and o["match"] != "uniprot")
    print("  BRENDA proteins matched without a shared accession: "
          + "; ".join(f"{m} {n}" for m, n in sorted(names.items())))
    unmatched = result["unmatched"]["brenda"]
    reasons = Counter(u["reason"] for u in unmatched)
    print(f"  BRENDA rows not assigned: " + ", ".join(f"{r} {n}" for r, n in reasons.items()))
    if result["unparsed"]:
        print(f"  substrate names without a category: {'; '.join(result['unparsed'])}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--data-dir", default=DEFAULT_DATA_DIR)
    parser.add_argument("--brenda", help="BRENDA JSON release (.json or .json.tar.gz)")
    parser.add_argument("--no-csv", action="store_true", help="skip activities_review.csv")
    args = parser.parse_args()
    structures = os.path.join(args.data_dir, "02_data")
    selection = pd.read_csv(os.path.join(structures, "pdb_selection.csv"), dtype=str, keep_default_na=False)
    cazy = pd.read_csv(os.path.join(args.data_dir, "01_proteins", "cazy_lpmos.csv"), dtype=str, keep_default_na=False)
    with open(os.path.join(structures, "uniprot_ids.json")) as f:
        identity = json.load(f)
    with open(os.path.join(structures, "uniprot_references.json")) as f:
        uniprot_references = json.load(f)
    with open(os.path.join(structures, "taxonomy.json")) as f:
        taxonomy = json.load(f)
    brenda = args.brenda or find_brenda_file(args.data_dir)
    print(f"2j: substrates and oxidised positions from CAZy, Swiss-Prot and BRENDA ({brenda})...")
    result = collect_activities(selection, cazy, identity, uniprot_references, taxonomy, brenda)
    report(result)
    out = os.path.join(structures, "activities.json")
    with open(out, "w") as f:
        json.dump(result, f, indent=1)
    print(f"Saved {out}")
    if not args.no_csv:
        csv_path = os.path.join(structures, "activities_review.csv")
        print(f"Saved {write_review_csv(result, selection, csv_path)} rows to {csv_path}")


if __name__ == "__main__":
    main()
