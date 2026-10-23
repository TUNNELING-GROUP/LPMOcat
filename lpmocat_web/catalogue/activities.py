"""Substrates, oxidised positions and the papers reporting them (Protein.activities, pipeline
step 2j), as the compare page's protein table shows them: the EC numbers, the substrates with
reported activity and the oxidised positions, each followed by the numbers of the papers behind
it in the row's References column.

The pipeline collects the claims from CAZy's EC numbers, the reviewed UniProt (Swiss-Prot)
entries and BRENDA without preferring a source; this module only arranges them and cites the
papers. Substrates tested without activity are not listed in the table -- they stay in the
hover text and in the CSV export.
"""
import json

# Substrate rows in this order (the rest follow alphabetically); the polymer before its
# oligosaccharides, cellulose and chitin first
SUBSTRATE_ORDER = ("cellulose", "cello-oligosaccharides", "cellulose derivatives", "chitin",
                   "chito-oligosaccharides", "chitosan", "xyloglucan", "xyloglucan oligosaccharides",
                   "mixed-linkage glucan", "beta-1,3-glucan", "glucomannan", "mannan", "xylan",
                   "xylo-oligosaccharides", "starch", "malto-oligosaccharides", "pectin", "galactan",
                   "peptidoglycan", "alginate", "biomass")
SUBSTRATE_LABELS = {"beta-1,3-glucan": "β-1,3-glucan", "biomass": "Plant biomass"}

SOURCE_LABELS = {"CAZy": "CAZy (EC number)", "Swiss-Prot": "Swiss-Prot", "BRENDA": "BRENDA",
                 "Literature": "the literature (curated)"}


def activity_data(protein):
    """The protein's stored activity record (step 2j), or None."""
    try:
        data = json.loads(getattr(protein, "activities", None) or "")
    except (ValueError, TypeError):
        return None
    return data or None


def substrate_label(name):
    return SUBSTRATE_LABELS.get(name) or name[:1].upper() + name[1:]


def _reference_url(meta):
    if meta.get("pmid"):
        return f"https://pubmed.ncbi.nlm.nih.gov/{meta['pmid']}/"
    if meta.get("doi"):
        return f"https://doi.org/{meta['doi']}"
    return ""


def _surname(author):
    """'Simmons T.J.' / 'Simmons, T.J.' / 'Simmons TJ' -> 'Simmons' (the three sources write the
    first author differently), keeping compound names ('Lo Leggio', 'van den Berg')."""
    words = (author or "").replace(",", " ").split()
    while words and (words[-1].endswith(".") or (len(words[-1]) <= 3 and words[-1].isupper())):
        words.pop()
    return " ".join(words) or (author or "").strip()


def _reference_label(meta, ref):
    """'Frandsen et al., 2017' from the stored metadata."""
    author = _surname(meta.get("author"))
    year = meta.get("year") or ""
    if author and year:
        return f"{author} et al., {year}"
    return author or str(year) or ref.split(":")[-1]


def activity_references(data):
    """The papers of a protein's activity record as entries for the table's References column:
    [{"key", "pmid", "doi", "author", "year", "title", "journal"}]; the papers already listed
    there (CAZy, RCSB, UniProt) are deduplicated by _format_references."""
    if not data:
        return []
    used = {ref for item in data["observations"] + data["notes"] for ref in item["refs"]}
    used |= {ref for sources in data["ec"].values() for refs in sources.values() for ref in refs}
    entries = []
    for ref in sorted(used):
        meta = (data.get("references") or {}).get(ref, {})
        entries.append({"key": ref, "pmid": meta.get("pmid"), "doi": meta.get("doi"),
                        "author": _surname(meta.get("author")), "n_authors": meta.get("n_authors"),
                        "year": meta.get("year"),
                        "title": meta.get("title"), "journal": meta.get("journal")})
    return entries


def _oxidation(data):
    """{position: [refs]} over the substrates and the statements without a substrate ("yielding
    C4 oxidation products"), and whether a paper says the enzyme oxidises only there."""
    positions, only = {}, False
    for summary in data["substrates"].values():
        for position, refs in summary["oxidation"].items():
            known = positions.setdefault(position, [])
            known += [r for r in refs if r not in known]
    for obs in data["observations"]:
        if obs["substrate"] is None:
            for position in obs["oxidation"]:
                known = positions.setdefault(position, [])
                known += [r for r in obs["refs"] if r not in known]
            only = only or (obs["only"] and len(obs["oxidation"]) == 1)
    return positions, only


def _numbers(refs, number_of):
    """All the papers of a claim as [{"n", "label"}], by the numbers of the row's References
    column; a paper the column does not list (no citable ID) is left out. The table shows them
    all as superscripts (user, 2026-09-23), separated by commas so 1,2 does not read as 12."""
    found = {}
    for ref in refs:
        entry = number_of(ref)
        if entry:
            found[entry["n"]] = entry
    return [found[n] for n in sorted(found)]


def _substrate_title(substrate, summary):
    """Hover text: which forms were tested, and from which sources."""
    forms = sorted(summary["forms"].items(), key=lambda item: item[0] != substrate)
    active, inactive = [], []
    for form, info in forms:
        label = form if form.lower() != substrate.lower() else substrate
        for result, names in (("active", active), ("inactive", inactive)):
            if info[result] and not any(label.lower() == seen.lower() for seen in names):
                names.append(label)
    parts = []
    if active:
        parts.append("activity on " + ", ".join(active))
    if inactive:
        parts.append("no activity reported on " + ", ".join(inactive))
    sources = ", ".join(SOURCE_LABELS.get(s, s) for s in summary["sources"])
    return "; ".join(parts) + (f" — from {sources}" if sources else "")


NO_PAPER = "from the EC numbers or the protein name — no paper in the sources reports it"


def _substrate_cells(data, number_of, substrate_names):
    """The protein's substrates in the table's order, each with the papers reporting activity
    on it (none for a substrate that only follows from an EC number or the protein name)."""
    summaries = {substrate_label(name): (name, summary)
                 for name, summary in (data or {}).get("substrates", {}).items()}
    order = {substrate_label(name): i for i, name in enumerate(SUBSTRATE_ORDER)}
    cells = []
    for label in sorted(set(substrate_names), key=lambda n: (order.get(n, len(order)), n)):
        name, summary = summaries.get(label, (None, None))
        if not summary or not summary["active"]:
            cells.append({"name": label, "refs": [], "conflict": False,
                          "title": NO_PAPER})
            continue
        cells.append({"name": label, "refs": _numbers(summary["active"], number_of),
                      "conflict": summary["result"] == "conflict",
                      "title": _substrate_title(name, summary)})
    return cells


# sources whose claims are a paper's own statement (not an EC number)
PAPER_SOURCES = ("Swiss-Prot", "BRENDA", "Literature")


def substrates_from_papers(protein):
    """The substrates (as the catalogue names them) a paper reports activity on -- a sentence of a
    Swiss-Prot entry, a BRENDA row or a curated claim; not a substrate that only follows from an
    EC number (CAZy's included), the protein name or the family."""
    data = activity_data(protein)
    if not data:
        return set()
    return {substrate_label(obs["substrate"]) for obs in data["observations"]
            if obs["substrate"] and obs["result"] == "active" and obs["source"] in PAPER_SOURCES}


def regioselectivity_from_papers(protein):
    """Whether the protein's regioselectivity is the oxidised positions papers report (step 2j)
    rather than what its EC numbers define (populate_db.reported_activity makes the same
    choice)."""
    data = activity_data(protein)
    if not data:
        return False
    positions, _ = _oxidation(data)
    return bool({"C1", "C4"} & set(positions))


def table_cells(data, number_of, substrate_names=()):
    """The activity part of a protein's row in the compare page's table:

    {"ec": [{"text", "refs": [{"n", "label"}], "title"}], "ec_differs": bool,
     "substrates": [{"name", "refs": [...], "title", "conflict"}],
     "regio": {"refs": [...], "title"} or None}   # the papers behind the regioselectivity

    The substrates are the protein's own (Protein.substrate, what the browse filters and the
    overview charts use): those a paper reports activity on, filled up from the EC numbers and
    the protein name. Only the reported ones carry papers; the others say so in their hover
    text. Substrates tested without activity are not listed -- they appear in the hover text of
    their substrate and in the CSV export.
    """
    substrates = _substrate_cells(data, number_of, substrate_names)
    if not data:
        return {"ec": [], "ec_differs": False, "substrates": substrates, "regio": None}

    by_source = {}
    for number, sources in data["ec"].items():
        for source, refs in sources.items():
            by_source.setdefault(source, {})[number] = refs
    ec = []
    for number in sorted({n for numbers in by_source.values() for n in numbers}):
        refs = [ref for numbers in by_source.values() for ref in numbers.get(number, [])]
        sources = [SOURCE_LABELS.get(s, s) for s, numbers in sorted(by_source.items()) if number in numbers]
        ec.append({"text": number, "refs": _numbers(refs, number_of),
                   "title": f"{number} in " + ", ".join(sources)})
    differs = len({frozenset(numbers) for numbers in by_source.values()}) > 1

    # the regioselectivity the catalogue shows comes from these papers when they give a
    # position, and from the EC numbers otherwise (populate_db.reported_activity)
    positions, _ = _oxidation(data)
    carbons = [position for position in ("C1", "C4") if position in positions]
    regio = None
    if carbons:
        refs = [ref for position in carbons for ref in positions[position]]
        regio = {"refs": _numbers(refs, number_of),
                 "title": "the oxidised positions these papers report; without them the EC "
                          "numbers decide"}
    return {"ec": ec, "ec_differs": differs, "substrates": substrates, "regio": regio}


def export_rows(proteins):
    """The card as CSV rows (protein, substrate, form, result, oxidised positions, sources,
    papers, and the sentence or reaction each claim comes from)."""
    rows = []
    for protein in proteins:
        data = activity_data(protein)
        if not data:
            continue
        code = protein.protein_code or protein.selected_structure
        references = data.get("references") or {}

        def papers(refs):
            return "; ".join(_reference_label(references.get(r, {}), r) for r in refs)

        for number, sources in sorted(data["ec"].items()):
            for source, srefs in sorted(sources.items()):
                rows.append([code, "(EC number)", number, "", "", source, papers(srefs), ""])
        for obs in data["observations"]:
            rows.append([code, obs["substrate"] or "(any substrate)", obs.get("form") or "",
                         obs["result"], "/".join(obs["oxidation"]), obs["source"],
                         papers(obs["refs"]), " ".join(obs["text"].split())])
    return rows


EXPORT_HEADER = ["protein", "substrate", "form", "result", "oxidised_positions", "source",
                 "papers", "reported_as"]
