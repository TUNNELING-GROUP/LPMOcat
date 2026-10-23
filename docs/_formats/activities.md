JSON object with five keys.

| Key | Content |
| --- | --- |
| `sources` | `checked` (date) and `brenda`: `{file, release}` of the BRENDA release read |
| `proteins` | One key per structure, see below |
| `references` | `{ref: {title, author, n_authors, year, journal, doi, pmid}}` for every paper cited |
| `unmatched` | `brenda`: BRENDA proteins and rows that could not be assigned to a catalogue protein, each with a `reason` |
| `unparsed` | Substrate texts no substrate category matched |

A `ref` is `PMID:<id>`, `DOI:<doi>` or, for a BRENDA paper without a usable PubMed ID,
`BRENDA:<ec>:<reference number>`.

Each protein has:

| Field | Meaning |
| --- | --- |
| `ec` | `{EC number: {source: [refs]}}` |
| `observations` | One item per claim: `source` (`CAZy`, `Swiss-Prot`, `BRENDA`, `Literature`), `substrate` (category, or null for a claim about the oxidised position alone), `form` (the substrate as written, e.g. `beta-chitin`), `result` (`active` / `inactive`), `oxidation` (list of `C1`, `C4`, `C6`), `only` (true when the text excludes the other position), `refs`, `text` (the sentence or reaction), and `accession` or `brenda_protein` with `match` (how a BRENDA protein was assigned) |
| `substrates` | Summary per substrate category: `result` (`active`, `inactive` or `conflict`), `forms` (`{form: {kind, active: [refs], inactive: [refs]}}`), `active`, `inactive` (refs), `oxidation` (`{position: [refs]}`), `sources`, and `check` (true when activity is reported only for the substrate in general while a named form is reported inactive) |
| `notes` | `{source, topic, text, refs}`: kinetics, pH and temperature optima, regulation, modifications, mutants |
