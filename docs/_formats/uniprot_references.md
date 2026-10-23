JSON object: structure identifier → list of papers, one per PubMed ID or DOI. Structures whose
entries cite no qualifying paper are absent.

| Field | Meaning |
| --- | --- |
| `pmid`, `doi` | Identifiers; at least one is present |
| `author`, `n_authors` | First author's surname and the number of authors |
| `year`, `title`, `journal` | Citation |
| `scopes` | What UniProt cites the paper for, e.g. `X-RAY CRYSTALLOGRAPHY (1.80 ANGSTROMS) OF 28-197` |
| `accessions` | The UniProt entries citing it |
