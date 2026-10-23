CSV with the columns of `cazy_lpmos.csv` (`in_cazy_characterized` and `in_cazy_structure` are
`no`; `references` are PubMed links) and five more:

| Column | Meaning |
| --- | --- |
| `list_source` | `Swiss-Prot`, `BRENDA`, `curated` or `local model` |
| `evidence` | Why the entry is in: the EC numbers and papers, or the curated reason, followed by how its family was decided |
| `variant_of` | Accession of the CAZy-listed protein this entry is a strain variant of, else empty |
| `struct_id` | 1e only: the identifier of the own model (`AF3-…`) |
| `local_file` | 1e only: path of the model file, relative to `lpmocat_pipeline/` |
