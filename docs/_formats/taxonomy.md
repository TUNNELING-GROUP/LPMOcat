JSON object, one key per organism name as it stands in `pdb_selection.csv`; the value is null
for a name NCBI Taxonomy does not know.

| Field | Meaning |
| --- | --- |
| `taxid`, `name`, `rank` | The NCBI taxon the name resolved to |
| `matched` | `name` (the full name), `species` (its first two words) or `genus` |
| `lineage` | List of `[rank, name]` from the root down to the taxon itself |
| `species` | `{taxid, name}` of the taxon's species, or null above species level |
