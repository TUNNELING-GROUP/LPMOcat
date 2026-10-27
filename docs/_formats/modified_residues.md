JSON object, one key per prepared structure: the modified amino acids of the deposited chain
that remain in the prepared structure, each replaced there by its parent (3f).

| Field | Meaning |
| --- | --- |
| `residue` | Residue number in the structure's own numbering |
| `name` | The modified residue's PDB code, e.g. `HIC` (N-methylated His) |
| `parent` | The standard amino acid it became, e.g. `HIS` |

An empty list means none: every model, and every crystal structure deposited with standard
residues only. The residues are read from the downloaded entry (first model, the selected
chain, first alternate conformer), not from the 3a output, whose rebuilt file already has the
parent residues.
