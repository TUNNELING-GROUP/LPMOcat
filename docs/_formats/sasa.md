JSON object, one key per structure; five parallel lists with one item per residue, in file
order.

| Field | Meaning |
| --- | --- |
| `residues` | [Residue numbers](../conventions.md#numbering) |
| `chains` | Chain identifier of each residue |
| `names` | Three-letter residue names |
| `sasa` | Solvent-accessible surface area in Å² |
| `rsa` | Relative SASA in %, two decimals; null for a residue type without a maximum |
