JSON object, one key per model sent to SignalP 6.

| Field | Meaning |
| --- | --- |
| `prediction` | SignalP's class: `SP`, `LIPO`, `TAT`, …, or `OTHER` (no signal peptide) |
| `sp_end` | Last position of the signal peptide (1-based position in the model's sequence), or null |
| `cleavage_probability` | SignalP's probability of that cleavage site |
| `organism_group` | `eukarya` or `other`, the SignalP organism option used |
| `first_residue_after` | Residue number following the signal peptide, the candidate N-His |
