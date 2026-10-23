JSON object, one key per structure; parallel lists with one item per residue (CA atom), in
file order.

| Field | Present | Meaning |
| --- | --- | --- |
| `residues` | always | [Residue numbers](../conventions.md#numbering) |
| `bfactor` | experimental structures | CA B-factor in Å², two decimals; null for residues rebuilt in step 3a |
| `plddt` | models | pLDDT (0–100) |
| `gnm` | when the model could be built | Mean-square fluctuation from the Gaussian network model, arbitrary units, six decimals |
| `cabsflex` | when a current CABS-flex run exists | RMSF in Å, four decimals |
