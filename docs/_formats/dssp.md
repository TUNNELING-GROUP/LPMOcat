JSON object, one key per structure.

| Field | Meaning |
| --- | --- |
| `residues` | Residue numbers of the standard amino acids, in file order |
| `ss8` | String, one DSSP state per residue (`H G I E B T S P`), `-` for none or for a residue DSSP skipped |
| `sheet` | String, one sheet label per residue (`A`, `B`, …), `-` outside sheets |
| `bridges` | List of `[residue, partner, ladder label]`, each pair once, by residue number. A lowercase label is a parallel ladder, uppercase antiparallel |
