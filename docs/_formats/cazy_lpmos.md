CSV, one row per protein, sorted by family, organism and protein name. Multi-valued cells are
separated by `; ` (PDB codes by `, `).

| Column | Meaning |
| --- | --- |
| `family` | CAZy family, `AA9` … `AA17` |
| `kingdom` | CAZy's kingdom heading the protein is listed under (`Bacteria`, `Eukaryota`, …) |
| `protein_name` | CAZy's protein name, with all the names it gives in brackets |
| `ec_numbers` | EC numbers CAZy assigns; empty for most structure-only entries |
| `organism` | Organism as CAZy writes it, with strain |
| `genbank_ids` | GenBank protein accessions |
| `uniprot_ids` | UniProt accessions, including ones UniProt has since deleted |
| `pdb_codes` | PDB codes CAZy lists |
| `references` | Links to the characterization papers (characterized page only) |
| `in_cazy_characterized`, `in_cazy_structure` | `yes` / `no`: which of the two CAZy pages list the protein |
