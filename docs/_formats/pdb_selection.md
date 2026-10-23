CSV, one row per entry of step 1, including entries without a structure (`source` `NONE`).
Later steps use only the rows with a structure. Multi-valued cells are separated by `; `.

| Column | Meaning |
| --- | --- |
| `family`, `protein_name`, `organism` | As in the step 1 lists; together they identify the entry |
| `source` | `PDB`, `AlphaFold`, `CAZyme3D`, `Local model` or `NONE` |
| `selected_pdb` | The structure identifier ([`struct_id`](../conventions.md#structure-identifiers)) used by every later step |
| `selected_chain` | Chain taken from the file (`A` for models) |
| `all_pdb_codes` | Every PDB code the entry lists |
| `uniprot_ids` | Every UniProt accession of the protein (2f), deleted ones included |
| `uniprot_active` | The accessions that are not deleted |
| `uniparc_ids` | UniParc identifiers of the protein's sequence |
| `all_af_models` | Accessions for which AlphaFold DB has a model |
| `pdb_code`, `chain_id`, `resolution`, `seq_length`, `modeled_count`, `coverage` | RCSB metadata of the selected entry and chain: resolution in Å, residues in the sample sequence, residues modelled, and their ratio. Empty for models (`coverage` 1.0) |
| `method` | `X-ray`, `NMR`, …; `AlphaFold`; `AlphaFold DB model of <accession>` or `CAZyme3D AlphaFold2 prediction`; `AlphaFold 3 (AlphaFold Server)` |
| `list_source` | `CAZy`, `Swiss-Prot`, `BRENDA`, `curated` or `local model` |
| `variant_of` | From step 1: the CAZy accession this entry is a strain variant of |
| `strain_variants` | Accessions of strain variants merged into this entry |
| `duplicate_of` | For `NONE` rows dropped as a duplicate or variant: the entry that holds the structure |
| `ec_numbers`, `ec_source` | EC numbers (2h) and where they come from: `CAZy`, `Swiss-Prot`, `BRENDA`, `UniProt <accessions>`, `family default` or empty |
| `genbank_id` | The GenBank protein a UniParc or CAZyme3D lookup started from |
| `uniprot_source` | `UniParc <id>` when the model's accession was found through UniParc |
| `plddt` | Mean pLDDT of a CAZyme3D model |
