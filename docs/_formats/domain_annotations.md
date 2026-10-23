JSON object, one key per UniProt accession (the active accessions of every protein, and the
accession of each AlphaFold DB model). Positions are UniProt sequence positions.

| Field | Meaning |
| --- | --- |
| `uniprot_id` | The accession |
| `domains` | List of `{source, accession, name, start, end}`: every InterPro member-database match (`pfam`, `interpro`, `cdd`, `cathgene3d`, `panther`, …) and UniProt's own `Domain` features (`source` `uniprot`) |
| `signal_peptide` | `{start, end}` of UniProt's `Signal` feature, or null |
| `chain` | `{start, end, description}` of UniProt's `Chain` feature (the mature protein), or null |
| `catalytic_domain` | `{start, end, pfam_id, name}` of the match to PF03443, PF03067 or PF22810, or null |
| `annotation_from` | Present when a deleted accession took the annotations of an active entry with the identical sequence: that entry's accession |
