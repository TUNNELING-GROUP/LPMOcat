JSON object, one key per structure.

```json
{
  "8RRY": {
    "checked": "2026-09-23",
    "uniparc": ["UPI00000B5D1F"],
    "accessions": {
      "P9WEM3": {"source": "CAZy", "status": "reviewed",
                 "organism": "Serratia marcescens", "ec_numbers": ["1.14.99.-"]},
      "A0A1G5JNT6": {"source": "identical sequence", "status": "unreviewed",
                     "organism": "Serratia nematodiphila", "ec_numbers": []}
    }
  }
}
```

| Field | Meaning |
| --- | --- |
| `checked` | Date of the lookup |
| `uniparc` | UniParc identifiers of the sequence |
| `accessions` | One key per UniProt accession, the entry's own first |
| `source` | How the accession was found: the entry's list (`CAZy`, `Swiss-Prot`, `BRENDA`, `curated`, `local model`), `UniParc` (the model's accession, found by sequence), `CAZyme3D model`, `identical sequence` or `strain variant` |
| `status` | `reviewed`, `unreviewed`, `deleted`, or `unknown` when the lookup failed |
| `organism`, `ec_numbers` | From the UniProt entry; empty for deleted entries. A deleted entry also has `deleted_reason` |
