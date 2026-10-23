JSON object, one key per aligned structure whose residue count matches its alignment row.
The keys `his1` and `his2` stand for the [N-His and the Internal-His](../conventions.md#brace-histidines).

```json
{
  "8RRY": {
    "his1": 28, "his2": 114,
    "his1_column": 2, "his2_column": 319,
    "his_n_distance": 2.98,
    "method": "metal", "metal": "CU",
    "pseudo_cu": [2.9026, -4.6502, -12.2004],
    "copper": [2.852, -4.596, -12.269], "copper_source": "metal"
  }
}
```

| Field | Type | Meaning |
| --- | --- | --- |
| `his1` | integer or null | [Residue number](../conventions.md#numbering) of the N-His |
| `his2` | integer or null | Residue number of the Internal-His; always null when `his1` is null |
| `his1_column`, `his2_column` | integer or null | Their alignment columns, 1-based |
| `his_n_distance` | number or null | Smallest distance in Å between any nitrogen of the N-His and any nitrogen of the Internal-His (backbone N included); null unless both were found |
| `method` | string | `metal`: the two histidines closest to the active-site metal; `alignment`: found at the brace columns of the metal-bound structures |
| `metal` | string or null | Element of the active-site metal (`CU`, `ZN`, `NI`, …), null for `method` `alignment`. The element, not the PDB ligand code: `CU1` and `CU` are both `CU` |
| `pseudo_cu` | list or null | `[x, y, z]` of the pseudo-copper in the aligned coordinate frame, Å; given for every structure with a brace, also those with a metal |
| `copper` | list or null | `[x, y, z]` of the [copper position](../steps/06-catalytic-domain.md#copper-position) that steps 7c, 8 and 12 use, Å: the metal's coordinates or `pseudo_cu` |
| `copper_source` | string or null | `metal` or `pseudo-Cu`; null when there is neither |
