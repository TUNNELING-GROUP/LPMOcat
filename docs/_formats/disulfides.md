JSON object, one key per structure (every aligned structure, also those without a bond).

```json
{
  "8RRY": [[41, 49, 2.02], [145, 162, 2.03]],
  "2YOY": []
}
```

| Field | Type | Meaning |
| --- | --- | --- |
| key | string | Structure identifier (`struct_id`) |
| value | list | One item per disulfide bond, in the order of the first cysteine; empty when the structure has none |
| item `[0]`, `[1]` | integer | [Residue numbers](../conventions.md#numbering) of the two cysteines, the earlier one in the chain first |
| item `[2]` | number | SG–SG distance in Å, rounded to 0.01 |
