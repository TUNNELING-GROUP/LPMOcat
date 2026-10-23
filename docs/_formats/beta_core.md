JSON object with two keys.

```json
{
  "columns": {"b1": [3, 10], "b3": [242, 245], "b4": [255, 261]},
  "structures": {
    "8RRY": {
      "b1": {"range": "2-6", "sheet": "A", "pieces": ["2-6"]},
      "b3": {"range": "70-73", "sheet": "B", "pieces": ["70-73"]}
    }
  }
}
```

| Field | Meaning |
| --- | --- |
| `columns` | The alignment columns (1-based, first and last) of each core-strand block. `b2` has none: it is found by its pairing with `b1` |
| `structures` | Per structure, the core strands it has. A strand it lacks is absent |
| `range` | `"first-last"` [sequence positions](../conventions.md#numbering) from the strand's first piece to its last |
| `pieces` | The strand's separate runs of strand residues, as `"first-last"` sequence positions; usually one |
| `sheet` | DSSP sheet label of most of its residues, `-` if none |
