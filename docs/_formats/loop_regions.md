JSON object, one key per structure. All ranges are `"first-last"`
[sequence positions](../conventions.md#numbering), not residue numbers.

```json
{
  "8RRY": {
    "loop_regions": {"L2": ["7-69"], "L3": ["84-88"], "LS": ["122-129"], "L8": ["154-157"]},
    "binding_surface": {"other": ["1-1"], "L2": ["24-24", "27-28", "30-31", "33-33", "55-55"],
                        "L3": ["84-88"], "LS": ["126-126"], "L8": ["155-156"], "β9": ["158-158"]}
  }
}
```

| Field | Meaning |
| --- | --- |
| `loop_regions` | Loop name → a list with one range. A loop the structure lacks is absent |
| `binding_surface` | Region → the runs of consecutive binding-surface residues in it. Regions are loop names, core strands (`β1` … `β9`) and `other` for residues in neither (typically the N-His). Null when no binding surface was found |
