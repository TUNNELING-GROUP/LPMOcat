JSON object, one key per aligned structure. The catalytic domain runs from the first residue
of the prepared structure (the N-His) to `end_residue`, inclusive.

```json
{
  "8RRY": {
    "end_residue": 197,
    "family": "AA10",
    "source": "consensus"
  },
  "AF-W4K8M0": {
    "end_residue": 243,
    "family": "AA9",
    "source": "disulfide",
    "uncapped_end_residue": 230,
    "original_source": "annotation",
    "end_before_disulfide": 210
  }
}
```

| Field | Type | Present | Meaning |
| --- | --- | --- | --- |
| `end_residue` | integer | always | [Residue number](../conventions.md#numbering) of the last residue of the catalytic domain |
| `family` | string | always | CAZy family of the structure (`AA9` … `AA17`), from `pdb_selection.csv`; `unknown` if the structure is not listed there |
| `source` | string | always | The rule that set `end_residue` last: `annotation`, `consensus`, `fallback` (6b), `capped` (6c or 6d) or `disulfide` (6e) |
| `original_source` | string | when 6c, 6d or 6e moved the end | The 6b source (`annotation`, `consensus` or `fallback`) |
| `uncapped_end_residue` | integer | when 6c or 6d moved the end | The end 6b assigned |
| `end_before_disulfide` | integer | when 6e moved the end | The end after 6b–6d, before it was moved to the second cysteine |

In the example, 6b gave AF-W4K8M0 its annotated end 230, 6d moved it back to 210 (it lay in
a long coil), and 6e moved it forward to Cys243, which is bonded to Cys163.
