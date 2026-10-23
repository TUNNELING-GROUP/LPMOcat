JSON object, one key per aligned structure, sorted by identifier.

```json
{
  "8RRY": {
    "inputs": "5efe295eb797e73e",
    "paths": [
      {"rank": 1, "path": ["TYR39", "TRP108", "TRP178", "CU201"],
       "mean_residence_time": 0.042158403286,
       "donor_residue": "TYR39", "donor_distance": 16.31}
    ],
    "note": ""
  }
}
```

| Field | Meaning |
| --- | --- |
| `inputs` | Fingerprint of everything the result depends on; null for a structure that was not run |
| `paths` | Paths ranked by exact mean residence time, shortest first; null when none were found |
| `rank` | 1 = shortest mean residence time |
| `path` | Residues from the donor to the copper, as residue name + [residue number](../conventions.md#numbering). The last item is the copper: the metal's own name and number, or `CU999` for a pseudo-copper |
| `mean_residence_time` | EHPath's exact mean residence time of the charge on the path, in seconds |
| `donor_residue` | The first residue of the path |
| `donor_distance` | Distance from the donor's ring centroid to the copper, Å |
| `note` | Why there are no paths (`excluded: …`, `no histidine brace`), or remarks from the run (`no metal: pseudo-Cu at the brace nitrogens`) |

`<id>_paths.txt` is the human-readable form of the same result: all paths ranked by
approximate and by exact mean residence time. `EH_files/` holds EHPath's input (donor,
bridge and acceptor atoms as CSV), `EH_results/` its raw output per donor, and
`<id>_pseudocu.pdb` the structure with the pseudo-copper inserted.
