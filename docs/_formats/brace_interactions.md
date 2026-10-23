JSON object, one key per aligned structure. Every field is null when it could not be
determined. Residues are given by [residue number](../conventions.md#numbering).

| Field | Meaning |
| --- | --- |
| `second_sphere_aromatic` | 8a: `{name, residue, column, method}`: `PHE` or `TYR`, its residue number, its alignment column (1-based) and `metal` or `alignment` |
| `interactions` | 8b: three lists, see below |
| `copper` | The copper position 8c and 8d used: `metal` or `brace nitrogens` (the pseudo-copper) |
| `second_sphere_motif` | 8c: `{motif, short, residues}`: the motif with spacing (`H-X8-Q-X-Y`), without (`H-Q-Y`), and its residues, each `{role, name, residue}` with `role` `first`, `gln_glu`, `middle` or `aromatic`. The first residue also has `found_by` (`stacked on His-2`, `sequence motif`, `hydrogen bond`); a Gln/Glu found by its distance to the copper has `cu_distance` |
| `distal_axial` | 8d, the shielding residue: `{name, residue, atom, distance, angle}` (Å and degrees from the brace-plane normal), or `{"name": null}` when the position is open |
| `brace_motifs` | 8e: `n_his` and `internal_his`, each `{motif, sequence, residues}`: the motif name or null, the one-letter sequence from the histidine on (10 and 5 residues, shorter at a chain break), and the motif residues `{role, name, residue, offset}` |

Items of `interactions`:

| List | Fields |
| --- | --- |
| `hydrogen_bonds` | `brace_res_id`, `other_res_name`, `other_res_id`, `donor_atom`, `acceptor_atom`, `distance` (Å) |
| `pi_pi_stacking` | `brace_res_id`, `brace_res_name`, `other_res_name`, `other_res_id`, `distance` (ring centroids, Å), `angle` (between ring planes, degrees) |
| `cation_pi` | `brace_res_id`, `brace_res_name`, `other_res_name`, `other_res_id`, `cation_atom`, `distance` (Å), `projection` (0–1, onto the ring normal) |
