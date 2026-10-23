JSON object, one key per structure with a binding surface.

| Field | Meaning |
| --- | --- |
| `plane.normal`, `plane.offset` | The plane in the aligned coordinate frame: points **x** with `normal · x = offset`; the normal points away from the protein |
| `copper.source` | Residue name of the metal used (`CU`, `CU1`, `ZN`, `MG`, …) or `pseudo-Cu`, following step 6f's `copper_source` |
| `copper.coord` | Copper position used |
| `copper.depth_below_plane` | Distance of the copper below the plane, Å |
| `copper.shallowest_possible` | Smallest copper depth any tangent plane reaches |
| `copper.recessed` | True when no tangent plane comes within 3 Å of the copper |
| `contacts_around_copper` | 0–1: how evenly the binding residues surround the copper (1 = all around, 0 = all on one side) |
| `residues_above_plane` | Residues sticking through the plane (empty with the default settings) |
| `n_residues` | Number of binding-surface residues |
| `settings` | The parameter values used |
| `residues` | One item per binding-surface residue: `position` (sequence position), `residue_number`, `insertion_code`, `name`, `distance` (of its outermost atom below the plane, Å), `rsa` (relative SASA, %), `region` (loop or core strand, empty if neither) |
