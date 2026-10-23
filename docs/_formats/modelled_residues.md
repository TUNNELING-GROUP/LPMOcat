JSON object, one key per structure with rebuilt gaps.

| Field | Meaning |
| --- | --- |
| `source`, `chain`, `output` | Input file, chain and output file name |
| `method` | MODELLER version and settings, e.g. `MODELLER 10.8 AutoModel, observed atoms fixed, best DOPE of 10 models (seed -8123)` |
| `filled` | Rebuilt gaps: `{start, end, length, sequence}` with residue numbers |
| `skipped` | Gaps left open, with a `reason` (`longer than 30 residues`, `unknown residue type`) |
| `n_term_missing`, `c_term_missing` | Unobserved terminal residues, never rebuilt |
| `dope` | DOPE score of the model kept |
| `fixed_max_dev` | Largest CA shift of an observed residue in Å (0 expected) |
| `junctions` | C–N bond lengths in Å where rebuilt and observed residues meet (about 1.33) |
| `min_contact` | Closest heavy-atom distance in Å between a rebuilt residue and a non-neighbouring residue or ligand |
