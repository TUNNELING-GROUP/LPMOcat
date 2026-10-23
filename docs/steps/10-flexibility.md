# Step 10 · Flexibility

Step 10 collects three per-residue measures of flexibility, computed on the prepared
structures of step 3.

| Name for `--steps` | `flexibility` or `10` |
| --- | --- |
| Code | `pipeline.py` `step_flexibility` |
| Reads | [Step 3](03-preparation.md): `preprocessed/<id>.pdb` · [step 2](02-data.md): `pdb_selection.csv` |
| Writes | `10_flexibility/flexibility_raw.json`, `cabsflex_results/<id>/` |
| External tools | ProDy (required), CABS-flex (only with `--cabsflex`) |
| Options | `--cabsflex` runs the missing and outdated CABS-flex simulations |
| Run time | under a minute without simulations; CABS-flex about 3 min per structure, in parallel (about 55 min for 160 structures on 19 worker processes) |

## 10a · B-factors and pLDDT

The B-factor column of the CA atoms is read as it is: **B-factors** for experimental
structures, **pLDDT** for models (AlphaFold DB, CAZyme3D, own models). Residues rebuilt by
MODELLER in step 3a (occupancy 0) get no B-factor (null), because their value is not a
measurement.

The two are different quantities and are stored under different keys. B-factors depend on
resolution and refinement and are not comparable between crystals; pLDDT is a confidence, not
a motion.

## 10b · Gaussian network model

ProDy's `GNM` on the CA atoms with its default settings (cutoff 10 Å, spring constant 1),
using the 20 slowest modes (fewer for a chain shorter than 21 residues). The value stored is
the mean-square fluctuation of each residue from those modes. The scale is arbitrary: it
compares residues within a protein, not proteins with each other.

## 10c · CABS-flex

A coarse-grained simulation of each structure:

```text
python3 -m CABS flex -i <id>.pdb --work-dir cabsflex_results/<id> --csv-output A --no-progress-bar
```

with CABS-flex's default protocol. The RMSF is computed from the CA trajectory CABS-flex
writes (`output_pdbs/replica.pdb`): all frames are superposed with ProDy and the fluctuation
of each CA is taken. This is the one measure every protein has on a comparable scale (Å),
which is why the web app uses it for per-loop flexibility.

**Reuse.** A run is tied to the structure it simulated by a fingerprint stored in
`cabsflex_results/<id>/input.sha1`: the SHA-1 of the atom names, residue names and numbers
and coordinates of the `ATOM` records. A run whose fingerprint matches the current prepared
structure is reused. Missing or outdated runs are only *reported* unless `--cabsflex` is
given; then they are rerun, on all cores but one. The fingerprint ignores atom serial
numbers, B-factors and the metal, so changes that do not affect the simulated protein do not
trigger a rerun.

## Outputs

### flexibility_raw.json

--8<-- "_formats/flexibility_raw.md"

`cabsflex_results/<id>/` keeps each CABS-flex run (trajectory, CSV output, `input.sha1`).

## Known limits

- CABS-flex's final all-atom reconstruction may fail when its helper program is missing;
  this does not matter, since the RMSF comes from the CA trajectory written before that.
- A structure without a current CABS-flex run has no `cabsflex` key; the log lists these.
- A structure for which the GNM cannot be built has no `gnm` key.
