# Step 5 · DSSP and SASA

Step 5 computes the solvent accessibility and the secondary structure of every residue, on
the aligned structures.

| Name for `--steps` | `dssp_sasa` or `5` |
| --- | --- |
| Code | `sasa_dssp.py` |
| Reads | [Step 4](04-alignment.md): `alignment.fasta`, `aligned/<id>.pdb` |
| Writes | `05_dssp_sasa/sasa.json`, `secondary_structure.json`, `dssp.json` |
| External tools | mkdssp 4.4.11 (required: the step stops without it), Biopython |
| Run time | 1–2 min |

## 5a · Solvent accessibility

The solvent-accessible surface area (SASA) of each residue is computed with the
Shrake–Rupley algorithm (Biopython `ShrakeRupley`, default probe radius 1.4 Å and 100 points
per atom) on the **protein atoms only**. The metal is removed first, so that it does not bury
part of the histidine brace in some structures and not in others.

The **relative SASA** is the SASA divided by the maximum for that residue type, in per cent.
The maxima are the theoretical values of Tien et al. 2013 for the residue in a Gly-X-Gly
tripeptide (Biopython's `Wilke` scale). It is a whole-residue measure: backbone and side
chain together.

Two thresholds are applied to it elsewhere: **≥ 25 %** counts as surface exposed (the web
app; `SURFACE_RSA`, used here only for the log), and **≥ 10 %** makes a Tyr, Trp or His a
hole-hopping donor ([step 12](12-hole-hopping.md)).

## 5b · Secondary structure

mkdssp is run on the standard amino acids of each structure. The full DSSP result is kept
(`dssp.json`): the eight states, the sheet labels and the bridge partners with their ladder
labels, which [step 7](07-core-fold.md) needs to identify the core strands.

A three-state reduction is stored separately: **H** for H, G and I; **E** for E and B; **-**
for everything else. This is the reduction the 3D viewer (Mol\*) draws, so plots and
structures agree.

mkdssp is used, and not a reimplementation, because the DSSP β-bulge rule matters here: it
joins the two parts of a strand across a bulge, where simpler assignments break the strand
into pieces.

mkdssp is found through `MKDSSP`, then `PATH`, then the conda environment `dssp`. It runs in
a temporary directory that holds links to its mmCIF dictionary and to the bundled
`lpmocat_pipeline/tools/dssp/components.cif`, because the conda-forge build looks for these data files
only in the working directory.

## Parameters

| Constant | Value | Meaning |
| --- | --- | --- |
| `MAX_ASA` (`sasa_dssp.py`) | Biopython `residue_sasa_scales["Wilke"]` | Maximum SASA per residue type |
| `THREE_STATE` | H, G, I → H; E, B → E | Three-state reduction |
| `SURFACE_RSA` (`pipeline.py`) | `25` % | Surface-exposed threshold reported in the log |
| `DSSP_ENV` | `dssp` | Conda environment searched for mkdssp |

## Outputs

### sasa.json

--8<-- "_formats/sasa.md"

### secondary_structure.json

--8<-- "_formats/secondary_structure.md"

### dssp.json

--8<-- "_formats/dssp.md"

## Known limits

- A residue with an incomplete backbone is skipped by DSSP and gets `-`.
- The step fails if SASA or DSSP fails for any structure; the log names it.
