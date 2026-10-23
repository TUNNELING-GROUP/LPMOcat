# Hydrophobicity scales

Hydrophobicity is the one property that is not a pipeline step. It needs no calculation on
the structure: every residue type has a fixed value on each scale, so the web app looks the
values up when it shows a protein. The table is
`lpmocat_web/catalogue/data/hydrophobicity_scales.csv`.

## Where the values come from

The values are those tabulated in the
[CLC Genomics Workbench manual (QIAGEN)](https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/650/Hydrophobicity_scales.html),
which collects the published scales on one page:

| Scale | Reference |
| --- | --- |
| Kyte–Doolittle | Kyte & Doolittle, J Mol Biol 1982, [doi:10.1016/0022-2836(82)90515-0](https://doi.org/10.1016/0022-2836(82)90515-0) |
| Hopp–Woods | Hopp & Woods, Mol Immunol 1983, [doi:10.1016/0161-5890(83)90029-9](https://doi.org/10.1016/0161-5890(83)90029-9) |
| Cornette | Cornette et al., J Mol Biol 1987, [doi:10.1016/0022-2836(87)90189-6](https://doi.org/10.1016/0022-2836(87)90189-6) |
| Eisenberg | Eisenberg et al., J Mol Biol 1984, [doi:10.1016/0022-2836(84)90309-7](https://doi.org/10.1016/0022-2836(84)90309-7) |
| Rose | Rose et al., Science 1985, [doi:10.1126/science.4023714](https://doi.org/10.1126/science.4023714) |
| Janin | Janin, Nature 1979, [doi:10.1038/277491a0](https://doi.org/10.1038/277491a0) |
| Engelman (GES) | Engelman et al., Annu Rev Biophys Biophys Chem 1986, [doi:10.1146/annurev.bb.15.060186.001541](https://doi.org/10.1146/annurev.bb.15.060186.001541) |

The scales differ in range and in what they measure. Hopp–Woods is a hydrophilicity scale:
its high values are the polar residues, the opposite of the others.

## How the app uses them

- **Per residue** (compare page): the Hydrophobicity plot and the 3D colouring give each
  residue its value on the chosen scale, Kyte–Doolittle by default.
- **Per protein** (catalogue overview): the mean over a set of residues: all, the
  surface-exposed ones (relative solvent accessibility ≥ 25 %, [step 5](../steps/05-dssp-sasa.md)),
  the [loop regions or the binding surface](../steps/07-core-fold.md). With Kyte–Doolittle over
  all residues this is the GRAVY value.
- **Similarity heatmap** (compare page): proteins compared position by position on the
  Kyte–Doolittle scale.

## File format

One row per amino acid, one column per scale:

| Column | Content |
| --- | --- |
| `aa` | One-letter code |
| `aa full` | Name |
| `aa classes` | Residue class (hydrophobic, polar, positive, negative, special); not used by the app |
| `Kyte-Doolittle` … `Engelman` | The value on each scale |
