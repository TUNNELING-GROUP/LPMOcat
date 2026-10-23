# Curation files

The pipeline handles structures by general rules. Where a rule gives the wrong answer for one
protein, or the databases lack something, the correction is a row in a CSV file under
`lpmocat_pipeline/curation/`. Each row carries its reason, so every exception can be traced.
The files are read from next to the code, independent of `--data-dir`. Lines starting with
`#` are comments in `extra_proteins.csv` and `local_models.csv`.

| File | Read by | Purpose |
| --- | --- | --- |
| [`extra_proteins.csv`](#extra_proteinscsv) | step 1d | Characterized LPMOs the automatic lists miss |
| [`local_models.csv`](#local_modelscsv) | steps 1e, 2 | In-house database entries, with their AlphaFold 3 models |
| [`structures.csv`](#structurescsv) | step 2a | PDB entries chosen by hand |
| [`n_terminus.csv`](#n_terminuscsv) | step 3h | First residue of a structure |
| [`c_terminus.csv`](#c_terminuscsv) | step 3i | Last residue of a structure |
| [`activities.csv`](#activitiescsv) | step 2j | Activity claims from the literature |
| [`protein_codes.csv`](#protein_codescsv) | Django import | Short names of proteins |
| [`substrate_reviews.csv`](#substrate_reviewscsv) | nothing | A reading list for curators |

## extra_proteins.csv

Unreviewed UniProt entries with published experimental results. Step 1b admits only
reviewed (Swiss-Prot) entries; this file is the exception list.

| Column | Meaning |
| --- | --- |
| `accession` | UniProt accession |
| `protein` | The name used in the literature; added to the entry's names |
| `family` | CAZy family |
| `pmid` | PubMed ID of the paper with the experimental results |
| `reason` | Why the entry is added |

Current entries:

| Accession | Protein | Family | PMID | Reason |
| --- | --- | --- | --- | --- |
| G0SDM5 | CtPMO1 | AA9 | 29991963 | Unreviewed entry with published results (Chen et al. 2018). It is the genome entry with the full sequence. The paper's own GenBank protein (AGY80105.1, UniProt V9PD45) is partial: it starts after the N-His, so a model of it would have no histidine brace |
| A0A896IH95 | CgAA15A | AA15 | 33421883 | Unreviewed entry with published results (Franco Cairo et al. 2021); the paper's GenBank protein QSC42436.1 |

## local_models.csv

Entries of the in-house database: characterized LPMOs that enter the catalogue with a model
made with AlphaFold 3. Each row points to a model file under `curation/models/`.

| Column | Meaning |
| --- | --- |
| `struct_id` | Identifier of the model, `AF3-<code>` |
| `protein` | Protein name |
| `family`, `organism`, `kingdom` | As in the CAZy list |
| `accession` | Optional: the UniProt entry the model stands in for. The protein then keeps that entry's EC numbers, papers and activities, and only the structure is replaced |
| `file` | Model file, relative to `lpmocat_pipeline/` (mmCIF or PDB) |
| `pmid` | Paper characterizing the protein |
| `reason` | Why the entry is added |

## structures.csv

Overrides the choice of the PDB entry in [step 2a](../steps/02-data.md#2a-choosing-the-pdb-entry):
when a protein lists one of these entries, only the listed entries compete.

| Column | Meaning |
| --- | --- |
| `pdb_entry` | PDB code |
| `reason` | Why |

| `pdb_entry` | Reason |
| --- | --- |
| 7PTZ | LsAA9A: the standard high-resolution structure (1.09 Å, copper). The rule would take 9TDI (1.06 Å), the highest-dose entry (2.8 MGy) of an X-ray dose series, whose copper site is likely photoreduced and radiation-damaged; the series' low-dose entries are 1.5–2.0 Å |

## n_terminus.csv

Overrides the N-terminal trim of [step 3h](../steps/03-preparation.md#3h-n-terminal-trim-to-the-n-his).

| Column | Meaning |
| --- | --- |
| `struct_id` | Structure |
| `first_residue` | Residue number to start at, or `keep` to leave the N-terminus as it is |
| `reason` | Why |

| `struct_id` | `first_residue` | Reason |
| --- | --- | --- |
| 6RS9 | `keep` | LsAA9B: the deposited chain starts at an arginine and has no N-His, so no histidine brace. The N-terminus stays as deposited instead of being cut at the first histidine further in |

## c_terminus.csv

Overrides the C-terminal tag removal of [step 3i](../steps/03-preparation.md#3i-c-terminal-expression-tags).

| Column | Meaning |
| --- | --- |
| `struct_id` | Structure |
| `last_residue` | Residue number to end at, or `keep` to leave the C-terminus as it is |
| `reason` | Why |

The file has no entries: the rules of 3i are right for every current structure.

## activities.csv

Substrate and regioselectivity claims that CAZy, Swiss-Prot and BRENDA do not hold
([step 2j](../steps/02-data.md#2j-substrates-and-oxidised-positions), source "Literature").
One row per claim. Only primary papers are cited, never reviews. Rows naming an unknown
structure, substrate, result or position are reported and skipped.

| Column | Meaning |
| --- | --- |
| `struct_id`, `protein` | The catalogue protein |
| `substrate` | Substrate category (as in `SUBSTRATE_PATTERNS`) |
| `form` | The substrate as the paper describes it |
| `result` | `active` or `inactive` |
| `oxidation` | `C1`, `C4`, or both |
| `pmid`, `doi` | The primary paper |
| `evidence` | Where in the paper the claim was checked |

| Structure | Protein | Substrate (form) | Result | Position | PMID | Evidence |
| --- | --- | --- | --- | --- | --- | --- |
| 6Z5Y | PiAA17C | pectin (demethylated homogalacturonan) | active | C4 | 34385392 | Cleaves demethylated pectin at C4; checked through Turella et al. 2025, which cites the paper for this |
| 6Z5Y | PiAA17C | pectin (methyl-esterified homogalacturonan) | active | C4 | 40947185 | Abstract: C4-oxidized oligogalacturonides generated |
| 8S6S | PsAA17A | pectin (demethylated) | active | C4 | 40216756 | Results and Supplementary Fig. 17: both LPMOs cleave the C4 position of demethylated pectin |
| 8S71 | PsAA17B | pectin (demethylated) | active | C4 | 40216756 | As for PsAA17A |
| AF-A0A2I6QAZ5 | PcAA14A | xylan (xylan-coated cellulose fibres) | active | C1 | 29377002 | C1 oxidation of the xylan coating cellulose fibres, checked in the paper. A later study (PMID 37418595) found no activity on re-evaluation |

## protein_codes.csv

A curated short name wins over the one `update_protein_codes` would derive
([Django import](../django-import.md)).

| Column | Meaning |
| --- | --- |
| `struct_id` | Structure |
| `code` | The protein code to use |
| `reason` | Why |

| `struct_id` | Code | Reason |
| --- | --- | --- |
| AF-A0A2I6QAZ5 | PcAA14A | The name used in the literature (Couturier et al. 2018); CAZy lists PcoAA14A first |
| AF-Q2US83 | AoAA9A | The name used in the literature (Frandsen et al. 2021); Swiss-Prot gives this *A. oryzae* entry the short name AgAA9A, which belongs to the *Armillaria gallica* protein |
| AF-G0SDM5 | CtPMO1 | The name used in the literature (Chen et al. 2018); a code built from UniProt's organism name would be TtAA9, like the *T. terrestris* and *T. thermophilus* proteins |

## substrate_reviews.csv

Not read by the pipeline. It lists reviews and one research paper with tables of LPMO
substrates and regioselectivity, as a starting point for finding primary papers to add to
`activities.csv`. Reviews themselves are not cited in the catalogue.

Columns: `first_author`, `year`, `journal`, `title`, `doi`, `open_access`, `table_checked`,
`what_it_holds`, `note`.

| Paper | Journal | DOI | Table checked |
| --- | --- | --- | --- |
| Zhou 2020 | Bioresources and Bioprocessing | 10.1186/s40643-020-0300-6 | yes |
| Frommhagen 2018 | Frontiers in Microbiology | 10.3389/fmicb.2018.01080 | yes |
| Frandsen 2016 | IUCrJ | 10.1107/s2052252516014147 | yes |
| Forsberg 2023 | Essays in Biochemistry | 10.1042/ebc20220162 | yes |
| Morgenstern 2014 | Briefings in Functional Genomics | 10.1093/bfgp/elu032 | yes |
| Frandsen 2021 (research paper) | Journal of Biological Chemistry | 10.1074/jbc.ra120.015545 | yes |
| Vandhana 2022 | New Phytologist | 10.1111/nph.17921 | no |
| Batka 2026 | Chemical Reviews | 10.1021/acs.chemrev.6c00006 | no |
| Munzone 2024 | Nature Reviews Chemistry | 10.1038/s41570-023-00565-z | no |
| Agostoni 2017 | Microbiology and Molecular Biology Reviews | 10.1128/mmbr.00015-17 | no |
| Courtade 2019 | Advances in Experimental Medicine and Biology | 10.1007/978-981-13-7318-3_6 | no |
| Li 2025 | Bioresource Technology | 10.1016/j.biortech.2025.133022 | no |
| Vu 2016 | Cellular and Molecular Life Sciences | 10.1007/s00018-016-2251-9 | no |
| Zhang 2020 | Applied Microbiology and Biotechnology | 10.1007/s00253-020-10467-5 | no |
| Monclaro 2017 | International Journal of Biological Macromolecules | 10.1016/j.ijbiomac.2017.04.077 | no |
