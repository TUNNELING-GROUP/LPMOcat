# Step 1 · Proteins

Step 1 lists the LPMOs the catalogue will hold: every protein CAZy lists as characterized or
with a structure, plus characterized LPMOs that CAZy's lists miss.

| Name for `--steps` | `scrape` or `1` |
| --- | --- |
| Code | `scrape_cazy.py` (1a), `other_sources.py` (1b–1e) |
| Reads | cazy.org, UniProt, UniParc, Europe PMC; the [BRENDA release file](../running.md#brenda); `curation/extra_proteins.csv`, `curation/local_models.csv` |
| Writes | `01_proteins/cazy_lpmos.csv`, `other_lpmos.csv` |
| Options | `--families AA9,AA10` limits both files to those families |
| Run time | about 20 s |

1b–1e can be rerun on their own against the existing CAZy table:

```bash
python lpmocat_pipeline/steps/other_sources.py [--data-dir DIR] [--brenda FILE]
```

## 1a · CAZy

For each family AA9, AA10, AA11, AA13, AA14, AA15, AA16 and AA17 the *characterized* and
*structure* pages of cazy.org are read. Entries of the two pages that share a GenBank
accession or a PDB code are one protein; `in_cazy_characterized` and `in_cazy_structure`
record where it was listed. Protein and organism names are reduced to ASCII.

CAZy is the starting list because it is the curated reference for these families: it
assigns the family and admits a protein only with published characterization or a structure.

## 1b–1e · Characterized LPMOs outside CAZy's lists

Each source needs experimental evidence tied to a paper, so that nothing enters on a
prediction alone.

| Sub-step | Source | An entry qualifies when |
| --- | --- | --- |
| 1b | Swiss-Prot | It is a reviewed UniProt entry with an LPMO EC number (1.14.99.53–56) whose catalytic activity carries the evidence code ECO:0000269 (experimental, published) with a PubMed reference |
| 1c | BRENDA | It is a protein of EC 1.14.99.53–56 with a UniProt accession and a substrate row citing a paper whose title or abstract (Europe PMC) names a genus of the protein's organism |
| 1d | `curation/extra_proteins.csv` | A curator lists it: an unreviewed UniProt entry with published experimental results |
| 1e | `curation/local_models.csv` | It is an entry of the in-house database, with its AlphaFold 3 model |

The genus check in 1c exists because BRENDA sometimes attaches a protein to a paper about
another organism. If Europe PMC does not answer, the paper counts as naming the organism and
the log lists it under `1c WARNING`, so a network failure never removes a protein silently.
Without a BRENDA file 1c is skipped.

A 1e row may name the UniProt `accession` its model stands in for. That accession is reserved
before 1b–1d run, so the protein is not added a second time, and the entry keeps the
accession's EC numbers, papers and activities while the model replaces only the structure.

### Family

In this order: the family given in the curation file; the family named in a reviewed entry's
similarity comment ("Belongs to the polysaccharide monooxygenase AA9 family"); the family of
the most similar CAZy-listed protein. Similarity is the best local alignment (BLOSUM62, gap
open −10, extend −0.5) against the sequences of the active UniProt entries CAZy lists, taking
the three best-scoring references per family. It must reach 40 % identity over at least 100
aligned residues and lie 5 percentage points above the best other family. An entry without a
clear family is reported for curation and not added.

Sequence similarity is used because domain annotations cannot separate the families: Pfam's
PF03067 matches AA10, AA11, AA13, AA15 and AA16 alike.

### What is dropped

- an entry CAZy already lists: by UniProt accession, PDB code or GenBank protein, or with
  the identical sequence (UniParc supplies the sequence of deleted accessions and of
  GenBank-only entries);
- a deleted UniProt entry, and a family outside AA9–AA17;
- within the other lists, an entry identical to one already added, or at least 95 %
  identical to one of the same species. Swiss-Prot entries are taken in the order of the
  number of papers they cite, so the better-documented entry stays.

### Strain variants

An entry at least 95 % identical (over at least 80 % of the shorter sequence) to a
CAZy-listed protein of the same species is the same protein from another strain. It is kept
with that protein's accession in `variant_of`. [Step 2](02-data.md) uses it only if the CAZy
entry gets no structure; otherwise it drops the variant and adds its accession to the CAZy
protein, so the variant's Swiss-Prot annotations and papers still reach it.

### Organism

UniProt's organism name without the strain. If that name, a synonym or UniProt's "common
name" matches a species among CAZy's organisms, CAZy's spelling is used, so that one species
has one name in the catalogue.

## Parameters

| Constant (`other_sources.py`) | Value | Meaning |
| --- | --- | --- |
| `LPMO_EC_NUMBERS` | 1.14.99.53, .54, .55, .56 | EC numbers searched in Swiss-Prot and BRENDA |
| `EXPERIMENTAL` | `ECO:0000269` | Evidence code required in 1b |
| `MIN_IDENTITY` | `0.40` | Lowest identity for a family by similarity |
| `MIN_ALIGNED` | `100` residues | Shortest alignment counted |
| `MIN_MARGIN` | `0.05` | Lead over the next family |
| `SAME_PROTEIN` | `0.95` | Identity from which two same-species entries are one protein |
| `LPMO_FAMILIES` (`scrape_cazy.py`) | AA9–AA11, AA13–AA17 | Families read from CAZy |

## Outputs

### cazy_lpmos.csv

--8<-- "_formats/cazy_lpmos.md"

### other_lpmos.csv

--8<-- "_formats/other_lpmos.md"

## What the log tells you

The number of proteins per CAZy page and family; for 1b–1e the number of candidates and
additions per source, every dropped entry grouped by reason, the strain variants, and the
entries whose family is unclear and need a curator's decision.

## Known limits

- cazy.org is read as HTML. A change of its page layout gives empty lists, reported as
  `0 proteins` or an error per family, not a failed step.
- With `--families`, both output files are overwritten with those families only.
