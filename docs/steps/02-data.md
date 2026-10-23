# Step 2 · Data

Step 2 gives every protein of step 1 one structure and collects what the databases hold
about it: UniProt accessions, domain annotations, EC numbers, taxonomy, and the substrates and
oxidised positions that papers report.

| Name for `--steps` | `fetch` or `2` |
| --- | --- |
| Code | `pipeline.py` `step_fetch`, `fetch_structures.py`, `taxonomy.py`, `activities.py` |
| Reads | [Step 1](01-proteins.md): `cazy_lpmos.csv`, `other_lpmos.csv` · RCSB PDB, AlphaFold DB, UniParc, UniProt, InterPro, CAZyme3D, NCBI Taxonomy, Europe PMC · the [BRENDA release file](../running.md#brenda) · `curation/activities.csv`, `curation/models/` |
| Writes | `02_data/`: structure files, `pdb_selection.csv`, `uniprot_ids.json`, `uniprot_references.json`, `domain_annotations.json`, `taxonomy.json`, `activities.json` |
| Options | `--families` processes only those families |
| Run time | 15–20 min for about 170 entries (network-bound) |

## 2a–2d · One structure per protein

Entries are taken in the order CAZy's lists, then the other lists. Each gets the first
structure found:

| | Source | Rule | `struct_id` |
| --- | --- | --- | --- |
| 1e | Own model | An entry of `curation/local_models.csv` takes its model file as it is | `AF3-<code>` |
| 2a | RCSB PDB | The best of the PDB entries the protein lists (below) | PDB code |
| 2b | AlphaFold DB | The model of the first UniProt accession the entry lists that has one | `AF-<accession>` |
| 2c | AlphaFold DB via UniParc | The model of a UniProt entry with the same sequence as the entry's GenBank protein (active entries first, Swiss-Prot first, the entry's species first; up to 20 tried) | `AF-<accession>` |
| 2d | CAZyme3D | The model for the entry's GenBank protein | `C3D-<GenBank ID>` |
| | none | The entry is kept in `pdb_selection.csv` with `source` `NONE` | |

An experimental structure comes first because it is a measurement; AlphaFold DB before
CAZyme3D because its models are tied to a UniProt entry, whose annotations then apply in the
model's own numbering. CAZyme3D covers the proteins without any UniProt entry.

Two cases give no structure of their own:

- **Duplicate (`DUP`)**: in 2c the GenBank protein has the sequence of a UniProt accession
  another CAZy entry lists. The entry is recorded with `duplicate_of`.
- **Strain variant (`VAR`)**: an entry of the other lists marked `variant_of` in step 1,
  whose CAZy entry already has a structure. Its accession is added to that protein
  (`strain_variants`). If the CAZy entry has no structure, the variant is processed normally.

An own model is converted on import: the first model, the protein as chain A with its
B-factor column (pLDDT), and metal ions moved into chain A, since step 3 keeps one chain.

### 2a · Choosing the PDB entry

An entry listed in `curation/structures.csv` wins whenever the protein lists it: only the
pinned entries then compete. This is for what the rule cannot see, such as the high-dose end
of an X-ray dose series winning on resolution. Otherwise every chain of every listed entry is
ranked by these criteria, in order:

1. **Wild type** before an entry with engineered mutations (`pdbx_mutation`), whatever the
   resolution. A mutant is not the protein; it is taken only when no wild-type entry exists.
2. **Method**: X-ray, then neutron, electron microscopy, NMR, other.
3. **Active-site metal**: copper (ligand `CU` or `CU1`), then another metal that can occupy
   the brace (Zn, Ni, Co, Mn, Fe, Mo, Mg), then none. This counts only for entries within
   **1.0 Å** (`MAX_RESOLUTION_PENALTY`) of the best resolution among the wild-type entries.
4. **Resolution**, lower first.
5. **Coverage**: modelled residues / residues of the sample sequence, higher first.
6. **Chain length**, longer first.

The metal ranks before resolution because the catalogue measures the active site relative
to the copper: brace geometry, second-sphere distances, the binding plane. The 1.0 Å window
keeps a much poorer structure from winning on its metal alone. The metal is read from the
entry's ligand list, which cannot tell where in the structure it sits; the metal actually
kept is decided in [step 3g](03-preparation.md) and recorded in [step 6f](06-catalytic-domain.md).

The file is downloaded in PDB format, or as mmCIF when the entry has no PDB-format file.

## 2e · Identical sequences

Structures whose sequences are identical (read from the ATOM records, from the first
histidine on, so that a signal peptide makes no difference) are one protein listed twice.
One is kept, by: source (PDB, AlphaFold, CAZyme3D, own model), resolution, an entry of
CAZy's lists over the other lists, an entry on CAZy's characterized page, a UniProt accession
the entry lists over one found through UniParc.

## 2f · UniProt accessions

Each protein gets every accession that belongs to it: those its entry lists, the accession of
its model, the accessions of merged strain variants, and every *active* UniProt entry with
the identical sequence (UniParc). Deleted accessions that the entry lists are kept, since
CAZy still shows them; deleted identical-sequence entries are not, as they come from other
genomes and carry no data. Each accession is stored with its status, organism and EC numbers
in `uniprot_ids.json`.

`uniprot_references.json` holds the papers the entries cite for a curated scope (function,
catalytic activity, structure, cofactor, mutagenesis, …). Genome-sequencing citations have no
such scope and are left out, as are citations marked as describing an incorrect function.

## 2g · Domain annotations

For every active accession and every AlphaFold DB model's accession: all InterPro
member-database matches, and UniProt's signal peptide, chain and domain features. The Pfam
match to PF03443, PF03067 or PF22810 is stored as the catalytic domain. An accession without
any annotation (a deleted entry) takes those of the first active entry with the identical
sequence, which has the same residue numbering; `annotation_from` records it.

These annotations drive the N-terminal trim ([3h](03-preparation.md)), the catalytic domain
end ([6b](06-catalytic-domain.md)) and the CBM assignment in the web app.

## 2h · EC numbers

The first that applies: the EC numbers of the entry's own list (`ec_source` `CAZy`,
`Swiss-Prot` or `BRENDA`); those of the protein's reviewed UniProt entries (`UniProt
<accessions>`); the family default, the most common EC number string among the family's CAZy
entries (`family default`). CAZy gives no EC number for most structure-only entries, which is
what the two fallbacks are for.

## 2i · Taxonomy

Each organism name is resolved in NCBI Taxonomy (E-utilities) by the first query that returns
exactly one taxon: the full name as a scientific name, then as any name (synonyms), then its
first two words (the species), then its genus. The lineage and the species' taxid are stored.
Requests are spaced 0.34 s apart (NCBI allows three per second without a key) and
`NCBI_API_KEY` is used when set.

```bash
python lpmocat_pipeline/steps/taxonomy.py [--data-dir DIR]     # refill taxonomy.json only
```

## 2j · Substrates and oxidised positions

Four sources are read, and every claim is stored with its source, its text and its papers.
No source overrides another here; what the catalogue shows is decided at the
[Django import](../django-import.md).

| Source | What is taken |
| --- | --- |
| CAZy | The EC numbers of the protein's CAZy entries, as activity on the EC's substrate (chitin, cellulose, starch). CAZy lists its papers per entry, so they are not tied to one EC number |
| Swiss-Prot | From reviewed entries: catalytic activities with ECO:0000269, and the function text sentence by sentence, each sentence with the PubMed IDs it cites. A sentence yields substrates, substrates without activity (negations), and oxidised positions. Kinetics, regulation, modifications and mutants become notes |
| BRENDA | Substrate → product rows of EC 1.14.99.53–56, "no activity" remarks and notes, each with its papers |
| Literature | `curation/activities.csv`: claims none of the three holds, one row per claim, primary papers only |

**Reading text.** Substrate names are matched against an ordered list of patterns
(`SUBSTRATE_PATTERNS`: oligosaccharides and derivatives before their polymers) and mapped to
a category, keeping the name as written as the `form`. Oxidised positions are read from
phrases such as "C1-oxidized", "aldonic acid" (C1), "4-keto", "gem-diol" (C4). UniProt's
family template "C1 or C4" is ignored, because it names no position.

**BRENDA proteins.** A BRENDA protein is the catalogue protein with its UniProt accession.
BRENDA also pools the enzymes of an organism under records without an accession; such a row
is assigned to the catalogue protein of that species whose names contain a specific BRENDA
synonym (e.g. `NcLPMO9C`), or, without such a synonym, to the only protein of the species
that cites the row's paper. Rows that cannot be assigned go to `unmatched`.

**Disagreement.** The per-substrate summary is `conflict` only when the same kind of
substrate is reported both active and inactive. Activity on β-chitin and none on α-chitin, or
on amorphous but not crystalline cellulose, is not a conflict: forms are classed by kind
(`form_class`) and compared within a kind.

**References.** A PubMed ID must be 1–8 digits; a malformed one in BRENDA keeps the claim
under `BRENDA:<ec>:<number>`, and such a reference is merged with the `PMID:` reference of the
same title.

```bash
python lpmocat_pipeline/steps/activities.py [--data-dir DIR] [--brenda FILE] [--no-csv]
```

Run on its own, the script rewrites `activities.json` and also writes
`02_data/activities_review.csv`, the same claims one per row for reading through.

## Parameters

| Constant | Value | Meaning |
| --- | --- | --- |
| `MAX_RESOLUTION_PENALTY` (`fetch_structures.py`) | `1.0` Å | Resolution window in which the metal decides (2a) |
| `COPPER_LIGANDS` | `CU`, `CU1` | Ligand codes counted as copper (2a) |
| `METAL_LIGANDS` | copper + `ZN`, `NI`, `CO`, `MN`, `FE`, `MO`, `MG` | Ligand codes counted as a brace metal (2a) |
| `LPMO_PFAM_IDS` | PF03443, PF03067, PF22810 | Pfam families of the catalytic domain (2g) |
| `CURATED_REFERENCE_SCOPES` | function, catalytic activity, X-ray crystallography, structure, biophysicochemical, cofactor, methylation, mutagenesis, substrate | Citation scopes kept (2f) |
| `REQUEST_INTERVAL` (`taxonomy.py`) | `0.34` s | Pause between NCBI requests (2i) |
| `LPMO_EC_NUMBERS` (`activities.py`) | 1.14.99.53–56 | EC numbers read from BRENDA (2j) |

## Outputs

Structure files are saved as `02_data/<struct_id>.pdb` (`.cif` where RCSB has no PDB-format
file).

### pdb_selection.csv

--8<-- "_formats/pdb_selection.md"

### uniprot_ids.json

--8<-- "_formats/uniprot_ids.md"

### uniprot_references.json

--8<-- "_formats/uniprot_references.md"

### domain_annotations.json

--8<-- "_formats/domain_annotations.md"

### taxonomy.json

--8<-- "_formats/taxonomy.md"

### activities.json

--8<-- "_formats/activities.md"

## What the log tells you

One line per entry with its outcome (`LOC`, `PDB`, `AF`, `C3D`, `DUP`, `VAR`, `SKIP`), then
per sub-step: the dropped duplicates, the accession counts by status and the proteins without
an active entry, the accessions without annotation, where the EC numbers came from, how
organism names were matched, and for 2j the coverage per source, the conflicts and the
BRENDA rows that could not be assigned.

## Known limits

- Every run queries the databases again and rebuilds all files; there is no record of what
  changed in UniProt since the last run.
- Without a BRENDA file the BRENDA source is skipped with a warning.
- A model whose UniProt entry is deleted and has no active identical-sequence entry has no
  annotation; step 3h then relies on SignalP.
