# Step 3 · Preparation

Step 3 turns every selected structure into one clean file: a single protein chain that starts
at the N-His, with its active-site metal and nothing else. All later steps work on these
files, so what is removed here is absent from every calculation.

| Name for `--steps` | `prepare` or `3` |
| --- | --- |
| Code | `fill_gaps.py`, `modeller_worker.py` (3a), `preprocess_structures.py` (3b–3i) |
| Reads | [Step 2](02-data.md): structure files, `pdb_selection.csv`, `domain_annotations.json` · `curation/n_terminus.csv`, `curation/c_terminus.csv` · RCSB Chemical Component Dictionary |
| Writes | `03_preparation/filled/<id>.pdb`, `filled/modelled_residues.json`, `signal_peptides.json`, `signalp/`, `modified_residues.json`, `preprocessed/<id>.pdb` |
| External tools | MODELLER (3a) and SignalP 6 (3h), both [optional](../installation.md#optional-licensed-tools) |
| Run time | 2–3 min |

## 3a · Rebuilding unresolved loops

In the selected chain of each PDB entry, residues the depositors list as unobserved
(`REMARK 465`, or `_pdbx_unobs_or_zero_occ_residues` in mmCIF) are rebuilt with MODELLER when
they form an **internal gap of at most 30 residues**. Missing termini and longer gaps stay
open, and a gap containing a non-standard residue without a standard parent is skipped.

Every observed atom is held fixed; only the gap residues are optimised (`AutoModel`, slow
refinement). Ten models are built with a fixed random seed and the one with the lowest DOPE
score is kept. Only the rebuilt residues are spliced into the downloaded file, so deposited
coordinates, B-factors and the metal stay exactly as they were. Rebuilt atoms get occupancy
0.00 and B-factor 0.00, which is how later steps recognise them (step 10 reports no B-factor
for them).

A chain with a break cannot be aligned or simulated as one polypeptide, and a missing loop
would look like a deletion in the alignment. Thirty residues is the limit because a longer
unresolved stretch is usually a disordered linker or domain, for which a built conformation
would be invention.

Without a working MODELLER the step continues: rebuilt files of an earlier run are reused
where they fill exactly the gaps found now, and the other structures go on with their gaps.

## 3b–3g · One chain, protein and metal only

| | What | Why |
| --- | --- | --- |
| 3b | mmCIF is converted to PDB-format lines (author chain and residue numbers; a multi-character chain ID is cut to its first character) | The tools downstream read PDB format |
| 3c | The first model of an NMR ensemble | One conformation per protein |
| 3d | The chain selected in step 2 | The alignment needs a single polypeptide |
| 3e | The first alternate conformer listed, chosen per residue, not per atom | Each residue stays one consistent conformer; the first is not always `A` |
| 3f | A modified residue of the chain becomes an `ATOM` record of its standard parent from the RCSB Chemical Component Dictionary (HIC → HIS, the N-methylated N-His of fungal LPMOs); which residues were modified is saved in `modified_residues.json` | Alignment, DSSP and force fields know only standard residues; the modification itself is a property of the deposited protein, so it is kept |
| 3g | Waters, ligands and ions are removed, except the active-site metal; `ANISOU` records too | See below |

**The active-site metal** is a Cu, Zn, Ni, Co, Mn, Fe, Mo or Mg `HETATM` with at least two
histidine nitrogens within 3.0 Å. A surface ion bound by one histidine does not qualify. If
several qualify, the one with the most, then the closest, histidine nitrogens is kept.

**No waters or ligands** take part in any calculation. They differ from one crystal to the
next and do not exist in models, so keeping them would make solvent accessibility, charges
and contacts depend on the crystallisation conditions instead of the protein.

## 3h · N-terminal trim to the N-His

The mature LPMO starts with the N-His. Residues before it (signal peptide, tags) are removed,
at the first of these that applies:

1. **Curated**: `curation/n_terminus.csv` gives the first residue, or `keep`.
2. **Annotation**, if the annotated residue is a histidine: the UniProt chain start after a
   signal peptide, then the start of the Pfam catalytic-domain match. The chain start is
   tried first because the Pfam match can begin a few residues later, on another histidine.
3. **SignalP 6**, for models without any annotated start: the residue after the predicted
   signal peptide, if it is a histidine.
4. **First histidine** among the first 35 residues. This also handles a crystal structure
   whose numbering differs from UniProt's.

The annotation is that of the model's own accession (AlphaFold DB), or of the first accession
the entry lists (PDB). SignalP runs in fast mode with the organism group `other` for bacteria
and archaea and `eukarya` for everything else. Viral proteins count as `eukarya`, since their
hosts secrete them.

Rule 4 alone is not reliable for models: a signal peptide can itself contain a histidine. A
model that ends up trimmed by rule 4, or not starting with a histidine, gets a `WARNING` in
the log and should be checked and, if needed, pinned in `n_terminus.csv`. Without SignalP
installed, models without annotation fall to rule 4.

## 3i · C-terminal expression tags

A crystal structure is of a construct; the catalogue describes the protein. A trailing
purification or epitope tag is removed, at the first of these that applies:

1. **Curated**: `curation/c_terminus.csv` gives the last residue, or `keep`.
2. **The entry's annotation**: the trailing run of residues that `SEQADV` records (PDB
   format) or `_struct_ref_seq_dif` rows (mmCIF) call an expression tag, cloning artifact or
   linker.
3. **A known tag** at the C-terminus (`KNOWN_TAGS`).

Only a run reaching the C-terminal end is removed. Annotated tag residues at the N-terminus
are ignored here: that end is 3h's task, and some entries annotate the N-His itself as a tag.

For rule 3 the sequence must *end* with a *prefix* of the tag motif, because the far end of a
tag is usually disordered and so missing from the structure. The prefix must be at least 3
residues for a poly-His tag and 4 for any other tag. A cloning scar directly before it
(`TAG_LINKERS`) goes too. Requiring a prefix at the very end, instead of searching for the
motif anywhere, avoids cutting natural sequence that happens to resemble a tag.

Two limits protect against a wrong cut: at most 30 residues are removed, and at least 100
must remain. Models go through the same rules, since a database sequence can carry a tag as
well; in practice only rules 1 and 3 can apply to them.

## Parameters

| Constant | Value | Meaning |
| --- | --- | --- |
| `MAX_GAP_LENGTH` (`fill_gaps.py`) | `30` residues | Longest internal gap rebuilt |
| `N_MODELS` | `10` | MODELLER models per structure |
| `RAND_SEED` | `-8123` | MODELLER seed, fixed so reruns give the same loops |
| `MODELLER_ENV` | `modeller` | Conda environment searched when `MODELLER_PYTHON` is unset |
| `BRACE_METALS` (`preprocess_structures.py`) | Cu, Zn, Ni, Co, Mn, Fe, Mo, Mg | Elements accepted as the active-site metal |
| `active_site_metal(max_dist)` | `3.0` Å | Metal–histidine nitrogen distance |
| `HIS_SEARCH_LIMIT` | `35` residues | Where the first histidine is looked for (3h rule 4) |
| `TAG_DETAILS` | expression tag, cloning artifact, linker | Annotation texts that mark tag residues |
| `KNOWN_TAGS` | His-tag `HHHHHHHHHH`, Strep-tag II `WSHPQFEK`, myc `(LE)EQKLISEEDL`, FLAG `DYKDDDDK`, HA `YPYDVPDYA`, V5 `GKPIPNPLLGLDST`, S-tag `KETAAAKFERQHMDS`, AviTag `GLNDIFEAQKIEWHE`, thrombin `LVPRGS`, TEV `ENLYFQG`, 3C `LEVLFQGP`, enterokinase `DDDDK`, Factor Xa `IEGR` | Tag motifs for rule 3 |
| `MIN_TAG_MATCH` / `MIN_TAG_MATCH_DEFAULT` | `3` (His-tag) / `4` | Shortest matching prefix |
| `TAG_LINKERS` | LE, GS, GG, SA, AAA, GSG, GGS, SG, GSA | Cloning scars removed with a tag |
| `MAX_TAG_RESIDUES` | `30` | Most residues removed as a tag |
| `MIN_KEPT_RESIDUES` | `100` | Fewest residues left after tag removal |

Environment variables: `MODELLER_PYTHON` (the Python of the MODELLER environment) and
`SIGNALP6` (the `signalp6` executable, default: the one on `PATH`).

## Outputs

### preprocessed/&lt;id&gt;.pdb

--8<-- "_formats/preprocessed_pdb.md"

### filled/modelled_residues.json

--8<-- "_formats/modelled_residues.md"

`filled/<id>.pdb` is the downloaded file with the rebuilt residues inserted; MODELLER's
working files are in `filled/modeller/<id>/`.

### signal_peptides.json

--8<-- "_formats/signal_peptides.md"

SignalP's own output is kept in `signalp/eukarya/` and `signalp/other/`.

### modified_residues.json

--8<-- "_formats/modified_residues.md"

## What the log tells you

For 3a every gap with its decision and, per rebuilt structure, the DOPE score and geometry
checks. For 3b–3i one line per change per structure: dropped chains, renamed residues,
removed ligands, the metal kept, where the N-terminus was cut and by which rule, and any tag
removed with its sequence. `WARNING` lines mark N-termini to check and tags not removed
because of the limits.

## Known limits

- The step fails if any selected structure cannot be prepared. Files of structures no longer
  selected are deleted from `preprocessed/`.
- Tag detection by motif (rule 3) knows only the tags listed; an unlisted tag in an entry
  without annotation stays.
