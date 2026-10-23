# Step 6 · Catalytic domain

Step 6 finds, in every aligned structure, its disulfide bonds, the last residue of its
catalytic domain, and the two histidines of the histidine brace with the copper position
they define. Every later step uses the domain end or the brace.

| Name for `--steps` | `catalytic_domain` or `6` |
| --- | --- |
| Code | `pipeline.py` `step_catalytic_domain` |
| Reads | [Step 2](02-data.md): `pdb_selection.csv`, `domain_annotations.json` · [step 4](04-alignment.md): alignment, aligned structures · [step 5](05-dssp-sasa.md): `secondary_structure.json` |
| Writes | `06_catalytic_domain/disulfides.json`, `catalytic_domain_boundaries.json`, `histidine_brace.json` |
| Run time | under 1 s |

Structures are **not cut** at the domain end. They keep their CBMs and linkers, and the end
is a stored number that later steps and the web app apply, so one alignment serves both the
catalytic-domain and the full-length view.

For the domain end the order of the sub-steps matters: 6b assigns it, 6c and 6d may shorten
it, and 6e, applied last, may lengthen it again and overrides them.

## 6a · Disulfide bonds

Two cysteines are bonded when their SG atoms are at most **3.0 Å** apart.

An intact S–S bond measures 2.03–2.05 Å, but radiation damage stretches bonds in crystal
structures, up to 2.72 Å in the catalogue (4YN2, 7ZJB). No other cysteine pair comes closer
than 3.5 Å, so 3.0 Å separates the two groups cleanly.

The bonds are found from the coordinates because that works for every structure: models have
no `SSBOND` records, some experimental files have none either, and `SSBOND` also lists bonds
to symmetry-related molecules in the crystal, which are not bonds within the chain.

## 6b · Catalytic domain end

A domain annotation applies directly only to AlphaFold DB models, whose residue numbers are
UniProt positions. The structural alignment carries those ends to all other structures.

1. **Votes.** Every AlphaFold DB model whose own UniProt entry has a Pfam catalytic-domain
   match (step 2g) votes with the alignment column of the match's last residue.
2. **Consensus.** Per family, the median vote. The median is robust against a few unusually
   long or short annotations.
3. **Assignment**, recorded in `source`:

| `source` | Applies to | End |
| --- | --- | --- |
| `annotation` | The models that voted | Their own annotated end |
| `consensus` | Other structures of a family with votes | The residue at the consensus column, or in the nearest occupied column within 50. If it lies at most **15** residues before the last residue, the last residue |
| `fallback` | Families without votes, or no residue within 50 columns | The last residue |

The 15-residue rule covers crystal structures of an isolated catalytic domain: there the
construct *is* the domain, and its last residues should not be cut off.

## 6c · Capping long tails

An end from 6b can still lie far outside the domain: a consensus column can fall into the
CBM of a multi-domain protein, an annotation can be over-long, and a fallback is the end of
the chain.

Per family a *cap column* is computed: the median end column of its `annotation` ends (of
all its ends if it has none), leaving out columns beyond 1.5 × the family's median end. An
end more than **100 alignment columns** past the cap column moves back to the last residue
at or before `cap column + 100`.

The margin is in columns, not residues. The alignment has thousands of columns, most of them
insertions, so 100 columns allow for C-terminal variation within a family, while a linker or
CBM lies much further out.

## 6d · Ends inside a linker

If the end lies in a coil stretch of at least **60 residues**, it moves back to the last
helix or strand residue before that stretch (coil = `-` in the three-state secondary
structure of step 5b). A coil that long without any helix or strand is a linker or a
disordered tail, not a loop of the β-sandwich. The rule only shortens the domain.

An end changed by 6c or 6d gets `source` `capped`.

## 6e · Never across a disulfide bond

If one cysteine of a bond lies inside the domain and its partner after the end, the end
moves forward to the partner; this repeats while the new end separates another bond.

A disulfide bond ties two parts of the chain into one folded unit, so a domain boundary
between its cysteines is wrong whatever an annotation says. This mainly concerns AA9 models,
whose Pfam match stops before a C-terminal cysteine that is bonded back to the β-sandwich.

A changed end gets `source` `disulfide`.

## 6f · Histidine brace

The copper is held by the N-terminal histidine (amino group and a ring nitrogen) and a
second histidine further along the chain: the **N-His** and the **Internal-His** (`his1` and
`his2` in the files).

**With an active-site metal** (the one step 3g kept): the two histidines with a nitrogen
closest to the metal; the earlier one in the chain is the N-His.

**Without a metal** (apo crystal structures and all models), the brace is found where the
metal-bound structures have it:

- N-His: the first residue if it is a histidine (step 3h trims to it); else the histidine at
  an N-His column of the metal-bound structures.
- Internal-His: the first histidine at the Internal-His columns of the metal-bound
  structures, most frequent column first, then at columns up to **3** away, nearest first.
  The window is needed because the alignment does not place every Internal-His in exactly
  the same column.

For every brace the smallest distance between a nitrogen of one histidine and a nitrogen of
the other is stored. It is about 3 Å in metal-bound structures; a brace above **7.0 Å** is
reported in the log but kept.

### Pseudo-copper { #pseudo-copper }

For every structure with a brace, the centroid of the five nitrogens that hold the copper is
stored as `pseudo_cu`: N, ND1 and NE2 of the N-His, ND1 and NE2 of the Internal-His (at
least two of them must be present). It stands in for the copper where a structure has no
metal, so that models and apo structures can be analysed like metal-bound ones. In
structures with copper it lies close to the real ion.

### Copper position { #copper-position }

Whether the metal or the pseudo-copper is the copper is decided here, once, and stored as
`copper` with `copper_source`; steps [7c](07-core-fold.md#7c-binding-surface),
[8c and 8d](08-second-sphere.md#the-copper-position-for-8c-and-8d) and
[12](12-hole-hopping.md#the-copper) all use it. The metal kept by step 3g (Cu, Zn, Ni, Co,
Mn, Fe, Mo or Mg) is the copper when a ring nitrogen (ND1 or NE2) of each brace histidine
lies within **3.0 Å** of it, i.e. when both rings hold it. The N-terminal amine is not
required: in a crystal structure whose N-His ring is modelled flipped, the metal sits on the
other ring nitrogen and away from the amine, yet the brace still holds it, and the deposited
coordinates are not changed. Otherwise, and in structures without a metal, the copper is the
pseudo-copper.

## Parameters

Constants of `pipeline.py` unless noted.

| Constant | Value | Used in | Meaning |
| --- | --- | --- | --- |
| `DISULFIDE_MAX_DISTANCE` | `3.0` Å | 6a | Largest SG–SG distance of a bond |
| (in `_consensus_domain_ends`) | `15` residues | 6b | A consensus end this close to the C-terminus becomes the last residue |
| (in `_nearest_position`) | `50` columns | 6b | Search range for a residue near the consensus column |
| (in `_cap_domain_ends`) | `1.5` × median | 6c | End columns beyond this are left out of the cap column |
| `CAP_MARGIN` | `100` columns | 6c | How far past the cap column an end may lie |
| `LINKER_MIN_COIL` | `60` residues | 6d | Shortest coil stretch taken as a linker |
| `BRACE_COLUMN_WINDOW` | `3` columns | 6f | Columns searched next to the known Internal-His columns |
| `BRACE_WARN_N_DIST` | `7.0` Å | 6f | Brace nitrogen distance above which a warning is printed |
| `BRACE_NITROGENS` | N-His N, ND1, NE2; Internal-His ND1, NE2 | 6f | Atoms whose centroid is the pseudo-copper |
| `COPPER_MAX_RING_N_DISTANCE` | `3.0` Å | 6f | Ring nitrogen–metal distance, for each brace histidine, for the metal to be the copper |
| `LPMO_PFAM_IDS` (`fetch_structures.py`) | PF03443, PF03067, PF22810 | 2g, 6b | Pfam families of the catalytic domain |

## Outputs

Residues in these files are [residue numbers](../conventions.md#numbering).

### disulfides.json

--8<-- "_formats/disulfides.md"

### catalytic_domain_boundaries.json

--8<-- "_formats/catalytic_domain_boundaries.md"

### histidine_brace.json

--8<-- "_formats/histidine_brace.md"

## What the log tells you

- **6b**: per family the number of votes, the consensus column, the range and the standard
  deviation. A wide range means the annotations disagree; with two votes the median is their
  midpoint.
- **6c–6d**: the cap column per family and every moved end. 6b prints columns 1-based and
  6c 0-based.
- **6e**: every moved end with the bond that moved it.
- **6f**: how many braces came from the metal and from the alignment, the Internal-His
  columns, a `WARNING` for every brace above 7.0 Å and a line for every structure without a
  complete brace, then in how many structures the copper position is the metal or the
  pseudo-copper, with a `NOTE` for every metal the brace rings do not hold. The log says
  `His-1` and `His-2` for the N-His and the Internal-His.

## Known limits

- **Terminal residues placed far away.** US-align sometimes puts the last one or two
  residues of a chain in a remote column. 6c then ends the domain just before them, so some
  single-domain structures lose their final residue from the catalytic domain.
- **Families without a Pfam match** (AA11, AA17) have no votes; their ends start as the last
  residue and are shaped by 6c–6e alone.
- **A wide brace** is a sign of a structure whose active site is not in the copper-bound
  conformation, or of a model that folds it wrongly. The assignment can still be the right
  pair of histidines; values that depend on the brace geometry are then unreliable.
- **No N-terminal histidine, no brace.** A chain deposited without its N-His has none, and
  steps that need a copper position skip it.
- **The metal is not always copper.** A structure with zinc, nickel or magnesium in the
  brace is treated like a copper one; `metal` names the element and `copper_source` says
  whether it was taken as the copper position.
