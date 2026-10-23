# Step 7 · Core fold

Step 7 names the parts of the catalytic domain: the β-strands of the sandwich core that all
LPMOs share (7a), the loops between them (7b), and the residues that form the flat
substrate-binding surface around the copper (7c).

| Name for `--steps` | `core_fold` or `7` |
| --- | --- |
| Code | `beta_core.py` (7a, 7b), `binding_surface.py` (7c), `pipeline.py` `step_core_fold` |
| Reads | [Step 4](04-alignment.md): alignment, aligned structures · [step 5](05-dssp-sasa.md): `dssp.json`, `sasa.json` · [step 6](06-catalytic-domain.md): `catalytic_domain_boundaries.json`, `histidine_brace.json` · [step 2](02-data.md): `pdb_selection.csv` (families) |
| Writes | `07_core_fold/beta_core.json`, `loop_regions.json`, `binding_surface.json` |
| Run time | about 1.5 min |

Everything is restricted to the catalytic domain (up to the end from step 6).

## 7a · β-sandwich core

The aim is one numbering for all LPMOs: β*n* is the structurally equivalent strand in every
protein, and a protein that lacks a strand skips its number.

**Blocks.** For every alignment column, the fraction of a family's structures with a strand
residue (DSSP `E`) there is computed and *averaged over the families*, so that AA9 and AA10
with their many structures do not decide alone. Columns with a value of at least **0.45**
form blocks; a block bridges up to 2 columns below the threshold. Two consecutive blocks are
merged when more structures have them in one continuous strand than in two.

In the catalogue's data the columns of the conserved strands score 0.52 or more and
family-specific strands 0.36 or less, so 0.45 lies in the gap between them. The step expects **8 blocks**
and stops otherwise, because the loop definitions below depend on this numbering.

**Naming.** In each structure a strand (a run of `E`) takes the block it overlaps most. The
blocks are β1, β3, β4, … β9 in sequence order. **β2** is not a block: it is the hairpin
partner of β1, a strand starting at most 4 residues after β1 with an antiparallel bridge to
it. It exists only in AA9, too rare for a block, but it is part of the sheet.

**Pairing check.** The alignment alone sometimes assigns the wrong strand, where a structure
has two strands near a block. The check uses the sheet topology instead:

- The *ladders* of the core are the strand pairs bridged, in one direction, in at least half
  of the structures that have both strands.
- A core strand that forms none of its ladders, although a partner is present, is replaced
  by the unnamed strand between its neighbouring core strands that has the most bridges of
  those ladders. A missing core strand is searched the same way.
- Unnamed strands after a core strand, up to the next one, that form one of its ladders
  become further *pieces* of it. Pieces at most 3 residues apart are merged.

Only runs of DSSP `E` can be core strands; an isolated bridge (`B`) never is. Strands of
family-specific insertions, such as the one inside the AA10 loop L2, stay unnamed.

## 7b · Loops

A loop is the segment between two core strands:

| Loop | From | To |
| --- | --- | --- |
| L2 | β2 | β3 |
| L3 | β4 | β5 |
| LS | β6 | β7 |
| L8 | β8 | β9 |
| LC | β9 | catalytic domain end, if at least 8 residues |

A loop starts after the *last piece* of its first strand. If a bounding strand is missing,
the nearest present core strand before (for the start) or after (for the end) takes its
place; L2 thus starts after β1 in every family without β2.

Defining loops from the core, and not from per-family templates, makes them comparable
across families and follows each structure's own strands.

## 7c · Binding surface

LPMOs bind a flat, crystalline polysaccharide surface. The binding surface is therefore
defined geometrically, as the residues a plane laid onto the copper face touches. Each
structure is treated on its own coordinates; the alignment is not used.

1. **Copper position**: the [copper position](06-catalytic-domain.md#copper-position) of
   step 6f: the metal when the brace holds it, otherwise the pseudo-copper.
2. **Tangent planes**: 20,000 plane orientations spread evenly on a sphere. For each, the
   plane is moved in from outside until it touches the outermost atom of the protein, so a
   plane never cuts the protein.
3. **Allowed planes** keep the copper at most **3.0 Å** below the plane.
4. **Score**: the number of residues the plane touches (outermost heavy atom within **3.0 Å**
   of it), multiplied by a balance term between 0 and 1 that is low when the touching
   residues all lie on one side of the copper and 1 when they surround it. The best
   orientation is refined in two cones around it (10° and 3°, 2,000 orientations each).
5. **Binding surface**: the touching residues with any solvent accessibility (relative
   SASA > 0, or unknown).

The balance term rejects planes that rest on one shoulder of the protein with the copper at
their edge: a substrate lies across the active site, not beside it.

Nothing is selected by residue type or by loop. Chemistry in the search would make the
composition of the result circular, and that the residues turn out to lie in the loops is a
finding, not an assumption.

**Recessed sites.** If no tangent plane comes within 3.0 Å of the copper, the structure is
marked `recessed` and fitted with the limit raised to the shallowest depth it allows plus
1.0 Å. Such a site is a groove or a bowl behind a rim of loops, and a plane describes it
poorly; the flag and the copper depth tell these cases apart.

The result per structure also carries two descriptors: the copper's depth below the plane,
and how evenly the contacts surround it.

```bash
python lpmocat_pipeline/steps/binding_surface.py 8RRY 7PU1 --pymol   # PyMOL scripts to inspect a plane
python lpmocat_pipeline/steps/binding_surface.py --all
```

Standalone options: `--data-dir`, `--out`, `--pymol`, `--slab`, `--max-cu-depth`,
`--balance`, `--allow-above`.

## Parameters

| Constant | Value | Meaning |
| --- | --- | --- |
| `CORE_MIN_FRACTION` (`beta_core.py`) | `0.45` | Family-averaged strand fraction of a core column |
| `BLOCK_GAP` | `2` columns | Gap a block bridges |
| `HAIRPIN_MAX_GAP` | `4` residues | Longest turn between β1 and β2 |
| `LADDER_MIN_FRACTION` | `0.5` | Share of structures in which a strand pair must be bridged to count as a ladder |
| `PIECE_GAP` | `3` residues | Pieces of a strand this close are merged |
| `CORE_STRANDS` | `9` | Strands of the core (8 blocks + β2) |
| `LC_MIN_LENGTH` | `8` residues | Shortest C-terminal segment called LC |
| `SLAB` (`binding_surface.py`) | `3.0` Å | A residue touches the plane within this distance |
| `MAX_CU_DEPTH` | `3.0` Å | Deepest the copper may lie below the plane |
| `DEPTH_MARGIN` | `1.0` Å | Extra depth allowed for a recessed copper |
| `BALANCE` | `1.0` | Weight of the balance term (0 switches it off) |
| `ALLOW_ABOVE` | `0` | Residues allowed to stick through the plane |
| `DIRECTIONS` | `20000` | Plane orientations in the first search |

## Outputs

Ranges in these files are [sequence positions](../conventions.md#numbering): 1 is the first
residue of the prepared structure.

### beta_core.json

--8<-- "_formats/beta_core.md"

### loop_regions.json

--8<-- "_formats/loop_regions.md"

### binding_surface.json

--8<-- "_formats/binding_surface.md"

## What the log tells you

For 7a: the block columns, a table of how many structures of each family have each strand,
the sheet composition, the ladders found, and every assignment the pairing check changed. For
7b: per family and loop, how many structures have it and its length range. For 7c: per
family the residue counts, the copper depth, the spread, the number of recessed sites and in
which regions the residues lie; then the recessed structures and those without a surface.

## Known limits

- The loop definitions assume the eight blocks found in the current data. A changed set of
  structures that gives another number of blocks stops the step, and `LOOPS` in
  `beta_core.py` must then be revised.
- A structure without a brace has no copper position and so no binding surface.
- For recessed sites the residues reported are those on the rim, not the site itself.
