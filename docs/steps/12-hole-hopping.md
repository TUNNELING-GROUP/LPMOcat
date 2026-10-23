# Step 12 · Hole hopping

Step 12 searches for chains of redox-active residues that could carry charge between the
protein surface and the copper, and ranks them with EHPath.

| Name for `--steps` | `hole_hopping` or `12` |
| --- | --- |
| Code | `hole_hopping_ehpath.py`, `lpmocat_pipeline/tools/hole_hopping/01_classify_resids.py`, `03_analyse_paths.py`, `EHPath/EHPath.py` |
| Reads | [Step 4](04-alignment.md): `aligned/<id>.pdb` · [step 5](05-dssp-sasa.md): `sasa.json` · [step 6](06-catalytic-domain.md): `histidine_brace.json`, `catalytic_domain_boundaries.json` |
| Writes | `12_hole_hopping/hole_hopping.json`, `<id>_paths.txt`, `EH_files/`, `EH_results/`, `<id>_pseudocu.pdb` |
| External tools | EHPath, bundled and modified (see [Installation](../installation.md#bundled-tools)) |
| Run time | about 6 min for 160 structures on 19 worker processes; seconds when nothing changed |

## The copper

The acceptor is the [copper position](06-catalytic-domain.md#copper-position) of step 6f.
When that is the structure's metal, the metal is used; when it is the pseudo-copper, a Cu is
inserted there into a copy of the structure as residue `CU 999` (`<id>_pseudocu.pdb`).

## Residue classes

Only Tyr, Trp and His inside the catalytic domain are considered, and the two brace
histidines are left out. Each residue is represented by its aromatic ring.

| Class | Residues |
| --- | --- |
| **Donor** (start of a path) | Tyr, Trp or His with a relative SASA of at least **10 %** (step 5's whole-residue value) |
| **Bridge** (intermediate) | Every buried Trp; a buried Tyr whose phenolic O is within **3.5 Å** of the copper or within **4.0 Å** of a proton acceptor; and the exposed Tyr and Trp donors. His is never a bridge |
| **Acceptor** | The copper |

A buried Tyr needs a proton acceptor because tyrosine oxidation is proton-coupled: without a
base to take the phenolic proton the hop is not feasible. Proton acceptors are the
carboxylate oxygens of Asp and Glu, the ring nitrogens of a histidine outside the brace, and
any backbone carbonyl oxygen.

The 10 % threshold is deliberately lower than the 25 % the catalogue calls "surface
exposed": a residue needs only some contact with the solvent to exchange charge with an
external reductant, and the stricter threshold leaves structures without any donor.

## EHPath

EHPath (Teo et al. 2019) enumerates the paths from one donor through the bridges to the
acceptor and computes each path's mean residence time from Marcus-theory hopping rates. It is
run once per donor, with these fixed settings:

| Setting | Value | Meaning |
| --- | --- | --- |
| `cutoff_num` | `5` | A path has at most 6 nodes |
| `total_paths` | `5` | Paths reported per donor |
| `type_da` | `electron` | Donor/acceptor type |
| `alpha_reorg` | `1` | Scaling of the reorganisation energy |
| Copper node | λ = 1.0 eV, E = 1.0 V, radius 1 Å | Assumed parameters |
| Histidine node | λ = 1.0 eV, E = 1.0 V, radius 1 Å | Assumed parameters |

Because the copper and histidine parameters are assumed, the residence times rank paths
against each other; they are not absolute rates.

## Ranking and filters

The paths of all donors are pooled and ranked by exact mean residence time, shortest first.
Three kinds are discarded:

- a path of only two nodes (a donor directly on the copper);
- a path starting at a histidine whose second residue is not a tyrosine;
- a path starting at a histidine whose second residue is a tyrosine with its phenolic O
  within 3.5 Å of the copper.

## Reuse

A structure is recalculated only when something it depends on changed. Its fingerprint
(`inputs`) is a SHA-256 over the aligned structure file, the two brace histidines, the
catalytic domain end, the donor threshold, the relative SASA values and the four scripts
(`01_classify_resids.py`, `02_get_EH_results.py`, `03_analyse_paths.py`, `EHPath.py`).
Structures with an unchanged fingerprint keep their result. Editing a script therefore
recalculates everything, and step 4 does too, since it rewrites the aligned files.

## Parameters

| Constant | Value | Meaning |
| --- | --- | --- |
| `EXPOSED_SASA_THRESHOLD` (`01_classify_resids.py`) | `10` % | Relative SASA of a donor |
| `TYR_CU_DISTANCE` | `3.5` Å | Tyr OH to copper |
| `PROTON_ACCEPTOR_DISTANCE` | `4.0` Å | Tyr OH to a proton acceptor |
| `IMPORTANT_ATOMS` | Tyr, Trp and His ring atoms | Residues considered and the atoms of each |
| `SASA_THRESHOLD_OVERRIDES` (`hole_hopping_ehpath.py`) | empty | Per-structure donor thresholds |
| `EXCLUDED_STRUCT_IDS` | `3UAM` | Structures not run: its N-His lies away from the Internal-His, so no copper position can be placed |
| EHPath timeout | `300` s per donor | A donor exceeding it is skipped |

## Output

### hole_hopping.json

--8<-- "_formats/hole_hopping.md"

## Known limits

- The paths are geometric candidates from a static structure with assumed copper
  parameters; they are hypotheses for experiments, not demonstrated electron-transfer routes.
- Cys and Met are not considered, although they can take part in hole hopping.
- A structure without a brace, or excluded by hand, has no paths; `note` says why.
