# Step 11 · Electrostatics

Step 11 assigns partial charges and radii to every atom, as a simple descriptor of the
charge distribution of each protein.

| Name for `--steps` | `electrostatics` or `11` |
| --- | --- |
| Code | `pipeline.py` `step_electrostatics` |
| Reads | [Step 4](04-alignment.md): `aligned/<id>.pdb` · [step 6](06-catalytic-domain.md): `histidine_brace.json` |
| Writes | `11_electrostatics/<id>.pqr` |
| External tools | pdb2pqr 2.1.1 |
| Run time | about 1 min |

## What it does

```text
pdb2pqr --ff=amber --chain <input>.pdb <id>.pqr
```

pdb2pqr adds hydrogens, optimises the hydrogen-bond network and assigns AMBER charges and
radii. Before it runs, two things are fixed in the input:

- **Brace histidines.** The N-His is renamed `HIE` (proton on NE2) and the Internal-His
  `HID` (proton on ND1). The N-His binds the copper through ND1 and the Internal-His through
  NE2, so those nitrogens cannot carry a proton; left to itself pdb2pqr would choose by
  hydrogen bonding and could protonate a copper-binding nitrogen.
- **Methylated N-His.** The methyl carbon (CZ) of an N-methylated histidine is removed,
  together with the `HIC` header records. Step 3f already renamed the residue to HIS, and
  the force field has no parameters for the methyl group.

All other residues get pdb2pqr's standard protonation states; no pKa prediction is run. The
`--chain` option keeps the chain identifier in the output.

## The copper is left out

The AMBER force field has no parameters for Cu, so pdb2pqr drops the metal. The charges are
therefore those of the protein without its metal. They describe the protein's surface and
loops well; within a few ångströms of the copper they are incomplete.

## Output

### &lt;id&gt;.pqr

--8<-- "_formats/pqr.md"

All PQR files are regenerated on every run, and files of structures no longer aligned are
removed.

## Known limits

- The command uses pdb2pqr 2.x syntax (`--ff=amber`). Version 3 names the options
  differently and writes slightly different output.
- A structure without a complete brace is run with pdb2pqr's own histidine states; the log
  lists it. The step fails if pdb2pqr fails for any structure.
