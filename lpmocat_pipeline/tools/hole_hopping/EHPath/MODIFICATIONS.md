# Changes to EHPath

EHPath (Teo, Wang, Smithwick, Migliore, Therien and Beratan, *Proc. Natl. Acad. Sci. U.S.A.*
2019, 116, 15811–15816, doi:10.1073/pnas.1906394116) is distributed under the GNU General
Public License v3 (`LICENSE`), and so is this modified version.

`EHPath.py` is EHPath modified for the LPMO active site and used by step 12 of the LPMOcat
pipeline (`../02_get_EH_results.py`). `original_EHPath.py` is the version we started from: the
copper node in it follows Sayler, Thomas, Rose and Marletta, *Proc. Natl. Acad. Sci. U.S.A.*
2025, 122, e2411229121 (doi:10.1073/pnas.2411229121), whose authors provided their code on
request. Modified in 2026 by the LPMOcat authors; the changes against the published EHPath:

- **The copper as a node**: residue type `CU`, with λ 1.0 eV, E 1.0 V and radius 1 Å; it may be
  a bridge and an acceptor in the hole- and electron-transfer terms. (In `original_EHPath.py`
  the copper is the type `Cu`, with radius 0.73 Å, left out of the donor nodes, with
  electron-transfer terms as donor and as bridge.)
- **Histidine**: node parameters added (λ 1.0 eV, E 1.0 V, radius 1 Å); its backbone and Cβ
  atoms (N, CA, C, O, CB) are left out of the donor nodes, so the node is its ring; in the
  electron-transfer terms the donor may be His, Cys, Met, Trp or Tyr besides a haem.
- **Command-line arguments instead of prompts**: the donor, bridge and acceptor files and the
  donor residue number are read from the command line (`sys.argv[1]`–`[4]`); the run parameters
  are fixed: `cutoff_num` 5, 5 paths, `type_da` "electron", α (`alpha_reorg`) 1.
- **Terms for other systems disabled**: the terms for haem, Fe–S cluster (SF4), glycine-radical
  and nucleotide donors/acceptors are commented out.
- A `__repr__` for the node class, and the program exits after writing its output.

The other files in this folder (`EHPath_multirun.py`, `EHPath_single_*.py`,
`get_input_files.py`, `consolidate_output.py`, `run_ehpath.sh`, the instructions and the example
CSV files) are those of the EHPath distribution, unchanged.
