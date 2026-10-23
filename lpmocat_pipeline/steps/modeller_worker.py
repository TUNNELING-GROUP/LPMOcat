"""
MODELLER side of the `fill_gaps` pipeline step (see fill_gaps.py).

Runs inside the MODELLER conda env. Builds the `target` sequence of a PIR alignment on the
single `template` structure, optimising only the residues listed in spec["select"] -- every
other atom keeps its template coordinates. Scores of all models go to result.json.

Usage:
    python modeller_worker.py spec.json      # spec written by fill_gaps.py
    python modeller_worker.py --version      # also checks the licence key
"""
import json
import os
import sys


def main():
    if sys.argv[1:] == ["--version"]:
        import modeller  # fails here if the licence key is missing or invalid
        print(getattr(modeller, "__version__", "unknown"))
        return

    with open(sys.argv[1]) as f:
        spec = json.load(f)

    from modeller import Environ, Selection, log
    from modeller.automodel import AutoModel, assess, refine

    os.chdir(spec["workdir"])  # AutoModel writes its models to the working directory
    log.minimal()
    env = Environ(rand_seed=spec["seed"])
    env.io.atom_files_directory = ["."]
    select = spec["select"]

    class GapModel(AutoModel):
        def select_atoms(self):
            return Selection(*[self.residues[i] for i in select])

    a = GapModel(env, alnfile=spec["alignment"], knowns="template", sequence="target",
                 assess_methods=(assess.DOPE,))
    a.starting_model, a.ending_model = 1, spec["n_models"]
    a.md_level = refine.slow
    a.make()

    models = [{"name": m["name"], "dope": m.get("DOPE score"), "molpdf": m.get("molpdf"),
               "failure": str(m["failure"]) if m["failure"] else None}
              for m in a.outputs]
    with open("result.json", "w") as f:
        json.dump({"models": models}, f, indent=1)


if __name__ == "__main__":
    main()
