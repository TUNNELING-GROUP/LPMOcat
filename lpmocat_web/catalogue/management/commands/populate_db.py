"""
populate_db.py — Django management command that loads the pipeline outputs (data/)
into the LPMO Catalogue DB. It only reads files: everything is calculated by
lpmocat_pipeline/pipeline.py (pipeline step in brackets).

Sources, relative to data/:
  01_proteins/cazy_lpmos.csv                  CAZy entry (1a)
  01_proteins/other_lpmos.csv                 entry of Swiss-Prot, BRENDA or the curated list (1b-1d)
  02_data/pdb_selection.csv                   selected structure, UniProt IDs, EC numbers (2)
  02_data/uniprot_ids.json                    every UniProt accession, source, status (2f)
  02_data/uniprot_references.json             papers cited by the UniProt entries (2f)
  02_data/domain_annotations.json             associated CBMs (2g)
  02_data/taxonomy.json                       NCBI Taxonomy lineage of the organisms (2i)
  02_data/activities.json                     reported substrates, oxidised positions and
                                              the papers reporting them (2j)
  03_preparation/modified_residues.json       modified residues of the deposited chain (3f)
  04_alignment/aligned/<struct_id>.pdb        structure and sequence (4)
  04_alignment/alignment.fasta                aligned sequence (4)
  05_dssp_sasa/sasa.json                      SASA and relative SASA per residue (5a)
  05_dssp_sasa/secondary_structure.json       secondary structure (5b)
  06_catalytic_domain/disulfides.json         disulfide bonds (6a)
  06_catalytic_domain/histidine_brace.json    histidine brace, active-site metal (6f)
  07_core_fold/beta_core.json                 beta-sandwich core strands (7a)
  07_core_fold/loop_regions.json              loop regions (7b), binding surface (7c)
  07_core_fold/binding_surface.json           the binding-surface plane and the copper (7c)
  08_second_sphere/brace_interactions.json    second-sphere aromatic, brace interactions,
                                              second-sphere motif, distal axial residue,
                                              motifs after the brace histidines (8)
  10_flexibility/flexibility_raw.json         flexibility per method (10)
  11_electrostatics/<struct_id>.pqr           atomic charges (11)
  12_hole_hopping/hole_hopping.json           hole hopping paths (12)

A missing file is reported and its fields stay empty. Substrates and regioselectivity are those
the papers report (step 2j, 02_data/activities.json), filled up from the EC numbers where
no paper says anything. Protein codes, clean other names and literature references come from
update_protein_codes, rename_lpmo_to_aa --apply, clean_other_names --apply,
backfill_references and backfill_cazy_references, run afterwards in that order.

Usage:
  python manage.py populate_db [--proteins ID1 ID2] [--force] [--dry-run] [--data-dir DIR]
"""

import csv
import json
import re
from pathlib import Path

from django.core.management.base import BaseCommand

from catalogue.models import (
    Clade, Electrostatics, Flexibility, HoleHopping,
    Organism, Protein, RegioSelecitivity, Substrate, SurfaceAccesability,
)
from catalogue.cbms import uniprot_cbms
from catalogue.overview import species_name

PROJECT_ROOT = Path(__file__).parents[4]
DATA_DIR = PROJECT_ROOT / "data"

FAMILY_SUBSTRATE = {
    "AA9":  "Cellulose",
    "AA10": "Chitin",
    "AA11": "Chitin",
    "AA13": "Starch",
    "AA14": "Xylan",
    "AA15": "Chitin",
    "AA16": "Cellulose",
    "AA17": "Pectin",
}

SUBSTRATE_KEYWORDS = {
    "chitin":     "Chitin",
    "cellulose":  "Cellulose",
    "xylan":      "Xylan",
    "glucan":     "Glucan",
    "starch":     "Starch",
    "amylase":    "Starch",
    "pectin":     "Pectin",
    "xyloglucan": "Xyloglucan",
    "mannan":     "Mannan",
}

# LPMO EC numbers (IUBMB ExplorEnz, enzyme-database.org): 1.14.99.53 lytic chitin monooxygenase
# oxidises C1 and C4; .54 lytic cellulose monooxygenase (C1-hydroxylating) and .55 lytic starch
# monooxygenase oxidise C1; .56 lytic cellulose monooxygenase (C4-dehydrogenating) oxidises C4.
# 1.14.99.- sets neither.
EC_SUBSTRATE = {"1.14.99.53": "Chitin", "1.14.99.54": "Cellulose", "1.14.99.55": "Starch",
                "1.14.99.56": "Cellulose"}
EC_OXIDISED_CARBONS = {"1.14.99.53": ("C1", "C4"), "1.14.99.54": ("C1",), "1.14.99.55": ("C1",),
                       "1.14.99.56": ("C4",)}
REGIOSELECTIVITY_NAMES = {(): "Unknown", ("C1",): "C1-hydroxylating", ("C4",): "C4-dehydrogenating",
                          ("C1", "C4"): "C1-hydroxylating/C4-dehydrogenating"}

# Flexibility methods of flexibility_raw.json; the first one present sets flex_min/flex_max
FLEXIBILITY_METHODS = ("cabsflex", "gnm", "bfactor", "plddt")


def activities_json(entry, references):
    """Protein.activities from a protein's activities.json entry (step 2j), with the metadata of
    the papers it cites; None when the protein has no activity data."""
    if not entry or not (entry.get("substrates") or entry.get("ec")):
        return None
    used = {ref for obs in entry["observations"] + entry["notes"] for ref in obs["refs"]}
    used |= {ref for sources in entry["ec"].values() for refs in sources.values() for ref in refs}
    return json.dumps({**entry, "references": {r: references[r] for r in sorted(used) if r in references}})


def second_sphere_json(site):
    """Protein.second_sphere_motif from a structure's brace_interactions.json entry (8c-8d):
    its motif with the solvent-facing axial residue under "distal_axial", or None."""
    motif, distal = site.get("second_sphere_motif"), site.get("distal_axial")
    if not motif and not distal:
        return None
    return json.dumps({**(motif or {}), "distal_axial": distal, "copper": site.get("copper")})


def brace_motifs_json(site):
    """Protein.brace_motifs from a structure's brace_interactions.json entry (8e), or None."""
    motifs = site.get("brace_motifs")
    return json.dumps(motifs) if motifs else None


def modified_residues_json(entry):
    """Protein.modified_residues from a structure's modified_residues.json entry (3f): the JSON
    list ('[]' = none), or None when the file has no entry for it."""
    return json.dumps(entry) if entry is not None else None


def binding_plane_json(entry):
    """Protein.binding_plane from a structure's binding_surface.json entry (7c), or None."""
    plane, copper = (entry or {}).get("plane"), (entry or {}).get("copper")
    if not plane or not copper or not copper.get("coord"):
        return None
    return json.dumps({"normal": plane["normal"], "copper": copper["coord"],
                       "depth": copper.get("depth_below_plane")})


def pipeline_files(data_dir):
    """Paths of the pipeline outputs under data_dir."""
    structures = data_dir / "02_data"
    alignment = data_dir / "04_alignment"
    return {
        "cazy_csv": data_dir / "01_proteins" / "cazy_lpmos.csv",
        "other_csv": data_dir / "01_proteins" / "other_lpmos.csv",
        "selection_csv": structures / "pdb_selection.csv",
        "uniprot_ids": structures / "uniprot_ids.json",
        "uniprot_references": structures / "uniprot_references.json",
        "domain_annotations": structures / "domain_annotations.json",
        "taxonomy": structures / "taxonomy.json",
        "activities": structures / "activities.json",
        "modified_residues": data_dir / "03_preparation" / "modified_residues.json",
        "aligned_dir": alignment / "aligned",
        "alignment_fasta": alignment / "alignment.fasta",
        "sasa": data_dir / "05_dssp_sasa" / "sasa.json",
        "secondary_structure": data_dir / "05_dssp_sasa" / "secondary_structure.json",
        "disulfides": data_dir / "06_catalytic_domain" / "disulfides.json",
        "braces": data_dir / "06_catalytic_domain" / "histidine_brace.json",
        "beta_core": data_dir / "07_core_fold" / "beta_core.json",
        "loops": data_dir / "07_core_fold" / "loop_regions.json",
        "binding_plane": data_dir / "07_core_fold" / "binding_surface.json",
        "brace_interactions": data_dir / "08_second_sphere" / "brace_interactions.json",
        "flexibility": data_dir / "10_flexibility" / "flexibility_raw.json",
        "pqr_dir": data_dir / "11_electrostatics",
        "hole_hopping": data_dir / "12_hole_hopping" / "hole_hopping.json",
    }


def taxonomy_fields(entry):
    """Organism.taxid and Organism.lineage (JSON) from a taxonomy.json entry (step 2i): those of
    the species (the lineage down to it), as an Organism is a species; the matched taxon's for
    an entry above species level."""
    if not entry:
        return {"taxid": None, "lineage": None}
    species = entry.get("species")
    if not species:
        return {"taxid": entry["taxid"], "lineage": json.dumps(entry["lineage"])}
    lineage = entry["lineage"]
    end = next(i for i, (rank, _) in enumerate(lineage) if rank == "species") + 1
    return {"taxid": species["taxid"], "lineage": json.dumps(lineage[:end])}


def organism_for(cazy_name, kingdom, taxonomy_entry):
    """The Organism of an organism name as CAZy gives it: one per species (species_name), with
    CAZy's kingdom (Archaea, Bacteria, Eukaryota, Viruses, unclassified) and the NCBI Taxonomy
    fields of the species."""
    organism, _ = Organism.objects.get_or_create(name=species_name(cazy_name)[:200] or "Unknown")
    updates = {"kingdom": kingdom} if kingdom else {}
    taxon = taxonomy_fields(taxonomy_entry)
    if taxon["taxid"] is not None:
        updates.update(taxon)
    changed = [field for field, value in updates.items() if getattr(organism, field) != value]
    if changed:
        for field in changed:
            setattr(organism, field, updates[field])
        organism.save(update_fields=changed)
    return organism


def _parse_protein_name(full_name):
    """Split 'name (aliases)(more aliases)' into (clean_name, other_names)."""
    parts = re.findall(r'\(([^)]+)\)', full_name)
    clean = re.sub(r'\s*\([^)]*\)', '', full_name).strip()
    other = "; ".join(parts) if parts else ""
    return clean, other


def _ec_list(ec_numbers):
    return [ec.strip() for ec in (ec_numbers or "").split(";") if ec.strip()]


def _infer_substrates(protein_name, family, ec_numbers, reported=()):
    """The protein's substrates for the catalogue (browse filters, overview charts and the
    compare table): those a paper reports activity on (step 2j), the substrates of the EC
    numbers and the one the protein name mentions; the family's usual substrate when none of
    them gives anything."""
    substrates = set(reported)
    substrates |= {EC_SUBSTRATE[ec] for ec in _ec_list(ec_numbers) if ec in EC_SUBSTRATE}
    name_lower = protein_name.lower()
    named = next((substrate for keyword, substrate in SUBSTRATE_KEYWORDS.items()
                  if keyword in name_lower), None)
    if named:
        substrates.add(named)
    return sorted(substrates) or [FAMILY_SUBSTRATE.get(family, "Unknown")]


def reported_activity(entry):
    """(substrates, regioselectivity) of a protein's activities.json entry (step 2j): the
    substrates a paper reports activity on, named as the catalogue names them, and the
    regioselectivity of the oxidised positions the papers give. ([], None) without a record."""
    from catalogue.activities import substrate_label

    if not entry:
        return [], None
    substrates = sorted(substrate_label(name) for name, summary in entry["substrates"].items()
                        if summary["active"])
    positions = {position for summary in entry["substrates"].values() for position in summary["oxidation"]}
    positions |= {position for obs in entry["observations"] if obs["substrate"] is None
                  for position in obs["oxidation"]}
    carbons = tuple(sorted(positions & {"C1", "C4"}))
    # no position in the papers -> None, so the EC numbers decide (REGIOSELECTIVITY_NAMES[()]
    # is "Unknown", which would otherwise overwrite them)
    return substrates, (REGIOSELECTIVITY_NAMES[carbons] if carbons else None)


def _regioselectivity(ec_numbers):
    """C1-hydroxylating, C4-dehydrogenating, both or Unknown, from the EC numbers."""
    carbons = {carbon for ec in _ec_list(ec_numbers) for carbon in EC_OXIDISED_CARBONS.get(ec, ())}
    return REGIOSELECTIVITY_NAMES[tuple(sorted(carbons))]


def _parse_sequence_from_pdb(pdb_path):
    """Extract amino acid sequence from PDB ATOM records (one letter per residue)."""
    THREE_TO_ONE = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        "SEC": "U", "PYL": "O", "MSE": "M", "HIC": "H",
    }
    seen = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            try:
                res_name = line[17:20].strip()
                res_seq = int(line[22:26].strip())
            except (ValueError, IndexError):
                continue
            if res_seq not in seen:
                seen[res_seq] = THREE_TO_ONE.get(res_name, "X")
    return "".join(seen[k] for k in sorted(seen))


def _sasa_csv(entry):
    """SASA as stored in SurfaceAccesability.sasa_file (read by catalogue.surface):
    residue_num,chain,residue_name,sasa,rsa -- A^2 and relative SASA in % (empty if unknown)."""
    if not entry:
        return "# no SASA: run pipeline step 5 (dssp_sasa)"
    rsa = entry.get("rsa") or [None] * len(entry["sasa"])
    rows = ["residue_num,chain,residue_name,sasa,rsa"]
    rows += [f"{n},{c},{r},{s:.4f},{'' if q is None else f'{q:.2f}'}" for n, c, r, s, q in
             zip(entry["residues"], entry["chains"], entry["names"], entry["sasa"], rsa)]
    return "\n".join(rows)


class Command(BaseCommand):
    help = "Load the pipeline outputs (data/) into the LPMO Catalogue DB"

    def add_arguments(self, parser):
        parser.add_argument("--proteins", nargs="+", metavar="ID",
                            help="Only process these struct_ids (default: all)")
        parser.add_argument("--force", action="store_true",
                            help="Re-load proteins that already exist in the DB")
        parser.add_argument("--dry-run", action="store_true",
                            help="Print what would be done without writing to DB")
        parser.add_argument("--data-dir", default=str(DATA_DIR),
                            help="Pipeline data directory (default: data/)")

    def handle(self, *args, **options):
        self.dry_run = options["dry_run"]
        self.force = options["force"]
        self.files = pipeline_files(Path(options["data_dir"]))

        self.stdout.write("Loading metadata...")
        meta = self._load_metadata()
        self.stdout.write(f"  {len(meta)} proteins in pdb_selection.csv")

        self.stdout.write("Loading alignment FASTA...")
        aligned_seqs = self._load_alignment()
        self.stdout.write(f"  {len(aligned_seqs)} sequences in alignment.fasta")

        self.stdout.write("Loading per-structure pipeline results...")
        self.data = {key: self._load_json(key) for key in (
            "uniprot_ids", "uniprot_references", "disulfides", "braces", "flexibility", "sasa",
            "secondary_structure", "brace_interactions", "loops", "binding_plane",
            "hole_hopping", "modified_residues")}
        self.data["beta_core"] = self._load_json("beta_core", part="structures")
        activities = self._load_json("activities", entries="activity records") or {}
        self.data["activities"] = activities.get("proteins", {})
        self.data["activity_references"] = activities.get("references", {})
        self.data["taxonomy"] = self._load_json("taxonomy", entries="organisms")

        target_ids = options["proteins"] or sorted(meta.keys())
        self.stdout.write(f"\nProcessing {len(target_ids)} proteins...\n")

        n_ok = n_skip = n_err = 0
        for i, sid in enumerate(target_ids, 1):
            if sid not in meta:
                self.stdout.write(self.style.WARNING(f"[{i}] {sid}: not in pdb_selection.csv, skipping"))
                n_skip += 1
                continue

            pdb_path = self.files["aligned_dir"] / f"{sid}.pdb"
            if not pdb_path.exists():
                self.stdout.write(self.style.WARNING(f"[{i}] {sid}: no aligned PDB, skipping"))
                n_skip += 1
                continue

            # Check both protein_code=sid (struct_id) AND selected_structure=sid so that
            # proteins already renamed by update_protein_codes are not re-added as duplicates.
            if not self.force and (
                Protein.objects.filter(protein_code=sid).exists() or
                Protein.objects.filter(selected_structure=sid).exists()
            ):
                self.stdout.write(f"[{i}] {sid}: already in DB (protein_code or selected_structure), skipping")
                n_skip += 1
                continue

            try:
                self._process_protein(i, sid, pdb_path, meta[sid], aligned_seqs.get(sid, ""))
                n_ok += 1
            except Exception as e:
                import traceback
                self.stdout.write(self.style.ERROR(f"[{i}] {sid}: ERROR — {e}"))
                traceback.print_exc()
                n_err += 1

        self.stdout.write(f"\nDone: {n_ok} added, {n_skip} skipped, {n_err} errors")

        if self.force and not options["proteins"]:
            self._remove_dropped_proteins(set(meta))

    def _remove_dropped_proteins(self, selected):
        """Proteins whose structure is no longer in pdb_selection.csv.

        A --force reload without --proteins loads the whole catalogue, so the DB must end up
        holding exactly the current selection. When step 2 picks another entry for a protein
        (2026-09-23: 4A02 -> 4ALS, 5LSV -> 4OPB, 8GUM -> 8GUL, 4OW5 -> 4X27, AF-D8Q364 ->
        AF3-ScAA9A) the new structure is loaded as a new row, while the old one matches
        neither protein_code nor selected_structure of the new struct_id and used to stay
        behind: the catalogue then had both, and the new row could not even take its code
        (VhAA10A was taken, so its reload became VhLPMO10A). The linked records go too --
        the OneToOne key sits on Protein, so nothing cascades."""
        dropped = Protein.objects.exclude(selected_structure__in=selected)
        if not dropped.exists():
            return
        for protein in dropped:
            self.stdout.write(self.style.WARNING(
                f"  removing {protein.protein_code} ({protein.selected_structure}): "
                f"no longer in pdb_selection.csv"))
            if self.dry_run:
                continue
            Flexibility.objects.filter(id__in=protein.flexibility.values_list("id", flat=True)).delete()
            for related, model in (("hole_hopping_id", HoleHopping),
                                   ("surface_accesability_id", SurfaceAccesability),
                                   ("electrostatics_id", Electrostatics)):
                pk = getattr(protein, related)
                if pk:
                    model.objects.filter(id=pk).delete()
            protein.delete()

    def _load_json(self, key, part=None, entries="structures"):
        """A per-structure JSON output of the pipeline (its `part` entry, if given); {} (with a
        warning) if missing."""
        path = self.files[key]
        if not path.exists():
            self.stdout.write(self.style.WARNING(f"  {path} not found - these fields stay empty"))
            return {}
        with open(path) as f:
            data = json.load(f)
        if part:
            data = data.get(part, {})
        self.stdout.write(f"  {len(data)} {entries} in {path.parent.name}/{path.name}")
        return data

    def _load_metadata(self):
        """Returns dict: struct_id → merged metadata of the selection and its step 1 entry
        ("cazy": the row of cazy_lpmos.csv, or of other_lpmos.csv for an entry of Swiss-Prot,
        BRENDA or the curated list)."""
        cazy_by_entry = {}
        cazy_by_uniprot = {}
        cazy_by_genbank = {}
        for listed in ("cazy_csv", "other_csv"):
            if not self.files[listed].exists():
                continue
            with open(self.files[listed]) as f:
                for row in csv.DictReader(f):
                    cazy_by_entry[(row["family"], row["protein_name"], row["organism"])] = row
                    for uid in re.split(r"[;,\s]+", row.get("uniprot_ids", "")):
                        uid = uid.strip()
                        if uid:
                            cazy_by_uniprot.setdefault(uid, row)
                    for gb in (row.get("genbank_ids") or "").split(";"):
                        if gb.strip():
                            cazy_by_genbank.setdefault(gb.strip(), row)

        # Domain annotations of the UniProt entries, for their CBMs
        domain_data = {}
        if self.files["domain_annotations"].exists():
            with open(self.files["domain_annotations"]) as f:
                domain_data = json.load(f)

        meta = {}
        with open(self.files["selection_csv"]) as f:
            for row in csv.DictReader(f):
                sid = row["selected_pdb"].strip()
                if not sid:
                    continue
                # Enrich with CAZy data and CBMs via uniprot match
                cazy = {}
                uids = [uid.strip() for uid in re.split(r"[;,\s]+", row.get("uniprot_ids", "")) if uid.strip()]
                for uid in uids:
                    if uid in cazy_by_uniprot:
                        cazy = cazy_by_uniprot[uid]
                # the CBM families annotated in any of the protein's UniProt entries
                cbms = uniprot_cbms(uids, domain_data)
                # A selection row copies family, name and organism from its own CAZy row, which
                # identifies that row also without UniProt IDs (7VFC, 9HDG); structures found
                # through UniParc or CAZyme3D also record its GenBank ID. A UniProt ID can
                # belong to another CAZy entry with the same sequence.
                cazy = (cazy_by_entry.get((row.get("family"), row.get("protein_name"),
                                           row.get("organism")))
                        or cazy_by_genbank.get((row.get("genbank_id") or "").strip()) or cazy)
                meta[sid] = {**row, "cazy": cazy, "cbms": cbms}
        return meta

    def _load_alignment(self):
        """Returns dict: struct_id → aligned sequence string."""
        result = {}
        try:
            from Bio import SeqIO
            for record in SeqIO.parse(str(self.files["alignment_fasta"]), "fasta"):
                # IDs like "7OKR.pdb:A" → "7OKR"
                sid = record.id.split(".pdb")[0]
                result[sid] = str(record.seq)
        except Exception as e:
            self.stdout.write(self.style.WARNING(f"  Could not load alignment: {e}"))
        return result

    def _process_protein(self, idx, struct_id, pdb_path, meta, aligned_seq):
        row = meta
        family = row.get("family", "").strip()
        protein_name_full = row.get("protein_name", "").strip()
        protein_name, other_names = _parse_protein_name(protein_name_full)
        organism_name = row.get("organism", "Unknown").strip()
        # Active UniProt accessions (step 2f) are shown and linked; a selection file of an older
        # pipeline run has only CAZy's accessions
        uniprot_ids = ((row.get("uniprot_active") if "uniprot_active" in row
                        else row.get("uniprot_ids")) or "").strip()
        identity = self.data["uniprot_ids"].get(struct_id)
        uniprot_all = (json.dumps({acc: {key: entry.get(key) for key in ("source", "status", "organism")}
                                   for acc, entry in identity["accessions"].items()})
                       if identity else None)
        uniprot_references = self.data["uniprot_references"].get(struct_id)
        all_pdb_codes = row.get("all_pdb_codes", "").strip()
        selected_pdb = struct_id

        # EC numbers of pdb_selection.csv (step 2h: CAZy's, else reviewed UniProt's, else the
        # family default); a selection file of an older pipeline run has none, so then those
        # of the CAZy row
        ec_numbers_str = ((row.get("ec_numbers") or "") if "ec_numbers" in row
                          else row.get("cazy", {}).get("ec_numbers", "")).strip()
        # Step 2j: what the papers report; the EC numbers fill the gaps (_infer_substrates,
        # _regioselectivity), so every page uses the same substrates and regioselectivity
        activity_entry = self.data["activities"].get(struct_id)
        reported_substrates, reported_regio = reported_activity(activity_entry)
        substrates = _infer_substrates(protein_name, family, ec_numbers_str, reported_substrates)
        regioselectivity = reported_regio or _regioselectivity(ec_numbers_str)
        cbms_list = row.get("cbms", [])
        associated_cbms_str = ", ".join(cbms_list) if cbms_list else None

        self.stdout.write(f"[{idx}] {struct_id} ({family}, {organism_name})")

        if self.dry_run:
            self.stdout.write("  [DRY RUN] would create protein record")
            return

        # --- Clade ---
        clade, _ = Clade.objects.get_or_create(name=family or "Unknown")

        # --- Organism: the species; the NCBI Taxonomy entry (2i) is keyed by CAZy's name ---
        organism = organism_for(organism_name, (row.get("cazy", {}).get("kingdom") or "").strip(),
                                self.data["taxonomy"].get(organism_name))

        # --- SASA (step 5a) ---
        sasa_text = _sasa_csv(self.data["sasa"].get(struct_id))
        surf, _ = SurfaceAccesability.objects.get_or_create(
            prot_code=struct_id,
            defaults={"sasa_file": sasa_text},
        )
        if self.force:
            surf.sasa_file = sasa_text
            surf.save()

        # --- Electrostatics (step 11) ---
        pqr_path = self.files["pqr_dir"] / f"{struct_id}.pqr"
        pqr_text = ""
        if pqr_path.exists():
            pqr_text = pqr_path.read_text()
        else:
            self.stdout.write(self.style.WARNING(f"  no PQR file for {struct_id}"))
        electro, _ = Electrostatics.objects.get_or_create(
            prot_code=struct_id,
            defaults={"pqr_file": pqr_text},
        )
        if self.force:
            electro.pqr_file = pqr_text
            electro.save()

        # --- Sequence and PDB content (step 4) ---
        sequence = _parse_sequence_from_pdb(pdb_path)
        pdb_content = pdb_path.read_text()

        # --- Disulfide bonds (6a) ---
        bonds = self.data["disulfides"].get(struct_id)
        disulfides = json.dumps(bonds) if bonds is not None else None

        # --- Histidine brace and metal (6f), second-sphere aromatic and interactions (8) ---
        brace = self.data["braces"].get(struct_id, {})
        brace_str = (f"({brace['his1']}, {brace['his2']})"
                     if brace.get("his1") is not None and brace.get("his2") is not None else "")
        site = self.data["brace_interactions"].get(struct_id, {})
        aromatic = site.get("second_sphere_aromatic")
        aromatic_str = str((aromatic["name"], aromatic["residue"])) if aromatic else None
        interactions = site.get("interactions")
        interactions_json = (json.dumps(interactions)
                             if interactions and any(interactions.values()) else None)
        motif_json = second_sphere_json(site)
        brace_motifs = brace_motifs_json(site)

        # --- Beta-sandwich core, loop regions and binding surface (step 7) ---
        regions = self.data["loops"].get(struct_id, {})
        loop_regions = (str(regions["loop_regions"])
                        if regions.get("loop_regions") is not None else None)
        binding_surface = (str(regions["binding_surface"])
                           if regions.get("binding_surface") is not None else None)
        core = self.data["beta_core"].get(struct_id)
        beta_core = (json.dumps({name: s["pieces"] for name, s in core.items()})
                     if core is not None else None)

        # Literature references are populated by backfill_references and
        # backfill_cazy_references; a force re-import carries over the references already stored.
        crystal_ref = other_refs_json = characterization_refs = None

        # --- Replace an existing record if force (match by struct_id OR selected_structure) ---
        if self.force:
            existing = (Protein.objects.filter(protein_code=struct_id).first() or
                        Protein.objects.filter(selected_structure=struct_id).first())
            old_flexibility, old_hole_hopping = [], None
            old_surface = old_electrostatics = None
            if existing:
                crystal_ref = existing.crystal_reference
                other_refs_json = existing.other_crystal_references
                characterization_refs = existing.characterization_references
                old_flexibility = list(existing.flexibility.values_list("id", flat=True))
                old_hole_hopping = existing.hole_hopping_id
                old_surface = existing.surface_accesability_id
                old_electrostatics = existing.electrostatics_id
            Protein.objects.filter(protein_code=struct_id).delete()
            Protein.objects.filter(selected_structure=struct_id).delete()
            # Deleting a protein leaves its flexibility and hole hopping records behind, and its
            # SASA and electrostatics records too: they are found above by prot_code = struct_id,
            # but update_protein_codes renames that to the protein code, so a reload after it
            # creates new ones (915 of each had piled up by 2026-09-21)
            Flexibility.objects.filter(id__in=old_flexibility).delete()
            if old_hole_hopping:
                HoleHopping.objects.filter(id=old_hole_hopping).delete()
            if old_surface and old_surface != surf.id:
                SurfaceAccesability.objects.filter(id=old_surface).delete()
            if old_electrostatics and old_electrostatics != electro.id:
                Electrostatics.objects.filter(id=old_electrostatics).delete()

        # --- Create Protein ---
        protein = Protein(
            protein_code=struct_id,
            protein_name=protein_name[:300],
            other_names=other_names[:300] if other_names else "",
            organism=organism,
            cazy_organism=organism_name[:200],
            clade=clade,
            pdb_code=all_pdb_codes or selected_pdb,
            uniprot=uniprot_ids or None,
            uniprot_all=uniprot_all,
            uniparc=(row.get("uniparc_ids") or "").strip() or None,
            uniprot_references=json.dumps(uniprot_references) if uniprot_references else None,
            selected_structure=selected_pdb,
            sequence=sequence,
            histidine_brace=brace_str[:30],
            metal=brace.get("metal"),
            second_sphere_aromatic=aromatic_str[:30] if aromatic_str else None,
            histidine_brace_interactions=interactions_json,
            second_sphere_motif=motif_json,
            brace_motifs=brace_motifs,
            pdb_structure=pdb_content,
            surface_accesability=surf,
            electrostatics=electro,
            ec_numbers=ec_numbers_str,
            aligned_sequence=aligned_seq,
            crystal_reference=crystal_ref,
            other_crystal_references=other_refs_json,
            characterization_references=characterization_refs,
            binding_surface=binding_surface,
            binding_plane=binding_plane_json(self.data["binding_plane"].get(struct_id)),
            loop_regions=loop_regions,
            beta_core=beta_core,
            disulfides=disulfides,
            modified_residues=modified_residues_json(self.data["modified_residues"].get(struct_id)),
            associated_CBMs=associated_cbms_str,
        )
        protein.save()

        # --- Substrates, oxidised positions and their papers (step 2j) ---
        protein.activities = activities_json(self.data["activities"].get(struct_id),
                                             self.data["activity_references"])
        if protein.activities:
            protein.save(update_fields=["activities"])

        # --- Substrates and regioselectivity, from the EC numbers ---
        for name in substrates:
            protein.substrate.add(Substrate.objects.get_or_create(name=name)[0])
        protein.regioselectivity.add(
            RegioSelecitivity.objects.get_or_create(name=regioselectivity)[0])

        # --- Flexibility (step 10), one record per method, with the secondary structure (5b) ---
        flex_entry = self.data["flexibility"].get(struct_id, {})
        residues = flex_entry.get("residues", [])
        ss = self.data["secondary_structure"].get(struct_id)
        sec_struct = json.dumps({"ss": ss}) if ss else None
        for method in FLEXIBILITY_METHODS:
            values = flex_entry.get(method)
            if not values:
                continue
            flex_obj = Flexibility.objects.create(
                prot_code=struct_id,
                method=method,
                rmsf=json.dumps({"residues": residues, "values": values}),
                sec_struct=sec_struct,
            )
            protein.flexibility.add(flex_obj)

        # --- flex_min/flex_max of the default flexibility method ---
        for method in FLEXIBILITY_METHODS:
            values = [v for v in flex_entry.get(method) or [] if v is not None]
            if values:
                protein.flex_min = float(min(values))
                protein.flex_max = float(max(values))
                break

        # --- Sequence length ---
        protein.sequence_length = len(sequence)

        # --- Hole hopping (step 12) ---
        paths = (self.data["hole_hopping"].get(struct_id) or {}).get("paths")
        if paths:
            protein.hole_hopping = HoleHopping.objects.create(prot_code=struct_id,
                                                              paths_json=json.dumps(paths))
        protein.save()

        n_flex = protein.flexibility.count()
        self.stdout.write(self.style.SUCCESS(
            f"  ✓ {struct_id}: {len(sequence)} aa, {n_flex} flex methods, "
            f"brace={brace_str[:20]}, {len(paths) if paths else 0} hole hopping paths"
        ))
