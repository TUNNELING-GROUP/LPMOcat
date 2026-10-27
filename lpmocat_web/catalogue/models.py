from django.db import models
from django.db.models import CASCADE


class Substrate(models.Model):
    name = models.CharField(max_length=100)
    def __str__(self):
        return self.name

class Clade(models.Model):
    name = models.CharField(max_length=30)
    def __str__(self):
        return self.name

class Organism(models.Model):
    # One per species: populate_db stores the species part of CAZy's name (overview.species_name)
    name = models.CharField(max_length=200)
    kingdom = models.CharField(max_length=50, blank=True, default='')
    taxid = models.IntegerField(null=True, blank=True, help_text="NCBI Taxonomy ID")
    lineage = models.TextField(null=True, blank=True,
                               help_text='NCBI Taxonomy lineage down to the organism, JSON: [["domain", "Eukaryota"], ...]')
    def __str__(self):
        return self.name

class RegioSelecitivity(models.Model):
    name = models.CharField(max_length=40)
    def __str__(self):
        return self.name

class SurfaceAccesability(models.Model):
    prot_code = models.CharField(max_length=40)
    sasa_file = models.TextField(max_length=10**10)
    def __str__(self):
        return self.prot_code

class Flexibility(models.Model):
    METHOD_CHOICES = [
        ('md_rmsf',  'MD RMSF'),
        ('bfactor',  'B-factor'),
        ('plddt',    'pLDDT'),
        ('gnm',      'GNM'),
        ('cabsflex', 'CABSflex'),
    ]
    prot_code  = models.CharField(max_length=40)
    method     = models.CharField(max_length=20, choices=METHOD_CHOICES, default='md_rmsf')
    rmsf       = models.TextField(max_length=10**10)
    sec_struct = models.TextField(max_length=10**10, blank=True, null=True,
                                  help_text="Simplified DSSP string (H/E/C per residue)")
    def __str__(self):
        return f"{self.prot_code} ({self.method})"

class Electrostatics(models.Model):
    prot_code = models.CharField(max_length=40)
    pqr_file = models.TextField(max_length=10**10)
    def __str__(self):
        return self.prot_code

class HoleHopping(models.Model):
    prot_code = models.CharField(max_length=40)
    paths_json = models.TextField(max_length=10**10, help_text="JSON with ranked hole hopping paths (exact mean residence time)")
    def __str__(self):
        return self.prot_code

class Protein(models.Model):
    protein_code = models.CharField(max_length=20)
    protein_name = models.CharField(max_length=300)
    other_names = models.TextField(max_length=300)
    organism = models.ForeignKey(Organism, on_delete=models.CASCADE)
    cazy_organism = models.CharField(max_length=200, blank=True, default='',
                                     help_text="The source organism as CAZy names it, with strain")
    clade = models.ForeignKey(Clade, on_delete=models.CASCADE)
    pdb_code = models.TextField(max_length=500)
    uniprot = models.TextField(max_length=500, blank=True, null=True, help_text="Active UniProt accessions (shown and linked), semicolon-separated")
    selected_structure = models.CharField(max_length=20)
    sequence = models.TextField(max_length=10000)
    substrate = models.ManyToManyField(Substrate)
    regioselectivity = models.ManyToManyField(RegioSelecitivity)
    binding_surface = models.TextField(blank=True, null=True)
    loop_regions = models.TextField(blank=True, null=True)
    beta_core = models.TextField(blank=True, null=True,
                                 help_text='Core beta strands as JSON {"b1": ["start-end", ...], ...} '
                                           '(pieces, sequence positions)')
    disulfides = models.TextField(blank=True, null=True,
                                  help_text='Disulfide bonds as JSON [[residue, residue, SG-SG '
                                            'distance (Å)], ...] (residue numbers of the structure)')
    modified_residues = models.TextField(
        blank=True, null=True,
        help_text='Modified amino acids of the deposited chain, replaced by their parent in the '
                  'structure (pipeline step 3f), JSON [{"residue": 1, "name": "HIC", "parent": "HIS"}]; '
                  '[] = none (every model)')
    histidine_brace = models.CharField(max_length=30)
    second_sphere_aromatic = models.CharField(max_length=30, blank=True, null=True)
    histidine_brace_interactions = models.TextField(max_length=10**10, blank=True, null=True)
    second_sphere_motif = models.TextField(
        blank=True, null=True,
        help_text='Second-sphere hydrogen-bonding motif and the residue on the solvent-facing axial '
                  'position of the copper (pipeline step 8c-8d), JSON: {"motif": "H-X8-Q-X-Y", '
                  '"short": "H-Q-Y", "residues": [{"role", "name", "residue", ...}], '
                  '"distal_axial": {"name", "residue", ...} or {"name": null} (open)}')
    brace_motifs = models.TextField(
        blank=True, null=True,
        help_text='Sequence motifs after the two brace histidines (pipeline step 8e), JSON: '
                  '{"n_his": {"motif": "H-G-X4-P-X2-R" or null, "sequence": "HGYVESPASR", '
                  '"residues": [{"role", "name", "residue", "offset"}]}, "internal_his": {...}}')
    binding_plane = models.TextField(
        blank=True, null=True,
        help_text='The plane laid on the copper face (pipeline step 7c), in the aligned '
                  'coordinates the viewer shows, JSON: {"normal": [x, y, z] (unit, towards the '
                  'solvent), "copper": [x, y, z], "depth": copper depth below the plane in A}')
    activities = models.TextField(
        blank=True, null=True,
        help_text='Substrates, oxidised positions and the papers reporting them (pipeline step 2j: '
                  'CAZy EC numbers, Swiss-Prot, BRENDA), JSON: {"ec": {ec: {source: [ref]}}, '
                  '"substrates": {substrate: {"result", "forms", "active", "inactive", "oxidation", '
                  '"sources", "check"}}, "observations": [...], "notes": [...], '
                  '"references": {ref: {"author", "year", "title", "pmid", "doi", "journal"}}}')
    metal = models.CharField(max_length=10, blank=True, null=True,
                             help_text="Metal ion in selected structure (e.g. CU, ZN); null if none")
    pdb_structure = models.TextField(max_length=10*10)
    surface_accesability = models.OneToOneField(SurfaceAccesability,on_delete=models.CASCADE)
    flexibility = models.ManyToManyField(Flexibility) #
    electrostatics = models.OneToOneField(Electrostatics, on_delete=models.CASCADE)
    hole_hopping = models.OneToOneField(HoleHopping, on_delete=models.CASCADE, null=True, blank=True)
    ec_numbers = models.CharField(max_length=200, blank=True, default='', help_text="EC numbers from CAZy, semicolon-separated")
    aligned_sequence = models.TextField(max_length=10000, blank=True, null=True)
    associated_CBMs = models.TextField(max_length=500, blank=True, null=True, help_text="Associated Carbohydrate-Binding Modules")
    crystal_reference = models.TextField(blank=True, null=True, help_text="DOI of the selected crystal structure")
    other_crystal_references = models.TextField(blank=True, null=True, help_text="Structure-paper citations, JSON: {pdb: {doi, pubmed, author, n_authors, year, title, selected}}")
    characterization_references = models.TextField(blank=True, null=True, help_text="Characterization-paper citations from CAZy, JSON list: [{pmid, doi, author, n_authors, year, title}]")
    uniprot_all = models.TextField(blank=True, null=True, help_text="Every UniProt accession of the protein, also deleted ones, JSON: {accession: {source, status, organism}}")
    uniparc = models.TextField(blank=True, null=True, help_text="UniParc IDs of the protein's sequences, semicolon-separated")
    uniprot_references = models.TextField(blank=True, null=True, help_text="Papers the protein's UniProt entries cite for function, activity or structure, JSON list: [{pmid, doi, author, n_authors, year, title, journal, scopes, accessions}]")
    # Pre-calculated fields for performance optimization
    flex_min = models.FloatField(null=True, blank=True, help_text="Minimum RMSF value from flexibility data")
    flex_max = models.FloatField(null=True, blank=True, help_text="Maximum RMSF value from flexibility data")
    sequence_length = models.IntegerField(null=True, blank=True, help_text="Sequence length excluding gaps")

    def __str__(self):
        return self.protein_code

    @property
    def uniprot_accessions(self):
        """Every UniProt accession of the protein, active and deleted, '; '-separated -- for
        search and for matching CAZy entries, which may still list deleted accessions."""
        import json
        import re
        try:
            accessions = list(json.loads(self.uniprot_all or '{}'))
        except (ValueError, TypeError):
            accessions = []
        accessions += [u for u in re.split(r"[;,\s]+", self.uniprot or '') if u and u not in accessions]
        return '; '.join(accessions)

class Protein_files(models.Model):
    filename = models.CharField(max_length=100)
    pdb_file = models.FileField(upload_to="pdb_files/")

    def __str__(self):
        return self.filename
    

