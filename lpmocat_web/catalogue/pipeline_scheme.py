"""The About page's scheme of how the catalogue is built (_pipeline_scheme.html, drawn by
pipeline_scheme.js): columns of boxes from the public resources through the
pipeline to the web app, joined by arrows (a box's `to`), with details shown on hover.

A resource's logo shows static/images/logos/<key>.svg (or .png, .webp) when that file exists,
its name otherwise, so a logo appears once its file is added."""
from django.contrib.staticfiles import finders

LOGO_DIR = 'images/logos'
LOGO_TYPES = ('svg', 'png', 'webp')

# key: (name, link); the key also names the logo file
RESOURCES = {
    'cazy': ('CAZy', 'https://www.cazy.org/'),
    'pdb': ('RCSB PDB', 'https://www.rcsb.org/'),
    'alphafold': ('AlphaFold DB', 'https://alphafold.ebi.ac.uk/'),
    'cazyme3d': ('CAZyme3D', 'https://pro.unl.edu/CAZyme3D/'),
    'uniprot': ('UniProt', 'https://www.uniprot.org/'),
    'interpro': ('InterPro', 'https://www.ebi.ac.uk/interpro/'),
    'ncbi': ('NCBI Taxonomy', 'https://www.ncbi.nlm.nih.gov/taxonomy'),
    'sifts': ('SIFTS', 'https://www.ebi.ac.uk/pdbe/docs/sifts/'),
    'explorenz': ('ExplorEnz', 'https://www.enzyme-database.org/'),
    'brenda': ('BRENDA', 'https://www.brenda-enzymes.org/'),
    'modeller': ('MODELLER', 'https://salilab.org/modeller/'),
    'signalp': ('SignalP 6.0', 'https://services.healthtech.dtu.dk/services/SignalP-6.0/'),
    'usalign': ('US-align', 'https://zhanggroup.org/US-align/'),
    'dssp': ('DSSP', 'https://github.com/PDB-REDO/dssp'),
    'biopython': ('Biopython', 'https://biopython.org/'),
    'pdb2pqr': ('PDB2PQR', 'https://github.com/Electrostatics/pdb2pqr'),
    'prody': ('ProDy', 'http://www.bahargroup.org/prody/'),
    'cabsflex': ('CABS-flex', 'https://bitbucket.org/lcbio/cabsflex'),
    'ehpath': ('EHPath', 'https://github.com/etransfer/EHPath'),
    'molstar': ('Mol*', 'https://molstar.org/'),
    'plotly': ('Plotly', 'https://plotly.com/javascript/'),
    'd3': ('D3.js', 'https://d3js.org/'),
    'django': ('Django', 'https://www.djangoproject.com/'),
}


# Page of the pipeline documentation (docs/, served at settings.DOCS_URL) for each pipeline
# step, and for the import into the database; a box's hover panel links to it
STEP_DOCS = {'1': 'steps/01-proteins/', '2': 'steps/02-data/', '3': 'steps/03-preparation/',
             '4': 'steps/04-alignment/', '5': 'steps/05-dssp-sasa/',
             '6': 'steps/06-catalytic-domain/', '7': 'steps/07-core-fold/',
             '8': 'steps/08-second-sphere/', '9': 'steps/09-distance-matrix/',
             '10': 'steps/10-flexibility/', '11': 'steps/11-electrostatics/',
             '12': 'steps/12-hole-hopping/'}
IMPORT_DOCS = 'django-import/'
HYDROPHOBICITY_DOCS = 'reference/hydrophobicity-scales/'


def _node(id, name, role, detail, to=(), resource=None, tools=(), step=None, own=False, docs=None):
    """A box: `resource` makes it that resource (its logo, a link), `tools` are the software it
    uses, `own` marks LPMOcat's own methods, `to` the boxes its arrows point to, `docs` its page
    of the pipeline documentation (default: the page of its `step`)."""
    return {'id': id, 'name': name, 'role': role, 'detail': detail, 'to': ' '.join(to),
            'resource': resource, 'tools': list(tools), 'step': step, 'own': own,
            'docs': docs or STEP_DOCS.get(step)}


SCHEME = [
    {'title': 'Public resources', 'note': 'what goes in',
     'nodes': [
         _node('lpmos', 'Characterized LPMOs', 'CAZy, UniProtKB/Swiss-Prot, in-house database',
               'Every LPMO that CAZy lists as characterized or with a solved structure, with its '
               'family, organism, EC numbers, GenBank and UniProt accessions and characterization '
               'papers. Added to them, the LPMOs CAZy misses: UniProtKB/Swiss-Prot entries with an '
               'LPMO EC number whose activity is experimentally supported (ECO:0000269) and cited '
               'to a paper; BRENDA proteins with a substrate record citing a paper that names the '
               'protein\'s genus; and entries of an in-house database.',
               to=['fetch'], step='1'),
         _node('pdb', 'RCSB PDB', 'experimental structures',
               'Experimental structures and their publications. A wild-type entry is preferred, '
               'then X-ray structures, then copper in the active site (within 1 Å of the best '
               'resolution), then resolution, completeness and chain length.',
               to=['fetch'], resource='pdb'),
         _node('alphafold', 'AlphaFold DB', 'predicted models',
               'For entries without an experimental structure: the model of the UniProt entry, or '
               'of an entry with the identical sequence (found through UniParc). Its pLDDT is kept '
               'as per-residue confidence.', to=['fetch'], resource='alphafold'),
         _node('cazyme3d', 'CAZyme3D', 'predicted models',
               'The last resort, for entries without an AlphaFold DB model.',
               to=['fetch'], resource='cazyme3d'),
         _node('uniprot', 'UniProt', 'accessions, annotations, papers',
               'Every accession of each protein (active and deleted, entries with the identical '
               'sequence found through UniParc and merged strain variants), EC numbers of reviewed '
               'entries, signal peptides and domains, the experimentally supported activity '
               'annotations of UniProtKB/Swiss-Prot entries, and the papers the entries cite.',
               to=['fetch'], resource='uniprot'),
         _node('brenda', 'BRENDA', 'reported substrates',
               'The substrate tables of the LPMO EC numbers: which substrates were tested, which '
               'gave no activity, which position was oxidised, and the paper for each of them.',
               to=['fetch'], resource='brenda'),
         _node('interpro', 'InterPro', 'domains',
               'Domain annotations of the active UniProt entries: carbohydrate-binding modules '
               'and the catalytic domain.', to=['fetch'], resource='interpro'),
         _node('ncbi', 'NCBI Taxonomy', 'lineages',
               'The lineage of each source organism, to group the proteins by taxon.',
               to=['fetch'], resource='ncbi'),
     ]},
    {'title': 'Structures', 'note': 'pipeline steps 2–4',
     'nodes': [
         _node('fetch', 'Collect data', 'structures and protein data',
               'One structure per protein: experimental (PDB), else AlphaFold DB, else CAZyme3D; '
               'entries of the in-house database come with their model. Entries with identical sequences are kept once, '
               'preferring an experimental structure, and strain variants of one species are merged. '
               'Then, for every protein: its UniProt accessions, domain annotations, EC numbers '
               '(CAZy, else reviewed UniProt entries, else the family\'s most common one), the '
               'lineage of its organism, and the substrates and oxidised positions each source '
               'reports, each with its paper; negative results are kept and disagreements flagged.',
               to=['prepare'], step='2'),
         _node('prepare', 'Prepare', 'clean and trim',
               'One model, chain and alternate location; modified residues become their parent '
               'amino acid (which ones were modified, e.g. a methylated N-His, is recorded); only the protein and its active-site metal are kept (no waters or '
               'ligands). Internal gaps of up to 30 residues are rebuilt. Every structure is trimmed '
               'to start at the N-terminal histidine of the brace (N-His), using the annotated '
               'signal peptide or, for models without one, a SignalP prediction; C-terminal '
               'expression tags are removed.', to=['align'], tools=['modeller', 'signalp'],
               step='3'),
         _node('align', 'Align', 'all structures superposed',
               'All structures aligned together (multiple-structure mode): every residue gets a '
               'column of one alignment, and all structures share one coordinate frame. LPMO '
               'families share only 2–19 % sequence identity, so the alignment is not refined by '
               'sequence; its gaps are real structural differences.',
               to=['annotation'], tools=['usalign'], step='4'),
     ]},
    {'title': 'Annotation', 'note': 'pipeline steps 5–8',
     'group': _node('annotation', 'Found in every structure', '',
                    'The parts of each structure, found on the aligned structures so that the '
                    'same part is found the same way in every protein.', to=['properties']),
     'nodes': [
         _node('surface', 'Secondary structure and surface accessibility', '',
               'The secondary structure (DSSP), and the solvent-accessible surface area (SASA) per residue '
               '(Shrake–Rupley, protein atoms only) and relative to the residue type\'s maximum, '
               '≥ 25 % counting as surface exposed.', tools=['dssp', 'biopython'], step='5'),
         _node('domain', 'Catalytic domain', 'disulfides, histidine brace',
               'The catalytic domain ends where UniProt/InterPro annotate it for AlphaFold models, '
               'and at the family consensus (the median alignment column of those ends) for all '
               'other structures; the end is moved out of long coil segments and never falls '
               'between the two cysteines of a disulfide bond (sulfur atoms within 3.0 Å). CBMs and '
               'linkers after it are left out of the tree and the catalytic-domain values. The '
               'N-His and the internal His come from the metal, or from the brace columns of '
               'metal-bound structures; without a metal, a pseudo-copper is placed at the centroid '
               'of the five brace nitrogens.', step='6', own=True),
         _node('core', 'β-sandwich core and loops', 'β1–β9; loops L2, L3, LS, L8, LC; binding surface',
               'The core strands β1–β9 are the alignment columns where a strand is common across '
               'the families (weighted so that large families do not dominate), checked by how the '
               'strands pair; β2, the hairpin partner of β1, exists only in AA9. The loops L2, L3, '
               'LS, L8 and LC are the segments between them. The binding surface is the solvent-'
               'accessible residues touched by a plane laid on the copper face, with the copper at '
               'most 3 Å below it.', step='7', own=True),
         _node('site', 'Second-sphere residues', 'around the copper',
               'The aromatic residue closest to the copper (axial position), the side chains '
               'forming hydrogen bonds, π-stacking or cation–π contacts with the histidine brace '
               '(3.5 Å), the second-sphere hydrogen-bonding motif (e.g. H-X8-Q-X-Y), the '
               'shielding residue over the copper\'s solvent side (the Ala of AA10s) and the motifs '
               'after the two brace histidines (e.g. H-G-X4-P-X2-R, H-X-G-P).',
               step='8', own=True),
     ]},
    {'title': 'Properties', 'note': 'pipeline steps 9–12',
     'group': _node('properties', 'Computed for every structure', '',
                    'Each property is computed the same way for every protein, on the aligned '
                    'structures and their catalytic domains.', to=['db']),
     'nodes': [
         _node('tree', 'Similarity tree', 'catalytic domains only',
               'BLOSUM62 distances between the catalytic domains in the structural alignment; the '
               'browse page draws the tree with average linkage (UPGMA). It groups proteins by '
               'sequence similarity and is not a phylogeny.',
               tools=['biopython'], step='9'),
         _node('flex', 'Flexibility', 'B-factor / pLDDT, fluctuations',
               'From the structure: B-factors of experimental structures, pLDDT of models. Computed '
               'for every structure: mean-square fluctuations in a Gaussian network model, and the '
               'RMSF of CABS-flex coarse-grained simulations. The measures have different scales.',
               tools=['prody', 'cabsflex'], step='10'),
         _node('charge', 'Charges', 'atomic partial charges',
               'Atomic partial charges and radii with the AMBER force field, as a simple descriptor '
               'of electrostatic properties. The brace histidines keep their copper-binding '
               'nitrogens unprotonated; other residues get standard protonation states (no pKa '
               'prediction). The copper is omitted, as the force field has no parameters for it.', tools=['pdb2pqr'], step='11'),
         _node('holes', 'Hole hopping', 'electron-transfer paths',
               'Residues that may carry electrons between the protein surface and the copper (or '
               'pseudo-copper): exposed Tyr, Trp and His are donors, buried Trp and Tyr (near the '
               'copper or a proton acceptor) are bridging residues. All paths of up to six residues '
               'are ranked by mean residence time, with assumed parameters for the copper, so the '
               'ranking compares paths rather than giving absolute rates.', tools=['ehpath'], step='12'),
         _node('hydro', 'Hydrophobicity', 'values per residue type',
               'Not computed: each residue type has a fixed value on seven published scales '
               '(Kyte–Doolittle by default, Hopp–Woods, Cornette, Eisenberg, Rose, Janin, '
               'Engelman), as tabulated in the CLC Genomics Workbench manual (QIAGEN). The web '
               'app looks them up per residue and averages them over residue sets. Hopp–Woods is '
               'a hydrophilicity scale.', docs=HYDROPHOBICITY_DOCS),
     ]},
    {'title': 'Catalogue', 'note': 'the database',
     'nodes': [
         _node('explorenz', 'ExplorEnz', 'EC definitions',
               'The IUBMB definitions of the four LPMO EC numbers, taken from ExplorEnz as a fixed '
               'table: they give a protein\'s substrates and regioselectivity only where no paper '
               'reports them.',
               to=['db'], resource='explorenz'),
         _node('db', 'LPMOcat database', 'one record per protein',
               'The pipeline results loaded without further calculation. On loading, substrates '
               'and regioselectivity are derived (papers first, EC numbers otherwise), references '
               'are merged and deduplicated, and organisms are grouped by species.',
               to=['browse', 'compare'], tools=['django'], own=True, docs=IMPORT_DOCS),
         _node('sifts', 'SIFTS', 'PDB–UniProt residue mapping',
               'Places the domains annotated in UniProt, including carbohydrate-binding modules, on '
               'the residues of the structures (or by sequence alignment where SIFTS has no mapping).',
               to=['db'], resource='sifts'),
     ]},
    {'title': 'Web app', 'note': 'what you see',
     'nodes': [
         _node('browse', 'Browse', 'overview, tree, filters, search',
               'Charts of what the catalogue contains, the similarity tree, filters and text or '
               'sequence search, to select proteins.', tools=['plotly', 'd3']),
         _node('compare', 'Compare', '3D viewer, plots, analysis',
               'The selected proteins side by side: 3D structures, plots along the structural '
               'alignment, and similarity and conservation analysis.',
               tools=['molstar', 'plotly']),
     ]},
]


def _logo(key):
    """Static path of the logo file for `key`, or None."""
    for ext in LOGO_TYPES:
        path = f'{LOGO_DIR}/{key}.{ext}'
        if finders.find(path):
            return path
    return None


def _resource(key):
    name, url = RESOURCES[key]
    return {'key': key, 'name': name, 'url': url, 'logo': _logo(key)}


def _resolved(node):
    return dict(node, resource=_resource(node['resource']) if node['resource'] else None,
                tools=[_resource(key) for key in node['tools']])


def pipeline_scheme():
    """SCHEME with the resources and tools of each box resolved to {key, name, url, logo}."""
    return [dict(column, group=_resolved(column['group']) if column.get('group') else None,
                 nodes=[_resolved(node) for node in column['nodes']])
            for column in SCHEME]
