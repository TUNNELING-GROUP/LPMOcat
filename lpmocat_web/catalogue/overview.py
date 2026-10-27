"""The browse page's "Catalogue overview" tab: one record per protein with the values its
charts count (catalogue_overview.js groups them in the browser, so a hovered bar can name its
proteins, highlight them in the dendrogram and select them)."""
import ast
import csv
import functools
import json
import re
from collections import Counter
from pathlib import Path

from Bio.SeqUtils import seq1

from .dendrogram import CLADE_COLORS, FAMILY_STACK_ORDER
from .brace_motifs import OTHER as OTHER_MOTIF, brace_motifs
from .cbms import structure_cbms
from .second_sphere import NOT_DETERMINED, OPEN_POSITION, second_sphere
from .activities import regioselectivity_from_papers, substrates_from_papers
from .protein_mass import residues_kda
from .surface import SURFACE_RSA, read_sasa_file

# Hydrophobicity scales of the compare page, in the file's column order with Kyte-Doolittle
# (GRAVY) first. Hopp-Woods is a hydrophilicity scale: higher = more hydrophilic.
HYDROPHOBICITY_CSV = Path(__file__).resolve().parent / 'data' / 'hydrophobicity_scales.csv'
DEFAULT_SCALE = 'Kyte-Doolittle'
HYDROPHILIC_SCALES = {'Hopp-Woods'}


@functools.lru_cache(maxsize=1)
def hydrophobicity_scales():
    """{scale name: {one-letter residue: value}}, Kyte-Doolittle first."""
    scales = {}
    try:
        with open(HYDROPHOBICITY_CSV, newline='') as f:
            reader = csv.DictReader(f)
            names = [n for n in reader.fieldnames if n not in ('aa', 'aa full', 'aa classes')]
            rows = list(reader)
        scales = {name: {row['aa'].strip().upper(): float(row[name]) for row in rows if row[name].strip()}
                  for name in names}
    except (OSError, ValueError, KeyError, AttributeError):
        scales = {}
    if DEFAULT_SCALE not in scales:
        from Bio.SeqUtils.ProtParamData import kd
        scales[DEFAULT_SCALE] = dict(kd)
    return {DEFAULT_SCALE: scales[DEFAULT_SCALE],
            **{name: values for name, values in scales.items() if name != DEFAULT_SCALE}}


def mean_hydropathy(sequence, scale):
    """Mean scale value of the residues of `sequence` the scale knows (gaps and unknown letters
    skipped), rounded to 0.001; None without any."""
    values = [scale[c] for c in (sequence or '').upper() if c in scale]
    return round(sum(values) / len(values), 3) if values else None


def gravy(sequence):
    """Kyte-Doolittle grand average of hydropathy (GRAVY) of `sequence`."""
    return mean_hydropathy(sequence, hydrophobicity_scales()[DEFAULT_SCALE])


def _pqr_residues(pqr_text, records=('ATOM', 'HETATM')):
    """[(residue number, insertion code), net charge] in file order, from the `records` (ATOM and
    HETATM) of a PQR file (residue number and insertion code from their fixed columns)."""
    charges = {}
    for line in (pqr_text or '').splitlines():
        if not line.startswith(records):
            continue
        try:
            key = (int(line[22:26]), line[26:27].strip())
            charge = float(line.split()[-2])
        except (ValueError, IndexError):
            continue
        charges[key] = charges.get(key, 0.0) + charge
    return list(charges.items())


def residue_profile(pqr_text, sasa_text, catalytic_end):
    """Values from the per-residue files, each as {'all': whole structure, 'cd': up to the
    catalytic domain end (None without an end)}:
        charge                  net charge (e, rounded to 0.1) of all residues (PQR atom charges)
        surface_charge          net charge of the surface-exposed residues
        surface_length          number of surface-exposed residues
        surface_mass            summed mass (kDa) of the surface-exposed residues
        surface_hydrophobicity  {scale: mean value of the surface-exposed residues}
    Surface exposed = relative SASA >= SURFACE_RSA (pipeline step 5a). A PQR residue gets the SASA
    row at its position when both files list the same residue numbers, else the first row with its
    number. The parts of a missing file are None."""
    residues = _pqr_residues(pqr_text)
    sasa = read_sasa_file(sasa_text)
    empty = {'all': None, 'cd': None}

    def scoped(items, number, value):
        return {'all': value(items),
                'cd': (value([item for item in items if number(item) <= catalytic_end])
                       if catalytic_end is not None else None)}

    def net(selected):
        return round(sum(charge for _, charge in selected), 1) + 0.0   # + 0.0: no -0.0

    profile = {'charge': empty, 'surface_charge': empty, 'surface_length': empty,
               'surface_mass': empty, 'surface_hydrophobicity': None}
    if residues:
        profile['charge'] = scoped(residues, lambda r: r[0][0], net)
    if sasa:
        exposed_rows = [row for row in sasa if row[3] >= SURFACE_RSA]   # NaN: never exposed
        profile['surface_length'] = scoped(exposed_rows, lambda row: row[0], len)
        profile['surface_mass'] = scoped(
            exposed_rows, lambda row: row[0],
            lambda rows: _kda(residues_kda(''.join(seq1(row[1]) for row in rows))))
        profile['surface_hydrophobicity'] = {
            name: scoped(exposed_rows, lambda row: row[0],
                         lambda rows, scale=scale: mean_hydropathy(''.join(seq1(row[1]) for row in rows), scale))
            for name, scale in hydrophobicity_scales().items()}
    if residues and sasa:
        if len(sasa) == len(residues) and all(row[0] == key[0] for row, (key, _) in zip(sasa, residues)):
            rsa = [row[3] for row in sasa]
        else:
            by_number = {}
            for row in sasa:
                by_number.setdefault(row[0], row[3])
            rsa = [by_number.get(key[0], float('nan')) for key, _ in residues]
        exposed = [r for r, value in zip(residues, rsa) if value >= SURFACE_RSA]
        profile['surface_charge'] = scoped(exposed, lambda r: r[0][0], net)
    return profile


def _kda(value):
    return round(value, 2) if value is not None else None


_profiles = {}   # residue_profile results of this process, by input fingerprint


def protein_residue_profile(protein, catalytic_end):
    """residue_profile of a Protein (its electrostatics PQR and SASA records); kept in memory
    while the files are unchanged (parsing all PQR files takes ~0.3 s per page load)."""
    pqr = getattr(getattr(protein, 'electrostatics', None), 'pqr_file', None) or ''
    sasa = getattr(getattr(protein, 'surface_accesability', None), 'sasa_file', None) or ''
    key = (protein.pk, catalytic_end, len(pqr), hash(pqr[-2000:]), len(sasa), hash(sasa[-2000:]))
    if key not in _profiles:
        if len(_profiles) > 2000:
            _profiles.clear()
        _profiles[key] = residue_profile(pqr, sasa, catalytic_end)
    return _profiles[key]


def structure_source(struct_id):
    """Where the selected structure comes from: 'PDB', 'AlphaFold DB', 'CAZyme3D' or 'own model'
    (a model made for the catalogue, pipeline step 1e/2)."""
    struct_id = (struct_id or '').upper()
    if struct_id.startswith('AF3-'):
        return 'own model'
    if struct_id.startswith('AF-'):
        return 'AlphaFold DB'
    if struct_id.startswith('C3D-'):
        return 'CAZyme3D'
    return 'PDB'


# Loop regions in sequence order (lpmocat_pipeline/steps/beta_core.LOOPS); others follow alphabetically
LOOP_ORDER = ['L2', 'L3', 'LS', 'L8', 'LC']
_RANGE = re.compile(r'^\s*(-?\d+)\s*-\s*(-?\d+)\s*$')


def _loop_pieces(loop_regions):
    """{loop: [(start, end), ...]} from Protein.loop_regions ("{'L2': ['6-80'], ...}", sequence
    positions); {} when unset or unreadable."""
    try:
        regions = ast.literal_eval(loop_regions) if loop_regions else {}
    except (ValueError, SyntaxError):
        return {}
    loops = {}
    for loop, ranges in (regions.items() if isinstance(regions, dict) else []):
        pieces = [_RANGE.match(str(r)) for r in (ranges if isinstance(ranges, list) else [ranges])]
        pieces = [(int(m.group(1)), int(m.group(2))) for m in pieces if m]
        if pieces:
            loops[loop] = pieces
    return loops


def loop_lengths(loop_regions):
    """{loop: [residues, 'start-end']} from Protein.loop_regions: the length as end - start + 1
    (the browse page's loop filters count the same way), summed over the pieces of a loop."""
    return {loop: [sum(end - start + 1 for start, end in pieces),
                   ', '.join(f'{start}-{end}' for start, end in pieces)]
            for loop, pieces in _loop_pieces(loop_regions).items()}


def loop_profile(loop_regions, sequence, charges, rsa, rmsf, binding_surface=None):
    """The residues of each loop region (Protein.loop_regions, sequence positions) described as
        residues        [residues of the loop, of its surface-exposed part, of its binding surface]
        rmsf            [mean CABS-flex RMSF (Å, rounded to 0.01) of the same three sets]
        charge          [net charge (e, rounded to 0.1), that of its surface-exposed residues,
                         that of its substrate-binding residues]
        hydrophobicity  [[mean value, surface-exposed, substrate-binding], ...] per scale, in
                        the order of hydrophobicity_scales() (lists keep the page's JSON small)
    Surface exposed = relative SASA >= SURFACE_RSA; substrate-binding = the residues the plane of
    pipeline step 7c touches (Protein.binding_surface), which are a subset of the loop.
    charges, rsa and rmsf hold one value per residue in the order of `sequence` (the structure's);
    a list of another length is ignored. A value without residues or data is None."""
    n = len(sequence or '')

    def per_residue(values):
        return values if values is not None and len(values) == n else None

    charges, rsa, rmsf = per_residue(charges), per_residue(rsa), per_residue(rmsf)
    scales = hydrophobicity_scales()

    def net(positions):
        if charges is None or not positions:
            return None
        return round(sum(charges[i] for i in positions), 1) + 0.0   # + 0.0: no -0.0

    def hydropathy(positions, scale):
        return mean_hydropathy(''.join(sequence[i] for i in positions), scale) if positions else None

    binding = _positions(binding_surface, n)

    def mean_rmsf(selected):
        values = [rmsf[i] for i in selected if rmsf is not None and rmsf[i] is not None]
        return round(sum(values) / len(values), 2) if values else None

    profile = {}
    for loop, pieces in _loop_pieces(loop_regions).items():
        positions = sorted({i for start, end in pieces for i in range(max(start, 1) - 1, min(end, n))})
        exposed = [i for i in positions if rsa is not None and rsa[i] >= SURFACE_RSA]   # NaN: never
        touching = [i for i in positions if i in binding]
        profile[loop] = {
            # residues of the loop, of its surface-exposed part and of its binding surface
            'residues': [len(positions),
                         len(exposed) if rsa is not None else None,
                         len(touching) if binding_surface else None],
            'rmsf': [mean_rmsf(positions), mean_rmsf(exposed) if rsa is not None else None,
                     mean_rmsf(touching) if binding_surface else None],
            'charge': [net(positions), net(exposed) if rsa is not None else None, net(touching)],
            'hydrophobicity': [[hydropathy(positions, scale),
                                hydropathy(exposed, scale) if rsa is not None else None,
                                hydropathy(touching, scale)]
                               for scale in scales.values()],
        }
    return profile


def _positions(regions, length):
    """0-based residue indices of a "{region: ['start-end', ...]}" field (sequence positions),
    clipped to a sequence of `length` residues."""
    return {i for pieces in _loop_pieces(regions).values() for start, end in pieces
            for i in range(max(start, 1) - 1, min(end, length))}


def residue_set_profile(loop_regions, binding_surface, sequence, charges, rsa):
    """The two residue sets the distribution charts offer besides all and surface-exposed
    residues, each as {'residues': count, 'mass': kDa (residues_kda), 'charge': e,
    'hydrophobicity': {scale: value}}:
        loops    every loop region taken together (Protein.loop_regions)
        binding  the substrate-binding residues -- what the plane of pipeline step 7c touches
                 (Protein.binding_surface)
    Both lie inside the catalytic domain, so neither changes with the region switch. A value
    without residues or without the underlying file is None."""
    n = len(sequence or '')
    charges = charges if charges is not None and len(charges) == n else None
    scales = hydrophobicity_scales()
    sets = {'loops': (loop_regions, sorted(_positions(loop_regions, n))),
            'binding': (binding_surface, sorted(_positions(binding_surface, n)))}
    profile = {}
    for name, (field, positions) in sets.items():
        profile[name] = {
            # no field = not determined (6RS9 has no N-His, so no copper and no binding plane),
            # which the size chart must leave out rather than draw at 0 residues
            'residues': len(positions) if field else None,
            'mass': _kda(residues_kda(''.join(sequence[i] for i in positions))) if positions else None,
            'charge': (round(sum(charges[i] for i in positions), 1) + 0.0
                       if charges is not None and positions else None),
            'hydrophobicity': {
                label: (mean_hydropathy(''.join(sequence[i] for i in positions), scale)
                        if positions else None)
                for label, scale in scales.items()},
        }
    return profile


_loop_profiles = {}   # loop_profile results of this process, by input fingerprint


def protein_loop_profile(protein):
    """loop_profile of a Protein: its PQR charges (ATOM records), SASA records and CABS-flex RMSF
    (prefetched flexibility); kept in memory while these are unchanged."""
    pqr = getattr(getattr(protein, 'electrostatics', None), 'pqr_file', None) or ''
    sasa = getattr(getattr(protein, 'surface_accesability', None), 'sasa_file', None) or ''
    flex = next((f for f in protein.flexibility.all() if f.method == 'cabsflex'), None)
    rmsf_text = (flex.rmsf if flex else None) or ''
    key = (protein.pk, protein.loop_regions, protein.binding_surface, protein.sequence,
           len(pqr), hash(pqr[-2000:]), len(sasa), hash(sasa[-2000:]), hash(rmsf_text))
    if key not in _loop_profiles:
        if len(_loop_profiles) > 2000:
            _loop_profiles.clear()
        try:
            rmsf = json.loads(rmsf_text).get('values') if rmsf_text else None
        except (TypeError, ValueError, AttributeError):
            rmsf = None
        _loop_profiles[key] = loop_profile(
            protein.loop_regions, protein.sequence,
            [charge for _, charge in _pqr_residues(pqr, records=('ATOM',))] or None,
            [row[3] for row in read_sasa_file(sasa)] or None, rmsf, protein.binding_surface)
    return _loop_profiles[key]


_residue_sets = {}   # residue_set_profile results of this process, by input fingerprint


def protein_residue_sets(protein):
    """residue_set_profile of a Protein (its loops, binding surface and PQR charges); kept in
    memory while these are unchanged."""
    pqr = getattr(getattr(protein, 'electrostatics', None), 'pqr_file', None) or ''
    sasa = getattr(getattr(protein, 'surface_accesability', None), 'sasa_file', None) or ''
    key = (protein.pk, protein.loop_regions, protein.binding_surface, protein.sequence,
           len(pqr), hash(pqr[-2000:]), len(sasa))
    if key not in _residue_sets:
        if len(_residue_sets) > 2000:
            _residue_sets.clear()
        _residue_sets[key] = residue_set_profile(
            protein.loop_regions, protein.binding_surface, protein.sequence,
            [charge for _, charge in _pqr_residues(pqr, records=('ATOM',))] or None,
            [row[3] for row in read_sasa_file(sasa)] or None)
    return _residue_sets[key]


_GLYCOSITE = re.compile(r'(?=N[^P][ST])')
_STANDARD_RESIDUES = set('ACDEFGHIKLMNPQRSTVWY')


def _region_of(position, loops, core):
    """The loop ('L2') or core strand ('β9') holding a 1-based sequence position, else 'other'."""
    for name, pieces in loops.items():
        if any(start <= position <= end for start, end in pieces):
            return name
    for name, pieces in core.items():
        if any(start <= position <= end for start, end in pieces):
            return 'β' + name[1:] if name.startswith('b') else name
    return 'other'


# Aromatic residues on the binding surface (user decision 2026-10-22: His included; the two brace
# histidines are the copper ligands, on every binding surface, so they are left out)
AROMATIC_RESIDUES = 'FWYH'
_AROMATIC_NAMES = {'F': 'Phe', 'W': 'Trp', 'Y': 'Tyr', 'H': 'His'}


def brace_numbers(histidine_brace):
    """The brace histidines' residue numbers from Protein.histidine_brace ('(26, 136)'), or ()."""
    return tuple(int(n) for n in re.findall(r'-?\d+', histidine_brace or ''))


def modified_residues(protein):
    """Protein.modified_residues as a list ([{'residue', 'name', 'parent'}]), or None when not
    recorded."""
    try:
        return json.loads(protein.modified_residues) if protein.modified_residues else None
    except (TypeError, ValueError):
        return None


def n_his_methylation(protein):
    """'yes' when the deposited PDB entry has the N-His as HIC (N-methylated His, made by
    fungal expression hosts), 'no' when it has a plain His; None for a model (a prediction
    shows no modification) or when not recorded. The catalogue's structure has His either way
    (pipeline step 3f)."""
    found = modified_residues(protein)
    if found is None or structure_source(protein.selected_structure) != 'PDB':
        return None
    n_his = (brace_numbers(protein.histidine_brace) or (None,))[0]
    if n_his is None:
        return None
    return 'yes' if any(m['residue'] == n_his and m['name'] == 'HIC' for m in found) else 'no'


def binding_aromatics(sequence, residue_numbers, binding_surface, loop_regions=None, beta_core=None,
                      brace=()):
    """The aromatic residues (Phe, Trp, Tyr, His; not the brace histidines) among the
    substrate-binding residues (pipeline step 7c), in sequence order, each {'label': 'Tyr54',
    'aa': 'Y', 'residue': 54 (residue number), 'position': sequence position, 'region': 'L2' |
    'β9' | 'other'}; None without a binding surface."""
    if not binding_surface or not sequence:
        return None
    loops, core = _loop_pieces(loop_regions), _loop_pieces(beta_core)
    numbers, brace = list(residue_numbers or []), set(brace or ())
    out = []
    for i in sorted(_positions(binding_surface, len(sequence))):
        aa = sequence[i].upper()
        if aa not in AROMATIC_RESIDUES:
            continue
        number = numbers[i] if i < len(numbers) else i + 1
        if aa == 'H' and number in brace:
            continue
        out.append({'label': f'{_AROMATIC_NAMES[aa]}{number}', 'aa': aa, 'residue': number,
                    'position': i + 1, 'region': _region_of(i + 1, loops, core)})
    return out


def protein_binding_aromatics(protein):
    """binding_aromatics of a Protein (its sequence, SASA residue order, regions and brace)."""
    sasa = getattr(getattr(protein, 'surface_accesability', None), 'sasa_file', None) or ''
    return binding_aromatics(protein.sequence, [row[0] for row in read_sasa_file(sasa)],
                             protein.binding_surface, protein.loop_regions, protein.beta_core,
                             brace_numbers(protein.histidine_brace))


def structural_features(sequence, domain, residue_numbers, catalytic_end, disulfides, loop_regions,
                        beta_core, binding_surface, brace_interactions=None, aligned_sequence=None,
                        hole_hopping_paths=None, brace=()):
    """Per-protein values the selection overview's "What separates these groups?" scores, besides
    the measures and loop values the charts already carry. `domain` is the catalytic-domain
    sequence, `residue_numbers` the structure's residue numbers in sequence order (the SASA
    records), `catalytic_end` a residue number. Keys (a value that cannot be computed is None):
        disulfides         bonds with both Cys inside the catalytic domain
        disulfide_pattern  those bonds by the regions they join, 'L2–β9 + other–LC' ('none')
        b2                 'yes' / 'no': the AA9 hairpin strand after β1 (beta_core)
        pi                 isoelectric point of the catalytic domain (Biopython)
        glycosites         N-X-S/T sequons (X not Pro) in the catalytic domain: potential sites only
        binding_aromatic   Phe/Trp/Tyr/His (not the brace histidines, `brace` residue numbers),
                           binding_amide Asn/Gln, binding_charged Asp/Glu/Lys/Arg among the
                           substrate-binding residues (pipeline step 7c)
        binding_aromatic_residues  the aromatic ones themselves (binding_aromatics)
        brace_hbonds       residues whose side chain H-bonds a brace histidine (pipeline step 8b),
        brace_stacking     ... is pi-stacked on one, brace_cation_pi ... makes a cation-pi contact
        hh_positions       the residues of every hole-hopping path (pipeline step 12) as
                           'column:residue' of the structural alignment (1-based column, one-letter
                           code), sorted by column; [] = no path; the metal is not a residue"""
    features = dict.fromkeys(('disulfides', 'disulfide_pattern', 'b2', 'pi', 'glycosites',
                              'binding_aromatic', 'binding_aromatic_residues', 'binding_amide', 'binding_charged',
                              'brace_hbonds', 'brace_stacking', 'brace_cation_pi', 'hh_positions'))
    sequence = sequence or ''
    loops, core = _loop_pieces(loop_regions), _loop_pieces(beta_core)
    position_of = {number: i + 1 for i, number in enumerate(residue_numbers or [])}

    try:
        bonds = json.loads(disulfides) if disulfides else None
    except (TypeError, ValueError):
        bonds = None
    if bonds is not None and catalytic_end is not None:
        inside = [b for b in bonds if len(b) >= 2 and b[0] <= catalytic_end and b[1] <= catalytic_end]
        features['disulfides'] = len(inside)
        pairs = []
        for a, b in (sorted(bond[:2]) for bond in inside):
            if a in position_of and b in position_of:
                pairs.append(f'{_region_of(position_of[a], loops, core)}–{_region_of(position_of[b], loops, core)}')
        features['disulfide_pattern'] = ' + '.join(sorted(pairs)) if pairs else ('none' if not inside else None)
    if core:
        features['b2'] = 'yes' if 'b2' in core else 'no'

    residues = ''.join(c for c in (domain or '').upper() if c in _STANDARD_RESIDUES)
    if residues:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        features['pi'] = round(ProteinAnalysis(residues).isoelectric_point(), 2)
        features['glycosites'] = len(_GLYCOSITE.findall(residues))

    if binding_surface:
        touching = ''.join(sequence[i] for i in sorted(_positions(binding_surface, len(sequence))))
        aromatics = binding_aromatics(sequence, residue_numbers, binding_surface, loop_regions,
                                      beta_core, brace)
        features['binding_aromatic'] = len(aromatics)
        features['binding_aromatic_residues'] = aromatics
        features['binding_amide'] = sum(c in 'NQ' for c in touching)
        features['binding_charged'] = sum(c in 'DEKR' for c in touching)

    try:
        contacts = json.loads(brace_interactions) if brace_interactions else None
    except (TypeError, ValueError):
        contacts = None
    if isinstance(contacts, dict):
        for key, kind in (('brace_hbonds', 'hydrogen_bonds'), ('brace_stacking', 'pi_pi_stacking'),
                          ('brace_cation_pi', 'cation_pi')):
            features[key] = len({c.get('other_res_id') for c in contacts.get(kind) or []})

    try:
        paths = json.loads(hole_hopping_paths) if hole_hopping_paths else None
    except (TypeError, ValueError):
        paths = None
    if isinstance(paths, list) and aligned_sequence and residue_numbers:
        # residue number -> (column, one-letter code); the SASA records are the aligned residues in order
        column_of, k = {}, 0
        for c, aa in enumerate(aligned_sequence):
            if aa != '-':
                if k < len(residue_numbers):
                    column_of[residue_numbers[k]] = (c + 1, aa.upper())
                k += 1
        on_path = set()
        for path in paths:
            for label in path.get('path') or []:
                match = _PATH_RESIDUE.match(str(label))
                if match and int(match.group(2)) in column_of:
                    on_path.add(column_of[int(match.group(2))])
        features['hh_positions'] = [f'{c}:{aa}' for c, aa in sorted(on_path)]
    return features


# a path residue label such as 'TYR91' or 'TRP-3'; metals ('CU300', 'NI1231') do not match
_PATH_RESIDUE = re.compile(r'^(ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HIS|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL)(-?\d+)')


_structural = {}   # structural_features results of this process, by input fingerprint


def protein_structural_features(protein, catalytic_end, domain):
    """structural_features of a Protein (its disulfides, core strands, binding surface, brace
    contacts, hole-hopping paths, aligned sequence and SASA residue order); kept in memory while
    these are unchanged."""
    sasa = getattr(getattr(protein, 'surface_accesability', None), 'sasa_file', None) or ''
    paths = getattr(getattr(protein, 'hole_hopping', None), 'paths_json', None) or ''
    # populate_db stores no contacts as an empty field: with a brace, that is 0 contacts, not unknown
    contacts = protein.histidine_brace_interactions or ('{}' if protein.histidine_brace else None)
    key = (protein.pk, catalytic_end, domain, protein.disulfides, protein.loop_regions, protein.beta_core,
           protein.histidine_brace, protein.modified_residues, protein.selected_structure,
           protein.binding_surface, len(sasa), hash(sasa[-2000:]), contacts, hash(paths),
           hash(protein.aligned_sequence or ''))
    if key not in _structural:
        if len(_structural) > 2000:
            _structural.clear()
        _structural[key] = structural_features(
            protein.sequence, domain, [row[0] for row in read_sasa_file(sasa)], catalytic_end,
            protein.disulfides, protein.loop_regions, protein.beta_core, protein.binding_surface, contacts,
            protein.aligned_sequence, paths, brace_numbers(protein.histidine_brace))
        # from the PDB entry, not the structure: 'yes' / 'no' (n_his_methylation), None for a model
        _structural[key]['n_his_methylated'] = n_his_methylation(protein)
    return _structural[key]


# Taxon groups in chart order: eukaryotes, prokaryotes, viruses, then the unresolved ones
TAXON_GROUPS = ['Fungi', 'Oomycetes', 'Animals', 'Plants', 'Other eukaryotes', 'Eukaryota',
                'Bacteria', 'Archaea', 'Viruses', 'Metagenomes', 'Unclassified', 'Not recorded']
# Without an NCBI lineage: CAZy's kingdom
_KINGDOM_GROUPS = {'unclassified': 'Unclassified', '': 'Not recorded'}
# Lineage ranks named in the tooltip
_LINEAGE_RANKS = ('domain', 'superkingdom', 'acellular root', 'kingdom', 'phylum', 'class',
                  'order', 'family')


def _lineage(organism):
    """The organism's NCBI lineage [[rank, name], ...], or []."""
    try:
        return json.loads(organism.lineage or '[]')
    except (TypeError, ValueError):
        return []


def taxon_group(lineage, kingdom=''):
    """(group, subgroup) of an organism for the taxonomy chart: Fungi, Oomycetes, Animals, Plants
    and other eukaryotes, Bacteria, Archaea, Viruses, Metagenomes; the subgroup is the phylum
    (the class for arthropods, the family for viruses) or None. Without a lineage, CAZy's
    kingdom."""
    ranks = {}
    for rank, name in lineage:
        ranks.setdefault(rank, name)
    names = {name for _, name in lineage}
    if not lineage:
        return _KINGDOM_GROUPS.get(kingdom, kingdom), None
    if 'Viruses' in names:
        return 'Viruses', ranks.get('family')
    if 'metagenomes' in names:
        return 'Metagenomes', None
    kingdom_rank = ranks.get('kingdom')
    if kingdom_rank == 'Fungi':
        return 'Fungi', ranks.get('phylum')
    if kingdom_rank == 'Metazoa':
        if ranks.get('phylum') == 'Arthropoda':
            return 'Animals', ranks.get('class') or 'Arthropoda'
        return 'Animals', ranks.get('phylum')
    if kingdom_rank == 'Viridiplantae':
        return 'Plants', ranks.get('phylum')
    if 'Oomycota' in names:
        return 'Oomycetes', None
    domain = ranks.get('domain') or ranks.get('superkingdom')
    if domain in ('Bacteria', 'Archaea'):
        return domain, ranks.get('phylum')
    if domain == 'Eukaryota':
        return 'Other eukaryotes', ranks.get('phylum')
    return 'Unclassified', None


def species_name(organism_name):
    """The species part of an organism name as CAZy gives it: the first two words, three for an
    unnamed species ('Pestalotiopsis sp. NCi6', a species of its own in NCBI Taxonomy); strain,
    serovar, subspecies ... are dropped. A species name comes back unchanged. CAZy's names, not
    NCBI's: NCBI has renamed some (Lentinus similis is Panus similis there)."""
    words = (organism_name or '').split()
    return ' '.join(words[:3 if len(words) > 2 and words[1] in ('sp.', 'sp') else 2])


def lineage_text(lineage):
    """'Eukaryota › Fungi › Ascomycota › Sordariomycetes …' for tooltips: the main ranks above
    the organism (all names for a lineage without ranks, e.g. metagenomes)."""
    above = lineage[:-1]
    ranked = [name for rank, name in above if rank in _LINEAGE_RANKS]
    names = ranked or [name for _, name in above if name != 'cellular organisms']
    return ' › '.join(names)


def _family_order(name):
    """AA9 before AA10: families sorted by their number."""
    number = re.search(r'\d+', name)
    return (int(number.group()) if number else float('inf'), name)


def _scopes(catalytic, whole):
    return {'cd': catalytic, 'all': whole}


def catalogue_summary(proteins):
    """What the catalogue holds, for the About and How to Use pages: protein, experimental
    structure, predicted model and species counts, the families in number order and the taxon
    groups with their protein counts, largest first."""
    sources = Counter(structure_source(p.selected_structure) for p in proteins)
    groups = Counter(taxon_group(_lineage(p.organism), p.organism.kingdom or '')[0] for p in proteins)
    return {
        'proteins': len(proteins),
        'experimental': sources['PDB'],
        'predicted': len(proteins) - sources['PDB'],
        'alphafold': sources['AlphaFold DB'],
        'cazyme3d': sources['CAZyme3D'],
        'own_models': sources['own model'],
        'species': len({species_name(p.organism.name) for p in proteins}),
        'families': sorted({p.clade.name for p in proteins}, key=_family_order),
        'taxon_groups': groups.most_common(),
    }


# Regioselectivity as the app shows it (browse and compare pages): C1, C4, C1/C4 or Unknown; the
# catalogue's RegioSelecitivity names are the IUBMB's (C1-hydroxylating ...)
REGIOSELECTIVITY_SHORT = {
    'C1-hydroxylating': 'C1',
    'C4-dehydrogenating': 'C4',
    'C1-hydroxylating/C4-dehydrogenating': 'C1/C4',
    'Unknown': 'Unknown',
}


def regioselectivity_label(names):
    """'C1', 'C4', 'C1/C4' or 'Unknown' for a protein's regioselectivity names (several: the
    positions of all of them)."""
    positions = set()
    for name in names:
        short = REGIOSELECTIVITY_SHORT.get(name, name)
        if short != 'Unknown':
            positions.update(short.split('/'))
    return '/'.join(sorted(positions)) or 'Unknown'


def regioselectivity_of(protein):
    """[short name, source] of the protein's regioselectivity: 'papers' when it is the oxidised
    positions the literature reports (pipeline step 2j), 'ec' when it follows from the EC
    numbers, 'none' when neither gives one."""
    label = regioselectivity_label(r.name for r in protein.regioselectivity.all())
    if label == 'Unknown':
        return [label, 'none']
    return [label, 'papers' if regioselectivity_from_papers(protein) else 'ec']


def catalogue_overview(proteins, catalytic_sequence):
    """The panel's data for `proteins` (the browse page's list, with cat_domain_length,
    sequence_length, cat_domain_mass, sequence_mass and residue_profile set).
    `catalytic_sequence(protein)` gives the aligned sequence masked past the protein's catalytic
    domain end, so the catalytic-domain values describe the same residues. Every measure has a
    catalytic-domain ('cd') and a whole-structure ('all') value."""
    scales = hydrophobicity_scales()
    records = []
    for protein in proteins:
        lineage = _lineage(protein.organism)
        domain = ''
        if protein.cat_domain_length is not None:
            domain = ''.join(c for c in catalytic_sequence(protein) if c.isalpha())
        residues = getattr(protein, 'residue_profile', None) or {}
        surface_hydrophobicity = residues.get('surface_hydrophobicity') or {}
        uniprot_cbms = [c.strip() for c in (protein.associated_CBMs or '').split(',')
                        if c.strip() and c.strip().lower() not in ('nan', 'none', 'null')]
        site = second_sphere(protein)
        brace = brace_motifs(protein)
        records.append({
            'code': protein.protein_code,
            'family': protein.clade.name,
            'species': species_name(protein.organism.name),
            'kingdom': protein.organism.kingdom or '',
            'taxon': list(taxon_group(lineage, protein.organism.kingdom or '')),
            'lineage': lineage_text(lineage),
            'substrates': sorted(s.name for s in protein.substrate.all()),
            # those of them a paper reports (the rest follow from EC numbers or the name)
            'substrate_papers': sorted(substrates_from_papers(protein)),
            # [C1 | C4 | C1/C4 | Unknown, papers | ec | none]
            'regio': regioselectivity_of(protein),
            # CBM families the protein's UniProt entries annotate, and those inside its structure
            # uniprot_known: the protein has an active UniProt entry, so an empty uniprot list
            # means no CBM rather than nothing to look in (deleted entries, own models, CAZyme3D)
            'cbms': {'uniprot': uniprot_cbms, 'structure': structure_cbms(protein.selected_structure),
                     'uniprot_known': bool((getattr(protein, 'uniprot', None) or '').strip())},
            # second-sphere motif [short, with spacing] and the shielding residue
            'motif': [site['short'] or NOT_DETERMINED, site['motif']],
            'distal': site['distal'] or NOT_DETERMINED,
            # motifs after the N-terminal and the internal His: [motif, residues ('H-H-G-P')]
            'brace': {which: [m['motif'], m['residues']] for which, m in brace.items()},
            'structure': protein.selected_structure or '',
            'source': structure_source(protein.selected_structure),
            'measures': {
                'length': _scopes(protein.cat_domain_length, getattr(protein, 'sequence_length', None) or None),
                'mass': _scopes(getattr(protein, 'cat_domain_mass', None), getattr(protein, 'sequence_mass', None)),
                'charge': residues.get('charge') or _scopes(None, None),
                'surface_charge': residues.get('surface_charge') or _scopes(None, None),
                'surface_length': residues.get('surface_length') or _scopes(None, None),
                'surface_mass': residues.get('surface_mass') or _scopes(None, None),
                'hydrophobicity': {
                    name: _scopes(mean_hydropathy(domain, scale) if domain else None,
                                  mean_hydropathy(protein.sequence, scale))
                    for name, scale in scales.items()},
                'surface_hydrophobicity': {
                    name: surface_hydrophobicity.get(name) or _scopes(None, None) for name in scales},
                # the loop regions together and the substrate-binding residues (step 7c): both
                # lie inside the catalytic domain, so they carry one value, not a scope pair
                'sets': getattr(protein, 'residue_sets', None) or {},
            },
            'loops': loop_lengths(protein.loop_regions),
            # per loop: rmsf, charge [all, surface], hydrophobicity [[all, surface] per scale]
            'loop_values': getattr(protein, 'loop_profile', None) or {},
            # disulfides, b2, pI, glycosylation sequons, binding-surface make-up, hole hopping
            # (structural_features): scored by the selection overview's grouping analysis
            'features': getattr(protein, 'structural_features', None) or {},
        })
    families = sorted({r['family'] for r in records}, key=_family_order)
    return {
        'families': [{'name': f, 'color': CLADE_COLORS.get(f, '#569976')} for f in families],
        # stacked segments follow the tree's family order (see FAMILY_STACK_ORDER)
        'stack_order': ([f for f in FAMILY_STACK_ORDER if f in families]
                        + [f for f in families if f not in FAMILY_STACK_ORDER]),
        'taxon_groups': TAXON_GROUPS,
        'hydrophobicity_scales': [{'name': name, 'hydrophilic': name in HYDROPHILIC_SCALES}
                                  for name in scales],
        'surface_rsa': SURFACE_RSA,
        'not_determined': NOT_DETERMINED,
        'open_position': OPEN_POSITION,
        'other_motif': OTHER_MOTIF,
        'loop_names': (LOOP_ORDER + sorted({loop for r in records for loop in r['loops']}
                                           - set(LOOP_ORDER))),
        'proteins': records,
    }
