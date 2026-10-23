"""The second-sphere hydrogen-bonding motif of a protein and its shielding residue, the residue
over the copper's solvent side (Protein.second_sphere_motif, pipeline step 8c-8d), as the
browse filters, the overview chart and the compare page's Sequence plot show them."""
import json
from collections import Counter

# Labels of the proteins without a motif (no histidine brace) and of an open axial position
NOT_DETERMINED = 'Not determined'
OPEN_POSITION = 'None (open)'

# What each motif residue is, for hover texts; the aromatic has its own category (axial position)
ROLE_LABELS = {
    'first': 'holds the Gln/Glu',
    'gln_glu': 'Gln/Glu near the copper',
    'middle': 'two residues before the axial aromatic',
    'aromatic': 'axial aromatic',
}
FOUND_BY_LABELS = {
    'stacked on His-2': 'stacked on His-2',
    'sequence motif': 'published sequence motif',
    'hydrogen bond': 'hydrogen bond',
}


def _residue_label(name):
    """'ALA' -> 'Ala'."""
    return (name or '').capitalize()


def second_sphere(protein):
    """{'motif': 'H-X8-Q-X-Y', 'short': 'H-Q-Y', 'residues': [{'role', 'name', 'residue', ...}],
    'distal': 'Ala' | OPEN_POSITION | '', 'distal_residue': {...} or None} of a Protein (or of
    the field's JSON text); empty values when the protein has none."""
    text = protein if isinstance(protein, str) or protein is None else \
        getattr(protein, 'second_sphere_motif', None)
    try:
        data = json.loads(text) if text else {}
    except (TypeError, ValueError):
        data = {}
    distal = data.get('distal_axial') if isinstance(data, dict) else None
    if not isinstance(distal, dict):
        distal_label, distal_residue = '', None
    elif distal.get('name'):
        distal_label, distal_residue = _residue_label(distal['name']), distal
    else:
        distal_label, distal_residue = OPEN_POSITION, None
    return {
        'motif': data.get('motif') or '',
        'short': data.get('short') or '',
        'residues': [r for r in data.get('residues') or [] if isinstance(r, dict)],
        'distal': distal_label,
        'distal_residue': distal_residue,
    }


def motif_residue_text(motif, residue):
    """Hover text of a motif residue: 'Second-sphere motif H-X8-Q-X-Y: holds the Gln/Glu
    (stacked on His-2)'."""
    role = ROLE_LABELS.get(residue.get('role'), residue.get('role') or '')
    found_by = FOUND_BY_LABELS.get(residue.get('found_by'))
    if residue.get('cu_distance') is not None:
        found_by = f"{residue['cu_distance']:.1f} Å from the copper"
    return f"Second-sphere motif {motif}: {role}" + (f' ({found_by})' if found_by else '')


def filter_options(values):
    """[(value, count)] of the non-empty values, most common first, NOT_DETERMINED and
    OPEN_POSITION last."""
    counts = Counter(v for v in values if v)
    last = {NOT_DETERMINED: 2, OPEN_POSITION: 1}
    return sorted(counts.items(), key=lambda item: (last.get(item[0], 0), -item[1], item[0]))
