"""The sequence motifs after the two histidines of the brace (Protein.brace_motifs, pipeline step
8e) -- after the N-terminal His (N-His: H-G-X4-P-X2-R, H-G-X4-P, H-G) and after the internal His
(Internal-His: H-X-G-P, H-X2-G-P, H-X-G) -- as the browse filters, the overview chart and the
compare page's Sequence plot show them."""
import json
from collections import Counter

from .second_sphere import NOT_DETERMINED

# The two histidines. Not His-1/His-2: the N-terminal His is residue 1 only in the mature
# numbering (most structures keep their signal peptide's or their construct's numbers).
HIS_LABELS = {'n_his': 'N-His', 'internal_his': 'Internal-His'}
# a His with none of the motifs after it
OTHER = 'Other'
ROLE_LABELS = {'gly': 'Gly', 'pro': 'Pro', 'arg': 'Arg'}


def _residues_label(motif, sequence):
    """The His and the next three residues, 'H-G-Y-V' (the next four after an H-X2-G-P)."""
    return '-'.join(sequence[:5 if motif == 'H-X2-G-P' else 4])


def brace_motifs(protein):
    """{'n_his': {...}, 'internal_his': {...}} of a Protein (or of the field's JSON text), each
    {'motif': 'H-X-G-P' | OTHER | NOT_DETERMINED, 'residues': 'H-H-G-P' or '', 'sequence':
    'HHGPV' or '', 'motif_residues': [{'role', 'name', 'residue', 'offset'}]}."""
    text = protein if isinstance(protein, str) or protein is None else \
        getattr(protein, 'brace_motifs', None)
    try:
        data = json.loads(text) if text else {}
    except (TypeError, ValueError):
        data = {}
    result = {}
    for which in HIS_LABELS:
        entry = data.get(which) if isinstance(data, dict) else None
        if not isinstance(entry, dict) or not entry.get('sequence'):
            result[which] = {'motif': NOT_DETERMINED, 'residues': '', 'sequence': '',
                             'motif_residues': []}
            continue
        motif = entry.get('motif') or OTHER
        result[which] = {
            'motif': motif,
            'residues': _residues_label(motif, entry['sequence']),
            'sequence': entry['sequence'],
            'motif_residues': [r for r in entry.get('residues') or [] if isinstance(r, dict)],
        }
    return result


def motif_residue_text(which, motif, residue):
    """Hover text of a motif residue: 'Internal-His motif H-X-G-P: Gly (Internal-His +2)'."""
    role = ROLE_LABELS.get(residue.get('role'), (residue.get('name') or '').capitalize())
    position = f" ({HIS_LABELS[which]} +{residue['offset']})" if residue.get('offset') else ''
    return f'{HIS_LABELS[which]} motif {motif}: {role}{position}'


def filter_options(values):
    """[(value, count)] of the non-empty values, most common first, OTHER and NOT_DETERMINED
    last."""
    counts = Counter(v for v in values if v)
    last = {NOT_DETERMINED: 2, OTHER: 1}
    return sorted(counts.items(), key=lambda item: (last.get(item[0], 0), -item[1], item[0]))
