"""Protein mass from sequence (kDa), for the catalogue's size metadata and filters."""
from Bio.SeqUtils import molecular_weight

_STANDARD = set('ACDEFGHIKLMNPQRSTVWY')


def polypeptide_kda(sequence):
    """Average molecular mass of the polypeptide in kDa (Biopython average masses, one water for
    the termini), counting only the 20 standard residues; gaps, whitespace and unknown letters
    are ignored. None when no residue is left. Glycans, metal ions and other modifications are
    not included."""
    residues = ''.join(c for c in (sequence or '').upper() if c in _STANDARD)
    if not residues:
        return None
    return molecular_weight(residues, seq_type='protein') / 1000


_WATER = 18.01528   # Da, the one water molecular_weight adds for the chain termini


def residues_kda(residues):
    """Summed average mass in kDa of a set of residues that need not be one chain (the surface
    residues, the loop regions, the substrate-binding residues): the residue masses without the
    terminal water polypeptide_kda counts. None when no standard residue is given."""
    residues = ''.join(c for c in (residues or '').upper() if c in _STANDARD)
    if not residues:
        return None
    return (molecular_weight(residues, seq_type='protein') - _WATER) / 1000
