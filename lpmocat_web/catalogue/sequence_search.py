"""
Sequence search for the browse page's search bar: find catalogue LPMOs similar to a pasted
protein sequence.

The mode is picked by query length:
  * alignment (>= MIN_ALIGN_LENGTH residues): Smith-Waterman local alignment of the query
    against every catalogue sequence (Biopython PairwiseAligner, BLOSUM62, gap open 11 /
    extend 1 as in BLASTP). Hits are ranked by score. Significance is an approximate
    BLAST-style E-value from Karlin-Altschul statistics (lambda/K for BLOSUM62 11/1) with the
    whole catalogue as the database -- no edge-effect or composition correction.
  * motif (shorter queries): exact match anywhere in the sequence; X = any residue,
    B = D/N, Z = E/Q, J = I/L.

The catalogue is small (~125 sequences, ~33k residues), so a full scan takes ~0.1 s and no
index or external BLAST install is needed.
"""
import math
import re

from Bio.Align import PairwiseAligner, substitution_matrices

MIN_ALIGN_LENGTH = 20
MIN_MOTIF_LENGTH = 3
MAX_QUERY_LENGTH = 3000
# Hits at or below this E-value count as homologues (list/dendrogram highlight, "select all").
SIGNIFICANT_EVALUE = 1e-3
# When nothing is significant, list this many closest matches (flagged significant=False).
N_WEAK_SHOWN = 3

# Karlin-Altschul parameters for BLOSUM62, gap open 11 / extend 1 (NCBI BLASTP defaults).
_LAMBDA = 0.267
_K = 0.041

_BLOSUM62 = substitution_matrices.load('BLOSUM62')
_VALID = set('ACDEFGHIKLMNPQRSTVWYBZXUOJ')
# Selenocysteine / pyrrolysine / I-or-L are not in BLOSUM62 -> score them as unknown.
_TO_BLOSUM = str.maketrans('UOJ', 'XXX')
_MOTIF_CLASSES = {'X': '.', 'B': '[DN]', 'Z': '[EQ]', 'J': '[IL]'}


class QueryError(ValueError):
    """The pasted query is not a usable protein sequence."""


def parse_query(text):
    """Pasted text (plain or FASTA) -> (header, sequence).

    Only the first FASTA record is used. Whitespace, digits (GenBank-style numbering), gap
    characters and a trailing stop '*' are dropped. Raises QueryError otherwise."""
    lines = text.strip().splitlines()
    header = ''
    if lines and lines[0].startswith('>'):
        header = lines[0][1:].strip()
        body = []
        for line in lines[1:]:
            if line.startswith('>'):
                break
            body.append(line)
        lines = body
    seq = re.sub(r'[\s\d.\-]', '', ''.join(lines)).upper().rstrip('*')
    if not seq:
        raise QueryError('The query contains no residues.')
    bad = sorted(set(seq) - _VALID)
    if bad:
        raise QueryError('Not a protein sequence - unexpected character(s): ' + ' '.join(bad))
    if len(seq) > MAX_QUERY_LENGTH:
        raise QueryError(f'The query is too long ({len(seq)} residues; max {MAX_QUERY_LENGTH}).')
    if len(seq) >= MIN_ALIGN_LENGTH and sum(seq.count(c) for c in 'ACGTUN') >= 0.95 * len(seq):
        raise QueryError('This looks like a nucleotide sequence - paste the protein '
                         '(amino-acid) sequence instead.')
    return header, seq


def evalue(score, query_len, db_len):
    return _K * query_len * db_len * math.exp(-_LAMBDA * score)


def bit_score(score):
    return (_LAMBDA * score - math.log(_K)) / math.log(2)


def _aligner():
    # Biopython scores a length-k gap as open + (k-1)*extend, so open=-12 == BLAST's 11/1.
    return PairwiseAligner(mode='local', substitution_matrix=_BLOSUM62,
                           open_gap_score=-12, extend_gap_score=-1)


def align_search(seq, targets):
    """Rank targets [(code, sequence), ...] by local-alignment score against seq.

    Returns every significant hit (E <= SIGNIFICANT_EVALUE), best first; if there is none,
    the N_WEAK_SHOWN closest matches instead, flagged significant=False."""
    query = seq.translate(_TO_BLOSUM)
    targets = [(code, s.translate(_TO_BLOSUM)) for code, s in targets if s]
    db_len = sum(len(s) for _, s in targets)
    aligner = _aligner()
    ranked = sorted(((aligner.score(query, s), code, s) for code, s in targets),
                    key=lambda r: (-r[0], r[1]))
    significant = [r for r in ranked if evalue(r[0], len(query), db_len) <= SIGNIFICANT_EVALUE]
    chosen = significant or [r for r in ranked[:N_WEAK_SHOWN] if r[0] > 0]

    hits = []
    for score, code, target in chosen:
        aln = aligner.align(query, target)[0]
        counts = aln.counts()
        q_blocks, t_blocks = aln.aligned
        q_start, q_end = int(q_blocks[0][0]), int(q_blocks[-1][1])
        t_start, t_end = int(t_blocks[0][0]), int(t_blocks[-1][1])
        length = counts.identities + counts.mismatches + counts.gaps
        e = evalue(score, len(query), db_len)
        hits.append({
            'code': code,
            'score': score,
            'bits': round(bit_score(score), 1),
            'evalue': e,
            'significant': e <= SIGNIFICANT_EVALUE,
            'identity': round(100 * counts.identities / length, 1),
            'aligned_length': length,
            'query_start': q_start + 1, 'query_end': q_end,
            'target_start': t_start + 1, 'target_end': t_end,
            'query_coverage': round(100 * (q_end - q_start) / len(query), 1),
            'target_coverage': round(100 * (t_end - t_start) / len(target), 1),
        })
    return hits


def motif_search(seq, targets):
    """Targets [(code, sequence), ...] containing the motif, with 1-based match starts
    (overlapping matches included), ordered by protein code."""
    pattern = re.compile('(?=' + ''.join(_MOTIF_CLASSES.get(c, c) for c in seq) + ')')
    hits = []
    for code, s in targets:
        positions = [m.start() + 1 for m in pattern.finditer(s or '')]
        if positions:
            hits.append({'code': code, 'positions': positions, 'significant': True})
    return sorted(hits, key=lambda h: h['code'].lower())


def search(seq, targets):
    """Motif match for short queries, local-alignment ranking otherwise."""
    if len(seq) < MIN_ALIGN_LENGTH:
        if len(seq) < MIN_MOTIF_LENGTH:
            raise QueryError(f'Type at least {MIN_MOTIF_LENGTH} residues for a motif search.')
        return {'mode': 'motif', 'hits': motif_search(seq, targets)}
    return {'mode': 'alignment', 'hits': align_search(seq, targets)}
