"""Produce per-structure domain ranges (catalytic domain + CBMs) in the structure's own
auth residue numbering, for the domain bar above the sequence alignment.

Numbering frames
----------------
Structures keep their deposited auth_seq_id; we never renumber. InterPro/UniProt
annotations (`domain_annotations.json`) are in UniProt *precursor* numbering, which only
coincides with auth numbering for AlphaFold models and for crystal structures the
depositor happened to number that way (offsets of -1..-72 observed otherwise).

Rather than doing offset arithmetic, we take domain ranges from the source that is already
in the right frame:

  * crystal structures -> PDBe SIFTS `/mappings/{pdb}` InterPro section, which reports
    boundaries directly in author residue numbers, per chain.
  * AlphaFold models    -> `domain_annotations.json` (AF models are numbered by UniProt
    position, so auth == UniProt).
  * neither available   -> fall back to the per-clade MSA-consensus catalytic domain end
    already computed by the pipeline in `catalytic_domain_boundaries.json`.

A CBM that one of the protein's UniProt entries annotates but the source above does not give
(the model's own entry is deleted or unannotated while a strain variant's or an older entry's
is, or only UniProt names the module, as in GbpA) is placed by aligning that entry's sequence
with the structure's; it counts when at least half of its residues align. Each CBM carries its
family (`catalogue.cbms`), which the overview counts.

The catalytic domain always *starts* at the brace histidine (the first resolved residue) —
Pfam/InterPro boundaries sometimes clip a few N-terminal residues, which would cut off
His1. This is why the pipeline only ever stored `end_residue`.

Input: pipeline files `data/02_data/domain_annotations.json` and
`data/06_catalytic_domain/catalytic_domain_boundaries.json`.
Output: `data/django_outputs/domain_ranges_auth.json` (SIFTS responses and UniProt
sequences cached in `sifts_cache.json` and `uniprot_sequences.json` next to it)
    {struct_id: [{"type": "catalytic"|"CBM", "name": str, "accession": str,
                  "start_residue": int, "end_residue": int, "source": str,
                  "family": "CBM1"|... (CBMs only)}, ...]}
"""

import json
import time
import urllib.error
import urllib.request
from pathlib import Path

from Bio import Align
from Bio.Align import substitution_matrices
from django.core.management.base import BaseCommand

from catalogue.cbms import INTERPRO_FAMILIES, annotated_cbms
from catalogue.makePlots import apply_catalytic_boundary, parse_pdb_residue_numbers
from catalogue.models import Protein

SIFTS_URL = 'https://www.ebi.ac.uk/pdbe/api/mappings/{pdb}'
UNIPROT_FASTA_URL = 'https://rest.uniprot.org/uniprotkb/{accession}.fasta'
# an entry's CBMs are placed only on a structure of the same protein
MIN_ALIGNED_IDENTITY = 0.9

# InterPro accessions whose span coincides with the annotated catalytic domain.
# Derived empirically from domain_annotations.json (see the analysis in git history):
# the AA10 catalytic module is called "Cellulose/chitin-binding protein, N-terminal"
# (IPR004302) — a binding-sounding name that a naive regex would misfile as a CBM.
# Ordered by preference: family-specific domain first, generic fold last.
CATALYTIC_ACCESSIONS = [
    'IPR005103',   # Auxiliary Activity family 9, catalytic domain
    'IPR004302',   # Cellulose/chitin-binding protein, N-terminal  (AA10 catalytic module)
    'IPR054497',   # Lytic polysaccharide monooxygenase AA14
    'IPR014756',   # Immunoglobulin E-set (fold spanning the catalytic domain)
]

# Genuine carbohydrate-binding modules (never overlap the catalytic domain).
CBM_ACCESSIONS = {
    'IPR000254': 'Cellulose-binding domain (CBM1)',
    'IPR035971': 'Cellulose-binding domain superfamily',
    'IPR001919': 'Carbohydrate-binding type-2 domain (CBM2)',
    'IPR012291': 'CBM2, carbohydrate-binding domain superfamily',
    'IPR008965': 'CBM2/CBM3, carbohydrate-binding domain superfamily',
    'IPR001956': 'Carbohydrate-binding module 3',
    'IPR036966': 'Carbohydrate-binding module 3 superfamily',
    'IPR003610': 'Carbohydrate binding module family 5/12',
    'IPR036573': 'Carbohydrate-binding module superfamily 5/12',
    'IPR009031': 'Carbohydrate binding module family 10',
    'IPR002883': 'CBM10/dockerin domain',
    'IPR036601': 'CBM10 superfamily',
    'IPR002044': 'Carbohydrate binding module family 20',
    'IPR034836': 'Glucoamylase, CBM20 domain',
}


def _merge(items):
    """Merge overlapping/adjacent CBM intervals, keeping the most specific name."""
    if not items:
        return []
    items = sorted(items, key=lambda d: (d['start_residue'], d['end_residue']))
    out = [dict(items[0])]
    for cur in items[1:]:
        last = out[-1]
        if cur['start_residue'] <= last['end_residue'] + 1:
            last['start_residue'] = min(last['start_residue'], cur['start_residue'])
            last['end_residue'] = max(last['end_residue'], cur['end_residue'])
            if _specificity(cur['name']) > _specificity(last['name']):
                last['name'], last['accession'] = cur['name'], cur['accession']
            last['family'] = last.get('family') or cur.get('family')
        else:
            out.append(dict(cur))
    return out


def _specificity(name):
    n = (name or '').lower()
    return 0 if ('superfamily' in n or '-like fold' in n) else 1


_ALIGNER = Align.PairwiseAligner(mode='local', substitution_matrix=substitution_matrices.load('BLOSUM62'),
                                 open_gap_score=-10, extend_gap_score=-0.5)


def _position_map(reference, sequence):
    """{reference index: sequence index} of the aligned residues (local alignment), and the
    identity of those pairs."""
    alignment = _ALIGNER.align(reference, sequence)[0]
    mapping, same = {}, 0
    for (r0, r1), (s0, s1) in zip(*alignment.aligned):
        for i, j in zip(range(r0, r1), range(s0, s1)):
            mapping[i] = j
            same += reference[i] == sequence[j]
    return mapping, (same / len(mapping) if mapping else 0.0)


class Command(BaseCommand):
    help = "Map catalytic domain + CBM ranges onto each structure's auth residue numbering."

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true')

    def handle(self, *args, **opts):
        data_dir = Path(__file__).resolve().parents[4] / 'data'
        annotations = json.loads((data_dir / '02_data' / 'domain_annotations.json').read_text())
        cat_bounds = json.loads((data_dir / '06_catalytic_domain' / 'catalytic_domain_boundaries.json').read_text())
        out_dir = data_dir / 'django_outputs'
        out_dir.mkdir(exist_ok=True)
        cache_path = out_dir / 'sifts_cache.json'
        sequences_path = out_dir / 'uniprot_sequences.json'
        out_path = out_dir / 'domain_ranges_auth.json'
        cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
        sequences = json.loads(sequences_path.read_text()) if sequences_path.exists() else {}

        result = {}
        stats = {'sifts': 0, 'af': 0, 'fallback': 0, 'skipped': 0, 'aligned': 0}

        for prot in Protein.objects.exclude(pdb_structure='').exclude(pdb_structure__isnull=True):
            struct_id = (prot.selected_structure or '').strip()
            res_nums = parse_pdb_residue_numbers(prot.pdb_structure)
            if not struct_id or not res_nums:
                stats['skipped'] += 1
                continue

            first_res, last_res = res_nums[0], res_nums[-1]
            chain = next((l[21:22] for l in prot.pdb_structure.split('\n')
                          if l.startswith('ATOM')), 'A')

            if struct_id.upper().startswith('AF-'):
                domains = self._from_annotations(annotations.get(struct_id[3:]), first_res, last_res)
                if domains:
                    stats['af'] += 1
            elif struct_id.upper().startswith(('C3D-', 'AF3-')):
                # CAZyme3D model or a model made for the catalogue (step 1e): no UniProt
                # annotation and no SIFTS entry, so the family consensus below decides
                domains = []
            else:
                domains = self._from_sifts(struct_id, chain, first_res, last_res, cache)
                if domains:
                    stats['sifts'] += 1

            # A structure can pick up CBMs but no catalytic domain (missing annotation).
            # Fill the catalytic block from the per-clade MSA consensus in that case too,
            # not only when nothing at all was found.
            if not any(d['type'] == 'catalytic' for d in domains):
                end = cat_bounds.get(struct_id, {}).get('end_residue')
                if end is not None:
                    domains.insert(0, {
                        'type': 'catalytic', 'name': 'Catalytic domain', 'accession': None,
                        'start_residue': first_res, 'end_residue': int(end),
                        'source': 'consensus',
                    })
                    stats['fallback'] += 1

            # CBMs of the protein's UniProt entries that the sources above did not give
            known = [d for d in domains if d['type'] == 'CBM']
            for cbm in self._from_uniprot_entries(prot, res_nums, annotations, sequences):
                if not any(cbm['start_residue'] <= k['end_residue'] and k['start_residue'] <= cbm['end_residue']
                           for k in known):
                    known.append(cbm)
                    stats['aligned'] += 1
                    self.stdout.write(f"  {struct_id}: {cbm['family']} {cbm['start_residue']}-"
                                      f"{cbm['end_residue']} from {cbm['source']}")
            domains = ([d for d in domains if d['type'] != 'CBM']
                       + sorted(known, key=lambda d: d['start_residue']))

            if not domains:
                stats['skipped'] += 1
                continue

            # The catalytic domain ends where catalytic_domain_boundaries.json says (the assignment
            # the app uses); annotated CBMs are clipped to start after it.
            result[struct_id] = apply_catalytic_boundary(domains, cat_bounds.get(struct_id, {}).get('end_residue'))

        n_cbm = sum(1 for v in result.values() if any(d['type'] == 'CBM' for d in v))
        self.stdout.write(
            f"SIFTS (crystal): {stats['sifts']}  |  AF: {stats['af']}  |  "
            f"consensus fallback: {stats['fallback']}  |  skipped: {stats['skipped']}  |  "
            f"CBMs placed by alignment: {stats['aligned']}"
        )
        self.stdout.write(f"total structures: {len(result)}  |  with >=1 CBM: {n_cbm}")

        if opts['dry_run']:
            self.stdout.write(self.style.WARNING('Dry run — nothing written.'))
            for sid in list(result)[:6]:
                self.stdout.write(f'  {sid}:')
                for d in result[sid]:
                    self.stdout.write('     %-9s %4d-%-4d  %s (%s)' % (
                        d['type'], d['start_residue'], d['end_residue'], d['name'], d['source']))
            return

        out_path.write_text(json.dumps(result, indent=2))
        cache_path.write_text(json.dumps(cache))
        sequences_path.write_text(json.dumps(sequences, indent=0, sort_keys=True))
        self.stdout.write(self.style.SUCCESS(f'Wrote {out_path} ({len(result)} structures)'))

    # ── crystal structures: SIFTS InterPro, already in auth numbering ────────
    def _from_sifts(self, pdb_id, chain, first_res, last_res, cache):
        key = pdb_id.lower()
        if key not in cache:
            try:
                with urllib.request.urlopen(SIFTS_URL.format(pdb=key), timeout=30) as r:
                    cache[key] = json.loads(r.read().decode())
                time.sleep(0.1)
            except (urllib.error.URLError, ValueError, TimeoutError) as exc:
                self.stderr.write(f'  SIFTS failed for {pdb_id}: {exc}')
                cache[key] = None

        data = (cache.get(key) or {}).get(key, {})
        interpro = data.get('InterPro') or {}
        if not interpro:
            return []

        def spans(acc):
            info = interpro.get(acc)
            if not info:
                return []
            out = []
            for m in info.get('mappings', []):
                if m.get('chain_id') != chain:
                    continue
                s = (m.get('start') or {}).get('author_residue_number')
                e = (m.get('end') or {}).get('author_residue_number')
                if s is not None and e is not None:
                    out.append((s, e, info.get('name') or acc))
            return out

        domains = []
        for acc in CATALYTIC_ACCESSIONS:
            sp = spans(acc)
            if sp:
                _, end, name = sp[0]
                domains.append({
                    'type': 'catalytic', 'name': name, 'accession': acc,
                    # always start at the brace His — InterPro can clip N-terminal residues
                    'start_residue': first_res,
                    'end_residue': min(end, last_res),
                    'source': 'sifts',
                })
                break

        cbms = []
        for acc, label in CBM_ACCESSIONS.items():
            for s, e, name in spans(acc):
                s, e = max(s, first_res), min(e, last_res)
                if s <= e:
                    cbms.append({'type': 'CBM', 'name': name or label, 'accession': acc,
                                 'start_residue': s, 'end_residue': e, 'source': 'sifts',
                                 'family': INTERPRO_FAMILIES.get(acc)})
        domains.extend(_merge(cbms))
        return domains

    # ── AlphaFold models: auth == UniProt, so annotations apply directly ─────
    def _from_annotations(self, entry, first_res, last_res):
        if not entry:
            return []
        domains = []

        cat = entry.get('catalytic_domain')
        if cat and cat.get('end') is not None:
            end = min(int(cat['end']), last_res)
            if end >= first_res:
                domains.append({
                    'type': 'catalytic', 'name': cat.get('name') or 'Catalytic domain',
                    'accession': cat.get('pfam_id'),
                    'start_residue': first_res, 'end_residue': end, 'source': 'annotation',
                })

        cbms = []
        for d in entry.get('domains', []):
            acc = d.get('accession')
            if d.get('source') != 'interpro' or acc not in CBM_ACCESSIONS:
                continue
            if d.get('start') is None or d.get('end') is None:
                continue
            s, e = max(int(d['start']), first_res), min(int(d['end']), last_res)
            if s <= e:
                cbms.append({'type': 'CBM', 'name': d.get('name') or CBM_ACCESSIONS[acc],
                             'accession': acc, 'start_residue': s, 'end_residue': e,
                             'source': 'annotation', 'family': INTERPRO_FAMILIES.get(acc)})
        domains.extend(_merge(cbms))
        return domains

    # ── CBMs of any UniProt entry of the protein, placed by sequence alignment ─
    def _from_uniprot_entries(self, prot, res_nums, annotations, sequences):
        """The CBMs annotated in the protein's UniProt entries (catalogue.cbms.annotated_cbms)
        that lie in the structure: each entry's sequence is aligned with the structure's, and a
        module counts when at least half of its residues align. Auth numbering."""
        sequence = prot.sequence or ''
        if len(sequence) != len(res_nums):
            return []
        try:
            accessions = list(json.loads(prot.uniprot_all or '{}'))
        except ValueError:
            accessions = []
        found = []
        for accession in accessions:
            cbms = annotated_cbms(annotations.get(accession))
            reference = self._uniprot_sequence(accession, sequences) if cbms else None
            if not reference:
                continue
            mapping, identity = _position_map(reference, sequence)
            if identity < MIN_ALIGNED_IDENTITY:
                self.stderr.write(f'  {prot.selected_structure}: {accession} aligns with only '
                                  f'{identity:.0%} identity, its CBMs are not placed')
                continue
            for cbm in cbms:
                placed = [mapping[i] for i in range(cbm['start'] - 1, cbm['end']) if i in mapping]
                if 2 * len(placed) < cbm['end'] - cbm['start'] + 1:
                    continue
                found.append({'type': 'CBM', 'name': cbm['name'], 'accession': cbm['accession'] or None,
                              'start_residue': res_nums[min(placed)], 'end_residue': res_nums[max(placed)],
                              'source': f'UniProt {accession}, aligned', 'family': cbm['family']})
        # entries of one protein annotate the same modules
        return _merge(found)

    def _uniprot_sequence(self, accession, cache):
        """The entry's sequence (cached; a failed lookup is retried on the next run)."""
        if accession not in cache:
            try:
                with urllib.request.urlopen(UNIPROT_FASTA_URL.format(accession=accession), timeout=30) as r:
                    text = r.read().decode()
                time.sleep(0.1)
            except (urllib.error.URLError, ValueError, TimeoutError) as exc:
                self.stderr.write(f'  UniProt sequence failed for {accession}: {exc}')
                return None
            cache[accession] = ''.join(line.strip() for line in text.splitlines() if not line.startswith('>'))
        return cache[accession]
