"""Carbohydrate-binding modules (CBMs): the CBM family of a domain annotation, the CBMs a
protein's UniProt entries annotate, and those inside the catalogue's structure.

InterPro and UniProt name CBMs in several ways, and some names mislead: the InterPro superfamily
"CBM2/CBM3, carbohydrate-binding domain superfamily" (IPR008965) holds every CBM2 as well, so a
name search called all CBM2 proteins "CBM2/CBM3"; UniProt's "Chitin-binding type-3" is CBM5/12
(IPR003610, "Chitin-binding domain type 3" in SMART), while its "Chitin-binding type-4" is the
AA10 catalytic domain itself (the former CBM33). The family therefore comes from the entries
below only; a superfamily names one when it holds a single family.
"""
import functools
import json
import re
from pathlib import Path

# InterPro entries of one CBM family
INTERPRO_FAMILIES = {
    'IPR000254': 'CBM1',      # Cellulose-binding domain
    'IPR035971': 'CBM1',      # Cellulose-binding domain superfamily
    'IPR001919': 'CBM2',      # Carbohydrate-binding type-2 domain
    'IPR012291': 'CBM2',      # CBM2, carbohydrate-binding domain superfamily
    'IPR001956': 'CBM3',      # Carbohydrate-binding module 3
    'IPR036966': 'CBM3',      # Carbohydrate-binding module 3 superfamily
    'IPR003610': 'CBM5/12',   # Carbohydrate binding module family 5/12
    'IPR036573': 'CBM5/12',   # Carbohydrate-binding module superfamily 5/12
    'IPR009031': 'CBM10',     # Carbohydrate binding module family 10
    'IPR036601': 'CBM10',     # CBM10 superfamily
    'IPR002044': 'CBM20',     # Carbohydrate binding module family 20
    'IPR034836': 'CBM20',     # Glucoamylase, CBM20 domain
}
# UniProt domain features ('uniprot' source in domain_annotations.json)
UNIPROT_FAMILIES = {
    'CBM1': 'CBM1', 'CBM2': 'CBM2', 'CBM3': 'CBM3', 'CBM10': 'CBM10', 'CBM20': 'CBM20',
    'Chitin-binding type-3': 'CBM5/12',
}

DOMAIN_RANGES_JSON = (Path(__file__).resolve().parents[2] / 'data' / 'django_outputs'
                      / 'domain_ranges_auth.json')


def cbm_family(domain):
    """The CBM family of a domain annotation ({source, accession, name}), or None."""
    if domain.get('source') == 'uniprot':
        return UNIPROT_FAMILIES.get((domain.get('name') or '').strip())
    return INTERPRO_FAMILIES.get(domain.get('accession'))


def family_order(family):
    """Sort key: CBM1, CBM2, …, CBM5/12, CBM10, CBM20."""
    m = re.match(r'CBM(\d+)', family or '')
    return (int(m.group(1)) if m else 999, family or '')


def annotated_cbms(entry):
    """The CBMs of a domain_annotations.json entry (UniProt numbering), in sequence order:
    [{family, name, accession, start, end}]; overlapping annotations of one family (InterPro
    domain, superfamily, UniProt feature) are one module, a repeated module is two."""
    found = []
    for domain in (entry or {}).get('domains', []):
        family = cbm_family(domain)
        if not family or domain.get('start') is None or domain.get('end') is None:
            continue
        found.append({'family': family, 'name': domain.get('name') or family,
                      'accession': domain.get('accession') or '', 'start': int(domain['start']),
                      'end': int(domain['end']), 'specific': domain.get('accession') in INTERPRO_FAMILIES
                      and 'superfamily' not in (domain.get('name') or '').lower()})
    modules = []
    for cbm in sorted(found, key=lambda c: (c['family'], c['start'], c['end'])):
        last = modules[-1] if modules else None
        if last and last['family'] == cbm['family'] and cbm['start'] <= last['end']:
            last['start'], last['end'] = min(last['start'], cbm['start']), max(last['end'], cbm['end'])
            if cbm['specific'] and not last['specific']:
                last.update(name=cbm['name'], accession=cbm['accession'], specific=True)
        else:
            modules.append(dict(cbm))
    return [{k: v for k, v in m.items() if k != 'specific'}
            for m in sorted(modules, key=lambda m: (m['start'], m['end']))]


def uniprot_cbms(accessions, annotations):
    """The CBM families annotated in any of the UniProt entries (`annotations`: the contents of
    domain_annotations.json), in family order."""
    families = {cbm['family'] for acc in accessions for cbm in annotated_cbms(annotations.get(acc))}
    return sorted(families, key=family_order)


@functools.lru_cache(maxsize=1)
def _domain_ranges():
    try:
        return json.loads(DOMAIN_RANGES_JSON.read_text())
    except (OSError, ValueError):
        return {}


def structure_cbms(struct_id, ranges=None):
    """The CBM families inside the catalogue's structure (map_domain_ranges: InterPro domains
    through SIFTS or the AlphaFold model's own entry, and UniProt-annotated CBMs aligned onto the
    structure), in family order. `ranges`: that command's output (default: its file)."""
    domains = (ranges if ranges is not None else _domain_ranges()).get((struct_id or '').strip(), [])
    families = {d.get('family') or INTERPRO_FAMILIES.get(d.get('accession')) or 'CBM'
                for d in domains if d.get('type') == 'CBM'}
    return sorted(families, key=family_order)
