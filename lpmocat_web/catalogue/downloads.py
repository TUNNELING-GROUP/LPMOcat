"""What the site offers for download: the zip of one comparison (compare page) and the files of
the whole catalogue (the Downloads page), both made of the same parts:

  proteins.csv            one row per protein: the compare page's table as plain values
  substrate_evidence.csv  every activity claim with its source and papers (step 2j)
  residues_*.csv          one row per residue: the plots' data along the alignment, merged
  alignment_*.fasta       the structural alignment
  structures/             the PDB files, full length and catalytic domain

The catalogue files are built by `manage.py build_downloads` after a data reload into
data/django_outputs/downloads/ (so they travel in the data bundle), together with release.json,
the data release the footer of every page shows."""
import csv
import datetime
import html
import io
import json
import re
import zipfile
from pathlib import Path

DOWNLOADS_DIR = Path(__file__).resolve().parents[2] / 'data' / 'django_outputs' / 'downloads'
RELEASE_FILE = DOWNLOADS_DIR / 'release.json'

VIEWS = (('catalytic', 'catalytic_domain'), ('full', 'full_length'))

LICENCE_TEXT = ('Licence: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/). The data are '
                'built from public resources (RCSB PDB, AlphaFold DB, UniProt, InterPro, BRENDA, '
                'CAZy, CAZyme3D, NCBI Taxonomy, SIFTS, Europe PMC), which keep their own terms; '
                'please cite them as well. Structure models made with AlphaFold Server (structure '
                'ids starting with AF3-) are for non-commercial use only (AlphaFold Server Output '
                'Terms of Use).')


# --- the data release --------------------------------------------------------------------------

_release_cache = {'mtime': None, 'data': None}


def release_info():
    """release.json of the catalogue downloads ({} before build_downloads has run), re-read when
    the file changes."""
    try:
        mtime = RELEASE_FILE.stat().st_mtime
    except OSError:
        return {}
    if _release_cache['mtime'] != mtime:
        try:
            _release_cache['data'] = json.loads(RELEASE_FILE.read_text())
        except (OSError, ValueError):
            _release_cache['data'] = {}
        _release_cache['mtime'] = mtime
    return _release_cache['data'] or {}


def release_label():
    release = release_info().get('release')
    return f'Data release {release}' if release else ''


# --- tables ------------------------------------------------------------------------------------

PROTEIN_HEADER = [
    'protein', 'other_names', 'family', 'organism', 'taxonomy', 'structure', 'structure_source',
    'uniprot', 'uniparc', 'ec_numbers', 'substrates', 'regioselectivity',
    'active_site_metal', 'n_his', 'internal_his', 'n_his_motif', 'internal_his_motif',
    'second_sphere_motif', 'shielding_residue', 'binding_aromatics', 'cbms', 'sequence_length',
    'catalytic_domain_length', 'catalytic_domain_end', 'mass_kda', 'catalytic_domain_mass_kda',
    'disulfide_bonds', 'references',
]


def _references_text(references_html):
    """The compare table's References cell as plain text: 'label <link>' per paper, '; '-joined."""
    out = []
    for href, label in re.findall(r'<a href="([^"]+)"[^>]*>(.*?)</a>', references_html or '', flags=re.S):
        item = f'{html.unescape(re.sub("<[^>]+>", "", label)).strip()} <{html.unescape(href)}>'
        if item not in out:
            out.append(item)
    return '; '.join(out)


def protein_rows(proteins, group_of=None):
    """[PROTEIN_HEADER] + one row per protein, in the order given; with group_of ({code: group
    name}, a comparison with groups) a first column 'group'."""
    from . import views
    from .activities import activity_data, activity_references
    from .overview import (_lineage, lineage_text, protein_binding_aromatics, regioselectivity_label,
                           structure_source)

    rows = [(['group'] if group_of else []) + PROTEIN_HEADER]
    for p in proteins:
        names = [n.strip() for n in (p.other_names or '').split(';')
                 if n.strip() and n.strip() != p.protein_code]
        site = views._site_cells(p)
        brace = re.findall(r'\d+', p.histidine_brace or '')
        disulfides = []
        try:
            disulfides = [f'Cys{a}-Cys{b}' for a, b, *_ in json.loads(p.disulfides or '[]')]
        except (ValueError, TypeError):
            pass
        masses = views._protein_masses(p)
        record = activity_data(p)
        references_html, _ = views._format_references(p, activity_references(record))
        lineage = _lineage(p.organism) if p.organism else None
        rows.append(([group_of.get(p.protein_code, '')] if group_of else []) + [
            p.protein_code,
            '; '.join(names),
            p.clade.name if p.clade else '',
            p.organism.name if p.organism else '',
            (lineage_text(lineage) or p.organism.kingdom or '') if p.organism else '',
            p.selected_structure or '',
            structure_source(p.selected_structure),
            p.uniprot or '',
            p.uniparc or '',
            p.ec_numbers or '',
            '; '.join(s.name for s in p.substrate.all()),
            regioselectivity_label(r.name for r in p.regioselectivity.all()),
            p.metal or '',
            brace[0] if len(brace) > 0 else '',
            brace[1] if len(brace) > 1 else '',
            site.get('n_his', {}).get('text', ''),
            site.get('internal_his', {}).get('text', ''),
            site['second_sphere']['text'] or '',
            site['shielding']['text'] or '',
            '; '.join(f"{a['label']} ({a['region']})" for a in (protein_binding_aromatics(p) or [])),
            p.associated_CBMs or '',
            p.sequence_length if p.sequence_length is not None else '',
            views._catalytic_domain_length(p) or '',
            views._domain_boundaries.get(p.selected_structure or '', {}).get('end_residue') or '',
            views._kda_text(masses[1]),
            views._kda_text(masses[0]),
            '; '.join(disulfides),
            _references_text(references_html),
        ])
    return rows


def substrate_rows(proteins):
    from .activities import EXPORT_HEADER, export_rows
    return [EXPORT_HEADER] + export_rows(proteins)


def residue_csv(export_view, scale_name):
    """One per-residue CSV from the plots' export of one view (makePlots.makeExportData()[view]):
    the secondary structure, core strand, loop, binding-surface and disulfide columns, then the
    relative SASA, charge, hydrophobicity and every flexibility method, merged by protein and
    sequence position."""
    from .makePlots import SS_LOOPS_PLOT_TYPE

    parts = [(SS_LOOPS_PLOT_TYPE, None), ('Relative Surface Accessibility', None),
             ('Charge', None), ('Hydrophobicity', f'hydrophobicity_{_slug(scale_name)}'),
             ('Flexibility', 'flexibility')]
    merged, order, columns = {}, [], ['protein', 'residue_number', 'position', 'residue']
    for plot_type, rename in parts:
        item = export_view.get(plot_type)
        if not item:
            continue
        reader = csv.reader(io.StringIO(item['content']))
        header = next(reader)
        extra = header[4:]
        if rename == 'flexibility':
            extra_names = [f'flexibility_{_slug(m)}' for m in extra]
        elif rename:
            extra_names = [rename]
        else:
            extra_names = extra
        columns += extra_names
        for row in reader:
            key = (row[0], row[2])
            if key not in merged:
                merged[key] = dict(zip(columns[:4], row[:4]))
                order.append(key)
            merged[key].update(zip(extra_names, row[4:]))
    out = io.StringIO()
    writer = csv.writer(out, lineterminator='\n')
    writer.writerow(columns)
    for key in order:
        writer.writerow([merged[key].get(c, '') for c in columns])
    return out.getvalue()


def _slug(text):
    return re.sub(r'[^a-z0-9]+', '_', str(text).lower()).strip('_')


def csv_text(rows):
    out = io.StringIO()
    csv.writer(out, lineterminator='\n').writerows(rows)
    return out.getvalue()


# --- README ------------------------------------------------------------------------------------

def readme_text(title, proteins, extra_lines=(), built=False):
    """The README of a download: what it holds, the release, the columns, licence and citation.
    built=True for the catalogue files (made once, not at download)."""
    codes = ', '.join(p.protein_code for p in proteins)
    lines = [
        title,
        '=' * len(title),
        '',
        f'LPMOcat, https://lpmo-cat.polsl.pl -- {release_label() or "data release not recorded"}',
        f'{"Built" if built else "Downloaded"}: {datetime.date.today().isoformat()}',
        f'Proteins ({len(proteins)}): {codes}',
        *extra_lines,
        '',
        'Files',
        '-----',
        'proteins.csv              one row per protein: family, organism and taxonomy, structure and its',
        '                          source, UniProt/UniParc, EC numbers, substrates and regioselectivity',
        '                          (as papers report them, else from the EC numbers), active-site metal,',
        '                          brace histidines (residue numbers) and their motifs, second-sphere',
        '                          motif, shielding residue, aromatic binding residues (Phe, Trp, Tyr',
        '                          and His on the binding surface, the brace left out, with their',
        '                          loops), CBMs, lengths and masses, disulfide bonds',
        '                          and references',
        'substrate_evidence.csv    every activity claim (substrate, form, result, oxidised positions)',
        '                          with its source (CAZy, Swiss-Prot, BRENDA, literature) and papers',
        'residues_catalytic_domain.csv, residues_full_length.csv',
        '                          one row per residue: residue_number (the structure\'s own numbering),',
        '                          position (1 = first residue of the prepared structure, normally the',
        '                          N-terminal brace His), secondary structure (DSSP, 3 states), core',
        '                          strand (b1-b9), loop region, on_binding_surface, disulfide partner,',
        '                          relative SASA (%), charge (PDB2PQR, AMBER), hydrophobicity and',
        '                          flexibility (CABS-flex RMSF, GNM, B-factor or pLDDT, where available)',
        'alignment_catalytic_domain.fasta, alignment_full_length.fasta',
        '                          the structural alignment (US-align)',
        'structures/               PDB files of the structures as the catalogue shows them (aligned',
        '                          coordinates): full_length/ and catalytic_domain/',
        '',
        'How each value is computed: https://lpmo-cat.polsl.pl/docs/',
        '',
        LICENCE_TEXT,
        '',
        'Contact: lpmocat@tunnelinggroup.pl',
        '',
    ]
    return '\n'.join(lines)


# --- the comparison zip ------------------------------------------------------------------------

def comparison_zip(codes, groups=()):
    """The zip of one comparison (proteins in `codes` order): tables, per-residue data in both
    views, the alignment and the structures; `groups` (groups.with_ungrouped) add a group column
    to proteins.csv and their members to the README."""
    from . import views
    from .models import Protein

    by_code = {p.protein_code: p for p in Protein.objects.select_related(
        'clade', 'organism', 'electrostatics').filter(protein_code__in=codes)}
    proteins = [by_code[c] for c in codes if c in by_code]
    bundle = views._get_molstar_coloring([p.protein_code for p in proteins])
    scale = (bundle.get('scale_names') or ['Kyte-Doolittle'])[0]

    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, 'w', zipfile.ZIP_DEFLATED) as zf:
        group_of = {c: g['name'] for g in groups for c in g['codes']}
        group_lines = [f"Group {g['name']}: {', '.join(g['codes'])}" for g in groups]
        zf.writestr('README.txt', readme_text('LPMOcat comparison', proteins, extra_lines=group_lines))
        zf.writestr('proteins.csv', csv_text(protein_rows(proteins, group_of or None)))
        zf.writestr('substrate_evidence.csv', csv_text(substrate_rows(proteins)))
        for view, name in VIEWS:
            export = bundle['export'][view]
            zf.writestr(f'residues_{name}.csv', residue_csv(export, scale))
            zf.writestr(f'alignment_{name}.fasta', export['Sequence']['content'] + '\n')
            for filename, text in views.structure_files(proteins, catalytic=view == 'catalytic'):
                zf.writestr(f'structures/{name}/{filename}', text)
    return buffer.getvalue()
