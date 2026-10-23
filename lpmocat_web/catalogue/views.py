from django.shortcuts import render, redirect
from django.urls import reverse
from django.http import FileResponse, HttpResponse, Http404, JsonResponse
from django.template import loader
from django.conf import settings
from django.core.cache import cache, caches
from catalogue.plot_subset import assemble_from_pregen, compute_keep_indices
# from django.contrib.admin.views.decorators import staff_member_required
from requests import session
from scipy.version import short_version
from .models import Protein, Substrate, Clade, Organism, Flexibility, RegioSelecitivity
from .activities import activity_data, activity_references, table_cells, export_rows, EXPORT_HEADER
from collections import Counter, defaultdict
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import textwrap
from .sec_struct_plot import sec_struct_plot
from plotly.subplots import make_subplots
from .dendrogram import DendrogramHandler
from .overview import (REGIOSELECTIVITY_SHORT, TAXON_GROUPS, _lineage, catalogue_overview, catalogue_summary, lineage_text,
                       protein_loop_profile, protein_residue_profile, protein_residue_sets,
                       protein_binding_aromatics, protein_structural_features,
                       regioselectivity_label,
                       structure_source, taxon_group)
from .pipeline_scheme import pipeline_scheme
from . import compare_limits, downloads, groups as protein_groups
from .brace_motifs import brace_motifs, filter_options as brace_filter_options
from .second_sphere import NOT_DETERMINED, OPEN_POSITION, filter_options, second_sphere
from .cbms import family_order as cbm_family_order
from .makePlots import  makePlots, parse_pdb_residue_numbers, catalytic_masked_sequence, FLEX_METHOD_META, recolor_hole_hopping_figure
import plotly.io as pio
import traceback
import ast
import json
import time
import orjson
from django.db.models import Prefetch
import re
import io
import zipfile
from pathlib import Path as _Path

# Load catalytic domain boundaries once at module level
_domain_boundaries_path = _Path(__file__).parents[2] / 'data' / '06_catalytic_domain' / 'catalytic_domain_boundaries.json'
try:
    with open(_domain_boundaries_path) as _f:
        _domain_boundaries = json.load(_f)
except Exception:
    _domain_boundaries = {}

def _catalytic_domain_length(protein):
    """Residues of the catalytic domain present in the structure: those up to the catalytic-domain
    end in catalytic_domain_boundaries.json (the assignment used across the app), counted in the
    structure sequence, so unresolved residues don't count. None when no boundary is known."""
    if _domain_boundaries.get(protein.selected_structure or '', {}).get('end_residue') is None:
        return None
    return sum(1 for c in _catalytic_masked_seq(protein) if c.isalpha())

def activeUsersAPI(request):
    active_users = cache.get("active_users", {})
    return JsonResponse({"active_users": len(active_users)})

def _first_range(ranges):
    """(start, end) of the first 'start-end' string of a region list, or None."""
    try:
        start, end = ranges[0].split('-')
        return int(start), int(end)
    except (TypeError, ValueError, IndexError, AttributeError):
        return None


def _all_ranges(ranges):
    """[(start, end)] of every 'start-end' string of a region list. The binding surface of a
    loop is several short runs (pipeline step 7c), not one range."""
    out = []
    for piece in ranges if isinstance(ranges, (list, tuple)) else []:
        try:
            start, end = str(piece).split('-')
            out.append((int(start), int(end)))
        except (TypeError, ValueError):
            continue
    return sorted(out)


def _loop_flexibility(protein):
    """Browse-page loop data of a protein: ({loop: mean CABSflex RMSF}, {loop: 'start-end' part
    on the binding surface}, {loop: mean RMSF of that part}). Ranges are sequence positions;
    the RMSF values are stored in structure order (Flexibility.rmsf JSON)."""
    averages, overlap, overlap_averages = {}, {}, {}
    try:
        loops = ast.literal_eval(protein.loop_regions) if protein.loop_regions else {}
        surface = ast.literal_eval(protein.binding_surface) if protein.binding_surface else {}
    except (ValueError, SyntaxError):
        return averages, overlap, overlap_averages
    if not isinstance(loops, dict):
        return averages, overlap, overlap_averages
    if not isinstance(surface, dict):
        surface = {}
    values = []
    flex = next((f for f in protein.flexibility.all() if f.method == 'cabsflex'), None)
    if flex is not None:
        try:
            values = json.loads(flex.rmsf).get('values') or []
        except (TypeError, ValueError, AttributeError):
            values = []

    def mean(first, last):
        part = [v for v in values[first - 1:last] if v is not None]
        return round(sum(part) / len(part), 2) if part else None

    def mean_of(runs):
        part = [v for first, last in runs for v in values[first - 1:last] if v is not None]
        return round(sum(part) / len(part), 2) if part else None

    for group, ranges in loops.items():
        loop = _first_range(ranges)
        if loop is None:
            continue
        if mean(*loop) is not None:
            averages[group] = mean(*loop)
        # every run of the binding surface that falls inside this loop
        runs = [(max(loop[0], first), min(loop[1], last))
                for first, last in _all_ranges(surface.get(group))
                if max(loop[0], first) <= min(loop[1], last)]
        if runs:
            overlap[group] = ','.join(f"{first}-{last}" for first, last in runs)
            if mean_of(runs) is not None:
                overlap_averages[group] = mean_of(runs)
    return averages, overlap, overlap_averages


def index(request):
    # link_list removed - not used in templates

    # OPTIMIZATION: Use select_related and prefetch_related to avoid N+1 queries
    proteins = Protein.objects.select_related(
        'organism', 'clade', 'surface_accesability', 'electrostatics', 'hole_hopping'
    ).prefetch_related(
        'flexibility', 'substrate', 'regioselectivity'
    ).order_by('protein_name')
    
    DH = DendrogramHandler()
    
    dendrogram_graph = DH.static_dendrogram_circular()

    clads = list(map(lambda x: x.name, Clade.objects.all()))
    # Filter out invalid substrate names like 'nan' that can appear from imports/CSV parsing.
    substrates = sorted(Protein.objects.exclude(substrate__name__isnull=True)
        .exclude(substrate__name__exact='')
        .values_list('substrate__name', flat=True)
        .distinct())
    _all_org_names = Protein.objects.values_list('organism__name', flat=True).distinct()
    organisms = sorted({n.strip() for n in _all_org_names if n and n.strip()})
    import re as _re
    _ec_raw = Protein.objects.exclude(ec_numbers='').values_list('ec_numbers', flat=True)
    ec_options = sorted(set(
        ec.strip() for raw in _ec_raw
        for ec in _re.split(r'[;,]', raw)
        if ec.strip() and ec.strip() != '1.14.99.-'
    ))
    from django.db.models import Count as _Count
    has_no_substrate_proteins = Protein.objects.annotate(
        sub_count=_Count('substrate')
    ).filter(sub_count=0).exists()
    # Regioselectivity filter options - only show C1, C4, C1/C4
    regioselectivity_options = [(name, REGIOSELECTIVITY_SHORT[name]) for name in
                                ('C1-hydroxylating', 'C4-dehydrogenating', 'C1-hydroxylating/C4-dehydrogenating')]
    
    proteins_with_data = []
    for protein in proteins:
        # OPTIMIZATION: Use pre-calculated fields from database instead of calculating on-the-fly
        # If fields are None, they haven't been calculated yet (for old records)
        # In that case, calculate them and save (one-time migration)
        if protein.flex_min is None or protein.flex_max is None:
            try:
                flex_data = protein.flexibility.all()
                if flex_data.exists():
                    rmsf_text = flex_data.first().rmsf
                    from .datashell import calculate_flexibility_stats
                    flex_min, flex_max = calculate_flexibility_stats(rmsf_text)
                    protein.flex_min = flex_min
                    protein.flex_max = flex_max
                    protein.save(update_fields=['flex_min', 'flex_max'])
            except Exception as e:
                pass
        
        # Use pre-calculated values (or None if not available)
        protein.flex_min = protein.flex_min
        protein.flex_max = protein.flex_max
        
        try:
            if protein.binding_surface:
                binding_dict = eval(protein.binding_surface)
                protein.binding_surface_json = json.dumps(binding_dict)
            else:
                protein.binding_surface_json = '{}'
        except Exception as e:
            protein.binding_surface_json = '{}'

        try:
            if getattr(protein, 'loop_regions', None):
                loop_dict = eval(protein.loop_regions)
                if isinstance(loop_dict, dict):
                    protein.loop_regions_json = json.dumps(loop_dict)
                else:
                    protein.loop_regions_json = '{}'
            else:
                protein.loop_regions_json = '{}'
        except Exception:
            protein.loop_regions_json = '{}'
        
        # Average RMSF of each loop region, and the part of each loop on the binding surface
        loop_rmsf_averages, loop_binding_overlap, loop_binding_rmsf_averages = \
            _loop_flexibility(protein)
        protein.loop_rmsf_averages_json = json.dumps(loop_rmsf_averages)
        protein.loop_binding_overlap_json = json.dumps(loop_binding_overlap)
        protein.loop_binding_rmsf_averages_json = json.dumps(loop_binding_rmsf_averages)
        
        # First, extract axial position aromatic (from second_sphere_aromatic field)
        axial_residue = ''
        try:
            if protein.second_sphere_aromatic and str(protein.second_sphere_aromatic).lower() not in ['nan', 'none', 'null', '']:
                # The field stores tuples like ('PHE', 123) or just the residue name
                axial_value = str(protein.second_sphere_aromatic)
                if 'PHE' in axial_value.upper():
                    axial_residue = 'PHE'
                elif 'TYR' in axial_value.upper():
                    axial_residue = 'TYR'
        except Exception:
            pass
        
        protein.axial_residue = axial_residue
        
        # Extract and categorize second-sphere residues from histidine brace interactions
        # Hierarchy: Axial > Pi-Pi Stacking > Cation-Pi > Hydrogen Bonds
        # A residue is only counted in the highest priority category it belongs to
        hbond_residues = set()
        pi_stacking_residues = set()
        cation_pi_residues = set()
        
        try:
            interactions_raw = getattr(protein, 'histidine_brace_interactions', None)
            if interactions_raw and str(interactions_raw).lower() not in ['nan', 'none', 'null', '']:
                interactions_json = json.loads(interactions_raw)
                
                # First, collect pi-pi stacking residues (exclude axial residue)
                for pi in (interactions_json.get('pi_pi_stacking') or []):
                    res_name = pi.get('other_res_name')
                    # Only add if it's not the axial residue
                    if res_name and res_name != axial_residue:
                        pi_stacking_residues.add(res_name)
                
                # Then, collect cation-pi residues (exclude axial residue AND pi-stacking residues)
                for cp in (interactions_json.get('cation_pi') or []):
                    res_name = cp.get('other_res_name')
                    # Only add if it's not the axial residue and not already in pi-stacking
                    if res_name and res_name != axial_residue and res_name not in pi_stacking_residues:
                        cation_pi_residues.add(res_name)
                
                # Finally, collect hydrogen bond residues (exclude axial, pi-stacking, AND cation-pi residues)
                for hb in (interactions_json.get('hydrogen_bonds') or []):
                    res_name = hb.get('other_res_name')
                    # Only add if it's not the axial residue and not already in pi-stacking or cation-pi
                    if res_name and res_name != axial_residue and res_name not in pi_stacking_residues and res_name not in cation_pi_residues:
                        hbond_residues.add(res_name)
        except (json.JSONDecodeError, TypeError, ValueError, AttributeError):
            pass
        
        protein.hbond_residues = ','.join(sorted(hbond_residues))
        protein.pi_stacking_residues = ','.join(sorted(pi_stacking_residues))
        protein.cation_pi_residues = ','.join(sorted(cation_pi_residues))

        # Second-sphere motif and the shielding residue on the copper's solvent side (8c-8d)
        site = second_sphere(protein)
        protein.motif_short = site['short']
        protein.motif_full = site['motif']
        protein.distal_axial = site['distal']
        # Motifs after the N-terminal and the internal His of the brace (8e)
        brace = brace_motifs(protein)
        protein.n_his_motif = brace['n_his']['motif']
        protein.internal_his_motif = brace['internal_his']['motif']
        
        # OPTIMIZATION: Use pre-calculated sequence_length from database
        # If field is None, calculate it and save (one-time migration)
        if protein.sequence_length is None:
            try:
                seq_to_use = None
                if hasattr(protein, 'aligned_sequence') and protein.aligned_sequence:
                    seq_to_use = protein.aligned_sequence
                elif protein.sequence:
                    seq_to_use = protein.sequence
                
                from .datashell import calculate_sequence_length
                if seq_to_use:
                    protein.sequence_length = calculate_sequence_length(seq_to_use)
                else:
                    protein.sequence_length = 0
                protein.save(update_fields=['sequence_length'])
            except Exception:
                protein.sequence_length = 0
        
        # Calculate total charge from PQR file
        total_charge = None
        if hasattr(protein, 'electrostatics') and protein.electrostatics and protein.electrostatics.pqr_file:
            match = re.search(r"Total charge.*?:\s*([+-]?\d+\.?\d*)", protein.electrostatics.pqr_file)
            if match:
                try:
                    total_charge = float(match.group(1))
                except (ValueError, TypeError):
                    total_charge = None
        
        protein.total_charge = total_charge

        # the overview panel's values; the catalytic-domain charge is also the charge filter's
        _set_overview_values(protein)
        # filter values named as the overview names them: the taxon group of the source organism
        # (not CAZy's kingdom), the structure source and the disulfide bonds of the catalytic domain
        protein.taxon_group_name = taxon_group(_lineage(protein.organism), protein.organism.kingdom)[0]
        protein.structure_source_name = structure_source(protein.selected_structure)
        protein.disulfide_count = (protein.structural_features or {}).get('disulfides')
        # the aromatic binding residues for the criteria: their one-letter codes and their number
        aromatics = (protein.structural_features or {}).get('binding_aromatic_residues')
        protein.binding_aromatic_types = ''.join(sorted({a['aa'] for a in aromatics})) if aromatics is not None else None
        protein.binding_aromatic_count = len(aromatics) if aromatics is not None else None

        proteins_with_data.append(protein)
    
    # Collect unique residues by category
    hbond_residues_all = set()
    pi_stacking_residues_all = set()
    cation_pi_residues_all = set()
    
    for protein in proteins_with_data:
        if hasattr(protein, 'hbond_residues') and protein.hbond_residues:
            residues = protein.hbond_residues.split(',')
            hbond_residues_all.update(residues)
        
        if hasattr(protein, 'pi_stacking_residues') and protein.pi_stacking_residues:
            residues = protein.pi_stacking_residues.split(',')
            pi_stacking_residues_all.update(residues)
        
        if hasattr(protein, 'cation_pi_residues') and protein.cation_pi_residues:
            residues = protein.cation_pi_residues.split(',')
            cation_pi_residues_all.update(residues)
    
    # Motif filters: the short motifs most common first, the spaced ones grouped under them
    motif_options = filter_options(p.motif_short for p in proteins_with_data)
    motif_rank = {short: rank for rank, (short, _) in enumerate(motif_options)}
    short_of = {p.motif_full: p.motif_short for p in proteins_with_data if p.motif_full}
    motif_full_options = sorted(filter_options(p.motif_full for p in proteins_with_data),
                                key=lambda item: (motif_rank.get(short_of[item[0]], 0), -item[1]))
    distal_axial_options = filter_options(p.distal_axial for p in proteins_with_data)
    n_his_motif_options = brace_filter_options(p.n_his_motif for p in proteins_with_data)
    internal_his_motif_options = brace_filter_options(p.internal_his_motif for p in proteins_with_data)
    hbond_residues_list = sorted(list(hbond_residues_all))
    pi_stacking_residues_list = sorted(list(pi_stacking_residues_all))
    cation_pi_residues_list = sorted(list(cation_pi_residues_all))
    
    # Calculate min and max sequence lengths for slider initialization
    sequence_lengths = [p.sequence_length for p in proteins_with_data if hasattr(p, 'sequence_length') and p.sequence_length > 0]
    min_sequence_length = min(sequence_lengths) if sequence_lengths else 0
    max_sequence_length = max(sequence_lengths) if sequence_lengths else 500

    # Catalytic-domain length range (for the "catalytic domain only" slider mode)
    cat_lengths = [p.cat_domain_length for p in proteins_with_data if getattr(p, 'cat_domain_length', None)]
    min_cat_length = min(cat_lengths) if cat_lengths else min_sequence_length
    max_cat_length = max(cat_lengths) if cat_lengths else max_sequence_length

    # Mass ranges (kDa) for the same slider in "Mass (kDa)" mode, as plain "12.3" strings
    sequence_masses = [p.sequence_mass for p in proteins_with_data if p.sequence_mass]
    cat_masses = [p.cat_domain_mass for p in proteins_with_data if p.cat_domain_mass]
    min_sequence_mass = f'{min(sequence_masses):.1f}' if sequence_masses else '0.0'
    max_sequence_mass = f'{max(sequence_masses):.1f}' if sequence_masses else '100.0'
    min_cat_mass = f'{min(cat_masses):.1f}' if cat_masses else min_sequence_mass
    max_cat_mass = f'{max(cat_masses):.1f}' if cat_masses else max_sequence_mass
    
    # Extract unique CBM values from associated_CBMs field
    cbm_values = set()
    for protein in proteins_with_data:
        if protein.associated_CBMs:
            cbm_str = str(protein.associated_CBMs).strip()
            if cbm_str and cbm_str.lower() not in ['nan', 'none', 'null', '']:
                # Split by comma and clean up
                cbms = [c.strip() for c in cbm_str.split(',') if c.strip()]
                cbm_values.update(cbms)
    cbm_options = sorted(cbm_values, key=cbm_family_order)

    _metal_display = {'CU': 'Copper', 'ZN': 'Zinc', 'NI': 'Nickel', 'MN': 'Manganese', 'FE': 'Iron', 'CO': 'Cobalt', 'MG': 'Magnesium', 'MO': 'Molybdenum'}
    _metal_raw = sorted(Protein.objects.exclude(metal__isnull=True).exclude(metal__exact='').values_list('metal', flat=True).distinct())
    metal_options = [(m, _metal_display.get(m, m)) for m in _metal_raw]

    # Calculate min and max total charge for slider initialization
    total_charges = [p.total_charge for p in proteins_with_data if hasattr(p, 'total_charge') and p.total_charge is not None]
    min_total_charge = min(total_charges) if total_charges else -20
    max_total_charge = max(total_charges) if total_charges else 20

    cat_charges = [p.cat_domain_charge for p in proteins_with_data if hasattr(p, 'cat_domain_charge') and p.cat_domain_charge is not None]
    min_cat_charge = min(cat_charges) if cat_charges else -20
    max_cat_charge = max(cat_charges) if cat_charges else 20

    # "Catalogue at a glance" panel above the search bar
    # the Selection criteria's options for those three, in the overview's order ("kingdoms" keeps
    # its name in the template and the page's scripts, but holds the overview's taxon groups)
    _taxa = Counter(p.taxon_group_name for p in proteins_with_data if p.taxon_group_name)
    kingdoms = [g for g in TAXON_GROUPS if g in _taxa] + sorted(set(_taxa) - set(TAXON_GROUPS))
    _sources = {p.structure_source_name for p in proteins_with_data}
    structure_options = [(value.lower(), label) for value, label in (
        ('PDB', 'Experimental structure (PDB)'), ('AlphaFold DB', 'AlphaFold DB model'),
        ('CAZyme3D', 'CAZyme3D model'), ('own model', 'Own AlphaFold 3 model')) if value in _sources]
    _bonds = Counter(p.disulfide_count for p in proteins_with_data if p.disulfide_count is not None)
    disulfide_options = sorted(_bonds.items())
    aromatic_count_options = sorted({p.binding_aromatic_count for p in proteins_with_data
                                     if p.binding_aromatic_count is not None})
    overview = catalogue_overview(proteins_with_data, _catalytic_masked_seq)

    # Proteins to select on load: those of "Back to selection" on the compare page (?select=),
    # else, once, those of the last comparison (the browser's back button)
    requested = request.GET.get('select') or ''
    remembered = request.session.pop('selected_protein_ids', None)
    known = {p.protein_code for p in proteins_with_data}
    sessionCookie = [code for code in (requested.split(',') if requested else remembered or [])
                     if code in known] or None
    # the groups of "Back to selection" (?groups=), restored by the page's script
    groups_on_load = protein_groups.parse(request.GET.get('groups'), sorted(known))

    if request.method == "POST":
        # "Compare selected": the ticked proteins, in the order they were sent (a set's order
        # changes from one server process to the next, and the compare page's cache is keyed by
        # the order); a whole family in its fixed order
        selected = compare_limits.unique(request.POST.getlist('protein_ids'))
        # groups ("Add as group"): the comparison in group order, the groups in its address
        named = protein_groups.parse(request.POST.get('groups'), selected)
        selected = (protein_groups.ordered_codes(named, selected) if named
                    else compare_limits.canonical_order(selected))
        if compare_limits.too_many(selected):
            return _back_to_selection(selected)
        if selected:
            request.session['selected_protein_ids'] = selected
            url = reverse('compare_proteins', kwargs={'protein_ids': ','.join(selected)})
            if named:
                from urllib.parse import urlencode
                url += '?' + urlencode({'groups': protein_groups.to_param(named)})
            return redirect(url)

    context = {
        'proteins': proteins_with_data,
        'clads': clads,
        'substrates': substrates,
        'has_no_substrate_proteins': has_no_substrate_proteins,
        'organisms': organisms,
        'kingdoms': kingdoms,
        'max_compare': compare_limits.MAX_COMPARE_PROTEINS,
        'disulfide_options': disulfide_options,
        'aromatic_count_options': aromatic_count_options,
        'ec_options': ec_options,
        'structure_options': structure_options,
        'regioselectivity_options': regioselectivity_options,
        'cbm_options': cbm_options,
        'metal_options': metal_options,
        'min_total_charge': min_total_charge,
        'max_total_charge': max_total_charge,
        'min_cat_charge': min_cat_charge,
        'max_cat_charge': max_cat_charge,
        'hbond_residues': hbond_residues_list,
        'pi_stacking_residues': pi_stacking_residues_list,
        'cation_pi_residues': cation_pi_residues_list,
        'motif_options': motif_options,
        'motif_full_options': motif_full_options,
        'distal_axial_options': distal_axial_options,
        'n_his_motif_options': n_his_motif_options,
        'internal_his_motif_options': internal_his_motif_options,
        'min_sequence_length': min_sequence_length,
        'max_sequence_length': max_sequence_length,
        'min_cat_length': min_cat_length,
        'max_cat_length': max_cat_length,
        'min_sequence_mass': min_sequence_mass,
        'max_sequence_mass': max_sequence_mass,
        'min_cat_mass': min_cat_mass,
        'max_cat_mass': max_cat_mass,
        'dendrogram_graph': dendrogram_graph,
        'overview': overview,
        'sessionCookie': sessionCookie,
        'groups_on_load': [{'name': g['name'], 'codes': g['codes']} for g in groups_on_load],
        'max_groups': protein_groups.MAX_GROUPS,
    }
    return render(request, 'catalogue/index.html', context)


# When the pipeline data were last retrieved and loaded (update after a reload)
CATALOGUE_DATA_UPDATED = 'September 2026'


def _catalogue_summary():
    """The About and How to Use pages' counts (catalogue_summary)."""
    proteins = list(Protein.objects.select_related('clade', 'organism')
                    .only('selected_structure', 'clade__name', 'organism__name',
                          'organism__kingdom', 'organism__lineage'))
    return dict(catalogue_summary(proteins), updated=CATALOGUE_DATA_UPDATED)


def home(request):
    """Home/About page: what the catalogue holds, how it is built and where its data come from"""
    files = _download_files()
    return render(request, "catalogue/home.html",
                  {"summary": _catalogue_summary(), "scheme": pipeline_scheme(),
                   "docs_url": docs_url(), "download_files": files,
                   "downloads_built": any(f['size'] for f in files),
                   "release": downloads.release_info()})


def docs_url():
    """URL of the pipeline documentation (docs/, built by MkDocs and served by WhiteNoise, see
    settings.py), or None when it has not been built."""
    return settings.DOCS_URL if settings.DOCS_BUILT else None






def _reference_url(info):
    """Best link for a citation: DOI, else PubMed, else None."""
    doi = info.get('doi')
    if doi:
        return f'https://doi.org/{doi}'
    pmid = info.get('pubmed') or info.get('pmid')
    if pmid:
        return f'https://pubmed.ncbi.nlm.nih.gov/{pmid}/'
    return None


def _citation_keys(info):
    """Identity keys for a citation, used to dedupe across sources (DOI / PubMed)."""
    keys = set()
    doi = info.get('doi')
    if doi:
        keys.add(('doi', str(doi).strip().lower()))
    pmid = info.get('pubmed') or info.get('pmid')
    if pmid:
        keys.add(('pmid', str(pmid)))
    return keys


def _format_references(protein, extra=None):
    """Numbered 'Author et al., Year' citation links for the comparison table.

    Merges three sources, deduped by DOI/PubMed so a paper is shown once:
      1. Structure papers from RCSB (`other_crystal_references`, keyed by PDB,
         selected structure first) — populated by `backfill_references`.
      2. Characterization papers from CAZy (`characterization_references`, a PMID
         list) — populated by `backfill_cazy_references`.
      3. Papers the protein's UniProt entries cite for function, activity or structure
         (`uniprot_references`, pipeline step 2f) — loaded by `populate_db`.
    2 and 3 are the only sources for AlphaFold-only proteins.
    Falls back to the legacy {pdb: doi} / bare crystal_reference DOI formats.

    `extra` adds the papers of the activity record (step 2j, `activities.activity_references`)
    that none of the three sources lists; each carries the "key" the record cites it by.

    Returns (html, {citation key: {"n", "label"}}): what the EC, substrate and oxidised-position
    cells of the same row cite, keyed by ('pmid', id) / ('doi', doi) and by an activity ref.
    """
    import html as _html
    import json as _json

    # Match the green link style used by the other table columns (UniProt / PDB).
    link_style = 'color:#569976;text-decoration:underline;'

    def _render(info, fallback_label):
        author, year = info.get('author'), info.get('year')
        etal = ' et al.' if (info.get('n_authors') or 0) > 1 else ''
        if author and year:
            label = f"{author}{etal}, {year}"
        elif author:
            label = f"{author}{etal}"
        elif year:
            label = str(year)
        else:
            label = fallback_label
        labels.append(str(label))
        label = _html.escape(str(label))
        tip = _html.escape(str(info.get('title') or fallback_label), quote=True)
        url = _reference_url(info)
        if url:
            return (f'<a href="{_html.escape(url, quote=True)}" target="_blank" rel="noopener" '
                    f'style="{link_style}" title="{tip}">{label}</a>')
        return f'<span title="{tip}">{label}</span>'

    parts = []
    labels = []
    seen = set()
    index = {}

    def _index(keys, extra_key=None):
        """Number the entry just appended and index it by its citation keys."""
        entry = {'n': len(parts), 'label': labels[-1]}
        for key in keys | ({extra_key} if extra_key else set()):
            index.setdefault(key, entry)

    # 1) Structure papers (selected first, deduped by shared paper).
    raw = getattr(protein, 'other_crystal_references', None)
    if raw:
        try:
            refs = _json.loads(raw)
        except Exception:
            refs = {}
        def _sort_key(item):
            pdb, info = item
            return (0 if (isinstance(info, dict) and info.get('selected')) else 1, pdb)
        for pdb, info in sorted(refs.items(), key=_sort_key):
            if isinstance(info, str):                # legacy {pdb: doi}
                info = {'doi': info}
            if not isinstance(info, dict):
                continue
            # Skip unpublished structures ("To Be Published" — no DOI and no PubMed):
            # they have no citable paper, so don't show a link-less placeholder.
            if not info.get('doi') and not info.get('pubmed'):
                continue
            keys = _citation_keys(info)
            if keys & seen:
                continue
            seen |= keys
            parts.append(_render(info, pdb))
            _index(keys)
    else:
        crystal_ref = getattr(protein, 'crystal_reference', None)
        if crystal_ref:
            info = {'doi': crystal_ref}
            keys = _citation_keys(info)
            seen |= keys
            parts.append(_render(info, getattr(protein, 'selected_structure', '') or ''))
            _index(keys)

    # 2) Characterization papers from CAZy, then 3) papers cited by the UniProt entries
    # (skip any already shown).
    for field in ('characterization_references', 'uniprot_references'):
        craw = getattr(protein, field, None)
        if not craw:
            continue
        try:
            clist = _json.loads(craw)
        except Exception:
            clist = []
        for info in clist:
            if not isinstance(info, dict):
                continue
            keys = _citation_keys(info)
            if keys and (keys & seen):
                continue
            seen |= keys
            parts.append(_render(info, info.get('pmid') or 'ref'))
            _index(keys)

    # 4) Papers only the activity record (step 2j) cites, e.g. a BRENDA-only paper.
    for info in extra or []:
        keys = _citation_keys(info)
        if keys and (keys & seen):
            index.setdefault(info['key'], index.get(next(iter(keys & seen))))  # same paper
            continue
        seen |= keys
        parts.append(_render(info, info.get('title') or info['key']))
        _index(keys, info['key'])

    numbered = [f'<span class="ref-entry"><span class="ref-n">{i}.</span>&nbsp;{part}</span>'
                for i, part in enumerate(parts, 1)]
    return ('; '.join(numbered) if numbered else ''), index


_UNIPROT_SOURCE_TEXT = {
    'CAZy': 'listed by CAZy',
    'UniParc': 'same sequence as the CAZy GenBank protein',
    'CAZyme3D model': 'entry of the CAZyme3D model',
    'identical sequence': 'identical sequence',
    # entries of the other lists of pipeline step 1 (not in CAZy's lists)
    'Swiss-Prot': 'Swiss-Prot entry with experimentally shown activity',
    'BRENDA': 'BRENDA protein with a cited substrate',
    'curated': 'unreviewed entry with published experimental results',
    'strain variant': 'the same protein from another strain (>= 95 % identical)',
}


def _uniprot_notes(protein):
    """Hover text per active UniProt accession -- review status, how the entry belongs to the
    protein, organism -- as JSON {accession: text}, from Protein.uniprot_all."""
    import json as _json
    try:
        accessions = _json.loads(getattr(protein, 'uniprot_all', None) or '{}')
    except (ValueError, TypeError):
        return '{}'
    notes = {}
    for acc, info in accessions.items():
        if not isinstance(info, dict) or info.get('status') == 'deleted':
            continue
        parts = [{'reviewed': 'Swiss-Prot', 'unreviewed': 'TrEMBL'}.get(info.get('status'), ''),
                 _UNIPROT_SOURCE_TEXT.get(info.get('source'), info.get('source') or ''),
                 info.get('organism') or '']
        notes[acc] = ' · '.join(p for p in parts if p)
    return _json.dumps(notes)


# Version tag baked into the plot-bundle cache key. Bump this whenever the plot
# generation logic changes (colors, layout, new plot types, alignment trimming)
# so stale bundles are ignored instead of served.
COMPARE_PLOTS_CACHE_VERSION = 'v41'  # v41: longer binding-surface ticks for the aromatic residues; v40: final reload of 2026-10-21 (CsAA15 added, CAZy read again); v39: binding-surface ticks and disulfide dots act on the whole set of the protein (hover/click texts); v38: hole hopping donors from the whole-residue relative SASA of step 5 (>= 10 %, no per-structure thresholds); v37: binding-surface ticks in the Secondary Structure and Loops plot, every run of a region kept (the plane gives several per loop); v36: ScAA9A's own AF3 model (stale UniProt gene model), C-terminal expression tags removed (3i), structure selection prefers the wild type and the copper (4A02->4ALS, 5LSV->4OPB, 8GUM->8GUL, 4OW5->4X27), binding surface from the tangent plane at the copper face; v35: 'Shielding residue' (was solvent-facing axial position) in the Sequence plot and its 3D category; v34: CBMs of other UniProt entries placed on the structure by alignment (domain bar of TtAA9B, NcAA13, VhAA10A); v33: the Swiss-Prot/BRENDA/curated lists and the own AlphaFold 3 models (new proteins, so a new alignment); v32: brace histidine motifs (pipeline 8e) in the Sequence plot; v31: MtAA9 His161 in its H-X5-Q-X-Y motif; v30: second-sphere motif and solvent-facing axial residue in the Sequence plot; v29: disulfide-bonded Cys dots in the Secondary Structure and Loops plot, catalytic ends keep disulfide bonds (step 6e); v28: relative SASA from pipeline step 5a (was A^2 read as %); v27: core strands without isolated bridges; v26: core strand pieces, split b6 (LS after its last piece); v25: core strands b1-b9 (AA9 hairpin = b2), loops from the core; v24: Secondary Structure and Loops plot (core strands); v23: DSSP (mkdssp) secondary structure; v22: one catalytic-domain assignment (boundaries file) for bars and lengths; Sequence pregen rebuilt


def _compare_plots_cache_key(protein_records):
    """Cache key for a compare-page plot bundle.

    Plots are trimmed of all-gap alignment columns *across the selected set*, so
    the content depends on the exact combination. The per-protein arrays are also
    ordered to match the page (protein_records / URL order), and each plot's
    content depends on its protein's index — so the key must preserve order, not
    sort. Two different orderings of the same set are cached separately.
    """
    codes = [p.protein_code for p in protein_records]
    return 'compareplots:{}:{}'.format(COMPARE_PLOTS_CACHE_VERSION, ','.join(codes))


# Plot-type order and legend/description HTML for the compare page section
# skeletons (rendered without generating any plots).
_COMPARE_PLOT_TYPES = ['Sequence', 'Hydrophobicity', 'Charge', 'Secondary Structure and Loops',
                       'Relative Surface Accessibility', 'Flexibility', 'Hole Hopping']


# The section titles where they differ from the plot type (the type stays the key of the
# caches, the pregenerated files and the JS).
_COMPARE_PLOT_LABELS = {
    'Secondary Structure and Loops': 'Secondary Structure, Loop Regions, Binding Surface',
}

_COMPARE_PLOT_DESCRIPTIONS = {
    'Sequence': """
        <p>Amino acid sequence colored by residue properties</p>
        <div class='legend-item'><span class='legend-box' style='background:rgb(186, 3, 186);'></span>Histidine brace coordinating copper ion</div>
        <div class='legend-item'><span class='legend-box' style='background:rgb(255, 140, 0);'></span>Second-sphere: Axial position (to Cu)</div>
        <div class='legend-item'><span class='legend-box' style='background:rgb(38, 166, 154);'></span>Second-sphere motif: the Gln/Glu near the copper and the residue holding it (e.g. H-X8-Q-X-Y, E…N-X-F; the axial aromatic ends the motif)</div>
        <div class='legend-item'><span class='legend-box' style='background:rgb(141, 110, 99);'></span>Second-sphere: shielding residue (over the copper on its solvent side, opposite the aromatic; the Ala of AA10s)</div>
        <div class='legend-item'><span class='legend-box' style='background:rgb(68, 136, 252);'></span>Brace histidine motifs: the Gly, Pro and Arg after the N-terminal His (H-G-X4-P-X2-R, H-G-X4-P, H-G) and the Gly and Pro after the internal His (H-X-G-P, H-X2-G-P, H-X-G)</div>
        <div class='legend-item'><span class='legend-box' style='background:rgb(255, 165, 0);'></span>Second-sphere: Pi-pi stacking</div>
        <div class='legend-item'><span class='legend-box' style='background:rgb(255, 175, 50);'></span>Second-sphere: Cation-pi interaction</div>
        <div class='legend-item'><span class='legend-box' style='background:rgb(255, 200, 100);'></span>Second-sphere: H-bond</div>
        <div class='legend-item'><span class='legend-box' style='background:rgb(255, 180, 180);'></span>Positively charged residue</div>
        <div class='legend-item'><span class='legend-box' style='background:rgb(180, 180, 255);'></span>Negatively charged residue</div>
        <div class='legend-item'><span class='legend-box' style='background:rgb(180, 255, 180);'></span>Aromatic residue</div>
        <div class='legend-item'><span class='legend-box' style='background:rgb(255, 255, 180);'></span>Residue containing sulphur</div>

    """,
    'Flexibility': """
        <p>Per-residue flexibility. Use the <b>Method</b> selector to switch measure — each has its
        own scale/units and the plot is empty where a protein lacks that measure.</p>
        <div class='legend-item'><b>CABSflex</b> — RMSF [Å] from coarse-grained simulation</div>
        <div class='legend-item'><b>GNM</b> — relative mobility from the elastic network (dimensionless; higher = more mobile)</div>
        <div class='legend-item'><b>B-factor</b> — experimental crystallographic B-factors [Å²]</div>
        <div class='legend-item'><b>pLDDT</b> — AlphaFold per-residue confidence (0–100)</div>
        <p style='margin-top:6px;'>Colouring:</p>
        <div class='legend-item'><span class='legend-box-flexibility' style='background:linear-gradient(0.25turn, #ffc8c8, #ff7070, #ff0000);'></span>Rigid → flexible (CABSflex / GNM / B-factor, scaled per method)</div>
        <div class='legend-item'><span class='legend-box-flexibility' style='background:linear-gradient(0.25turn, #FF7D45, #FFDB13, #65CBF3, #0053D6);'></span>pLDDT: low → high confidence (AlphaFold colours; low = flexible)</div>
    """,
    'Hydrophobicity': """
        <p>Residue hydrophobicity. Use the scale selector to switch between different scales.</p>
        <div class='legend-item'><span class='legend-box-hydrophobicity' style='background:linear-gradient(0.25turn, #420dbd, #bcabf5, #f59d9d, #f72d2d);'></span>Hydrophilic - Hydrophobic</div>
        <table style='font-size:0.85em; margin-top:6px; border-collapse:collapse;'>
        <tr><th style='text-align:left; padding:2px 8px;'>Scale</th><th style='padding:2px 8px;'>Range</th></tr>
        <tr><td style='padding:2px 8px;'><a href='https://doi.org/10.1016/0022-2836(82)90515-0' target='_blank'>Kyte&ndash;Doolittle</a></td><td style='padding:2px 8px;'>&minus;4.5 to 4.5</td></tr>
        <tr><td style='padding:2px 8px;'><a href='https://doi.org/10.1016/0161-5890(83)90029-9' target='_blank'>Hopp&ndash;Woods</a></td><td style='padding:2px 8px;'>&minus;3.4 to 3.0</td></tr>
        <tr><td style='padding:2px 8px;'><a href='https://doi.org/10.1016/0022-2836(87)90189-6' target='_blank'>Cornette</a></td><td style='padding:2px 8px;'>&minus;3.1 to 5.7</td></tr>
        <tr><td style='padding:2px 8px;'><a href='https://doi.org/10.1016/0022-2836(84)90309-7' target='_blank'>Eisenberg</a></td><td style='padding:2px 8px;'>&minus;2.53 to 1.38</td></tr>
        <tr><td style='padding:2px 8px;'><a href='https://doi.org/10.1126/science.4023714' target='_blank'>Rose</a></td><td style='padding:2px 8px;'>0.52 to 0.91</td></tr>
        <tr><td style='padding:2px 8px;'><a href='https://doi.org/10.1038/277491a0' target='_blank'>Janin</a></td><td style='padding:2px 8px;'>&minus;1.8 to 0.9</td></tr>
        <tr><td style='padding:2px 8px;'><a href='https://doi.org/10.1146/annurev.bb.15.060186.001541' target='_blank'>Engelman (GES)</a></td><td style='padding:2px 8px;'>&minus;12.3 to 3.7</td></tr>
        </table>
    """,
    'Charge': """
        <p>Electrostatic charge per residue calculated from PQR files (atomic partial charges). Total protein charge is shown in the PQR file header (note: copper ions are not included in the total charge calculation).</p>
        <div class='legend-item'><span class='legend-box-charge' style='background:linear-gradient(0.25turn, #ff0000, #ffffff, #0000ff);'></span>Negative (Red) - Neutral (White) - Positive (Blue)</div>
    """,
    'Secondary Structure and Loops': """
        <p>Secondary structure along the sequence (DSSP), with the beta-sandwich core strands and the
        loop regions in the bars above it, and the binding-surface residues and the disulfide-bonded
        cysteines below it.</p>
        <div class='legend-item'><span class='legend-box' style='background:#C62828;'></span>Helix (alpha, 3<sub>10</sub>, pi)</div>
        <div class='legend-item'><span class='legend-box' style='background:#FBC02D;'></span>Beta strand (incl. isolated bridges)</div>
        <div class='legend-item'><span class='legend-box' style='background:#616161;'></span>Coil (turns, bends, loops)</div>
        <p style='margin-top:8px;'><strong>Core strands</strong> (lower bar): the strands of the beta-sandwich core
        common to the LPMO families, numbered β1–β9 along the sequence (β2 is the hairpin partner of β1, found
        in AA9); the same number is the same strand in every protein, and a protein lacking a strand skips
        its number. Other strands are not named.</p>
        <p style='margin-top:8px;'><strong>Loop regions</strong> (upper bar), the segments between core strands:
        L2 β2 (or β1)–β3, L3 β4–β5, LS β6–β7, L8 β8–β9, LC after β9 (when at least 8 residues long).</p>
        <div style='display: grid; grid-template-columns: repeat(2, 1fr); gap: 5px;'>
            <div class='legend-item'><span class='legend-box' style='background:rgb(255, 69, 0);'></span>L2 (loop 2)</div>
            <div class='legend-item'><span class='legend-box' style='background:rgb(255, 140, 0);'></span>L3 (loop 3)</div>
            <div class='legend-item'><span class='legend-box' style='background:rgb(138, 43, 226);'></span>LS (loop short)</div>
            <div class='legend-item'><span class='legend-box' style='background:rgb(30, 144, 255);'></span>L8 (loop 8)</div>
            <div class='legend-item'><span class='legend-box' style='background:rgb(34, 139, 34);'></span>LC (C-terminal loop)</div>
        </div>
        <div class='legend-item' style='margin-top:8px;'><span class='legend-box' style='background:rgb(0, 121, 107);'></span><span><strong>Binding surface</strong> (ticks below the row): the residues a plane laid on the copper face touches &mdash; the stand-in for a flat crystalline substrate. Longer ticks mark the aromatic residues (Phe, Trp, Tyr, His; the brace histidines left out). Hover a tick to name the residue and the loop it belongs to and to highlight all of the protein's binding-surface residues; click to show them all as sticks.</span></div>
        <div class='legend-item' style='margin-top:8px;'><span class='legend-box' style='background:rgb(184, 134, 11);'></span><span><strong>Disulfide bonds</strong> (dots below): the cysteines bonded to another cysteine (sulfur atoms within 3.0 Å). Hover a dot to see its partner and to highlight all of the protein's bonded cysteines; click to show them all as sticks.</span></div>
        <p style='margin-top:8px;'>Hover a core strand or loop to highlight it in 3D; click it to show it as sticks.</p>
    """,
    'Relative Surface Accessibility': """
        <p>Relative Surface Accessibility calculated with BioPython. Thresholded at 25%.</p>
        <div class='legend-item'><span class='legend-box' style='background:rgba(208, 100, 88, 0.7);'></span>(> 25% - Exposed Residue)</div>
        <div class='legend-item'><span class='legend-box' style='background:rgba(95, 140, 200, 0.7);'></span>(≤ 25% - Buried Residue)</div>
    """,
    'Binding Surface': """
        <p>Compartments building the substrate-binding surface of the protein</p>
        <div style='display: grid; grid-template-columns: repeat(4, 1fr); gap: 5px;'>
            <div class='legend-item'><span class='legend-box' style='background:rgb(220, 20, 60);'></span>L2a</div>
            <div class='legend-item'><span class='legend-box' style='background:rgb(255, 69, 0);'></span>L2</div>
            <div class='legend-item'><span class='legend-box' style='background:rgb(178, 34, 34);'></span>L2b</div>
            <div class='legend-item'><span class='legend-box' style='background:rgb(255, 140, 0);'></span>L3</div>
            <div class='legend-item'><span class='legend-box' style='background:rgb(138, 43, 226);'></span>LS</div>
            <div class='legend-item'><span class='legend-box' style='background:rgb(30, 144, 255);'></span>L8</div>
            <div class='legend-item'><span class='legend-box' style='background:rgb(34, 139, 34);'></span>LC</div>
            <div class='legend-item'><span class='legend-box' style='background:rgb(50, 205, 50);'></span>LCb</div>
        </div>
    """,
    'Hole Hopping': """
        <p>Residues potentially involved in electron transfer (hole hopping) pathways. Each protein has its own colour, the same in its plot row, the 3D viewer and the path overlay (shown here for the first protein).</p>
        <div class='legend-item'><span class='legend-box' style='background:rgb(70, 119, 173);'></span>Darker shade - Donor residue (surface-exposed electron entry point)</div>
        <div class='legend-item'><span class='legend-box' style='background:rgb(153, 197, 247);'></span>Lighter shade - Bridging residue (electron relay to Cu)</div>
    """

}


# Plot types served from pre-generated full-alignment figures (combination-independent,
# column-subset per selection). The 3 line-gradient types are regenerated live.
_PREGEN_MARKER_TYPES = {'Sequence', 'Charge', 'Hole Hopping', 'Hydrophobicity'}


def _section_cache_key(codes, suffix):
    """Per-combination, per-piece cache key (order-preserving, version-tagged)."""
    return 'compareplots:{}:{}:{}'.format(COMPARE_PLOTS_CACHE_VERSION, ','.join(codes), suffix)


def _plots_cache(codes):
    """The cache of a comparison: whole families in their fixed order have their own, which is
    never culled (see settings.CACHES and manage.py warm_compare_cache)."""
    codes = list(codes)
    if compare_limits.whole_family(codes) and codes == sorted(codes):
        return caches['family_plots']
    return caches['plots']


def _back_to_selection(codes):
    """To the browse page with the proteins still selected; the page then says, under its toolbar,
    how many proteins a comparison holds."""
    from urllib.parse import urlencode
    return redirect(reverse('index') + '?' + urlencode({'select': ','.join(compare_limits.unique(codes))}))


def _load_records(codes):
    qs = Protein.objects.select_related('hole_hopping').filter(protein_code__in=codes)
    by_code = {p.protein_code: p for p in qs}
    return [by_code[c] for c in codes if c in by_code]


def _catalytic_columns(protein):
    """Full-alignment 0-indexed columns within this protein's catalytic domain — the
    columns whose (auth) residue number is <= the annotated catalytic-domain end_residue,
    matching what the 3D viewer keeps. Cached per protein (combination-independent)."""
    plot_cache = caches['plots']
    ck = 'catcols:{}:{}'.format(COMPARE_PLOTS_CACHE_VERSION, protein.protein_code)
    cached = plot_cache.get(ck)
    if cached is not None:
        return cached

    aligned = protein.aligned_sequence or ''
    cat_end = _domain_boundaries.get(protein.selected_structure or '', {}).get('end_residue')
    res_nums = parse_pdb_residue_numbers(protein.pdb_structure)
    seq_len = len(protein.sequence or '')
    if len(res_nums) != seq_len:
        res_nums = list(range(1, seq_len + 1))

    cols, seq_pos = [], 0
    for col, ch in enumerate(aligned):
        if ch == '-':
            continue
        auth = res_nums[seq_pos] if seq_pos < len(res_nums) else seq_pos + 1
        seq_pos += 1
        if cat_end is None or auth <= cat_end:
            cols.append(col)
    plot_cache.set(ck, cols, timeout=None)
    return cols


def _catalytic_keep_indices(records):
    """Union (sorted) of the selected proteins' catalytic-domain columns."""
    cols = set()
    for r in records:
        cols.update(_catalytic_columns(r))
    return sorted(cols)


def _catalytic_masked_seq(protein):
    """This protein's aligned sequence with every residue *beyond* its own catalytic-domain
    end_residue replaced by a gap ('-'). Used for the catalytic-domain view so each protein
    is shown only up to its own end: a kept alignment column that is catalytic for a longer
    protein no longer borrows a real (non-catalytic) residue from a shorter one — it shows a
    gap instead. Auth numbering is monotonic, so only trailing residues are masked, which
    every per-residue build in makePlots (counter-advances-on-non-gap) skips cleanly. Cached
    per protein (combination-independent)."""
    plot_cache = caches['plots']
    ck = 'catmask:{}:{}'.format(COMPARE_PLOTS_CACHE_VERSION, protein.protein_code)
    cached = plot_cache.get(ck)
    if cached is not None:
        return cached

    masked = catalytic_masked_sequence(protein.aligned_sequence or '', protein.sequence,
                                       protein.pdb_structure, protein.selected_structure)
    plot_cache.set(ck, masked, timeout=None)
    return masked


def _protein_masses(protein):
    """(catalytic-domain mass, total mass) in kDa, rounded to 0.1, from the protein's structure
    sequence: the polypeptide only (no glycans, metal or signal peptide; for crystal structures
    only the resolved construct). The catalytic value is None when the domain boundary is
    unknown, like the catalytic length."""
    from .protein_mass import polypeptide_kda
    total = polypeptide_kda(protein.sequence)
    catalytic = None
    if _catalytic_domain_length(protein) is not None:
        catalytic = polypeptide_kda(_catalytic_masked_seq(protein))
    return (round(catalytic, 1) if catalytic is not None else None,
            round(total, 1) if total is not None else None)


def _set_overview_values(protein):
    """The per-protein values overview.catalogue_overview reads, set on `protein`: net charge of
    all and of the surface-exposed residues and surface hydrophobicity, whole structure and
    catalytic domain (PQR atom charges up to end_residue, insertion codes included); RMSF, charge
    and hydrophobicity of each loop region; catalytic-domain length (the full length is
    protein.sequence_length); polypeptide masses (kDa) of the catalytic domain and the whole
    structure sequence."""
    end_res = _domain_boundaries.get(protein.selected_structure or '', {}).get('end_residue')
    protein.residue_profile = protein_residue_profile(protein, end_res)
    protein.cat_domain_charge = protein.residue_profile['charge']['cd']
    protein.loop_profile = protein_loop_profile(protein)
    protein.residue_sets = protein_residue_sets(protein)
    protein.cat_domain_length = _catalytic_domain_length(protein)
    protein.cat_domain_mass, protein.sequence_mass = _protein_masses(protein)
    domain = ''.join(c for c in _catalytic_masked_seq(protein) if c.isalpha()) if end_res is not None else ''
    protein.structural_features = protein_structural_features(protein, end_res, domain)


def _catalogue_fingerprint():
    """A hash of what the catalogue overview is made of (the proteins' fields, the lengths of
    their structure, PQR, SASA and RMSF texts, their substrates and regioselectivity, the pipeline
    files it reads, the app's code), so its cached copy is rebuilt after any data or code change
    without a version bump. One query per table, ~60 ms."""
    import hashlib
    from django.db.models.functions import Length
    digest = hashlib.md5()

    def add(rows):
        for row in rows:
            digest.update(repr(row).encode())

    add(Protein.objects.annotate(
        pqr_length=Length('electrostatics__pqr_file'),
        sasa_length=Length('surface_accesability__sasa_file'),
        pdb_length=Length('pdb_structure'),
        hh_length=Length('hole_hopping__paths_json'),
    ).order_by('pk').values_list(
        'pk', 'protein_code', 'protein_name', 'selected_structure', 'clade__name', 'organism__name',
        'organism__kingdom', 'organism__lineage', 'associated_CBMs', 'sequence', 'aligned_sequence',
        'second_sphere_motif', 'brace_motifs', 'loop_regions', 'activities', 'sequence_length',
        'disulfides', 'beta_core', 'binding_surface', 'hh_length',
        'pdb_length', 'pqr_length', 'sasa_length'))
    add(Protein.substrate.through.objects.order_by('pk').values_list('protein_id', 'substrate__name'))
    add(Protein.regioselectivity.through.objects.order_by('pk')
        .values_list('protein_id', 'regioselecitivity__name'))
    add(Protein.flexibility.through.objects.order_by('pk').annotate(rmsf_length=Length('flexibility__rmsf'))
        .values_list('protein_id', 'flexibility__method', 'rmsf_length'))
    from .dendrogram import DISTANCE_MATRIX_CSV
    from .overview import HYDROPHOBICITY_CSV
    # the files it reads, and the app's code (a change of the overview's records rebuilds it)
    files = [DISTANCE_MATRIX_CSV, _domain_boundaries_path, HYDROPHOBICITY_CSV,
             _Path(__file__).parents[2] / 'data' / 'django_outputs' / 'domain_ranges_auth.json']
    for path in files + sorted(_Path(__file__).parent.glob('*.py')):
        try:
            digest.update(f'{path}:{path.stat().st_mtime_ns}'.encode())
        except OSError:
            digest.update(str(path).encode())
    return digest.hexdigest()


def _catalogue_overview_bundle():
    """The overview of the whole catalogue (overview.catalogue_overview, as the browse page's
    panel) and the browse page's tree (dendrogram.average_linkage_tree), for the compare page's
    selection overview. Building it takes ~5 s (the per-protein PQR and loop profiles), so it is
    cached until the data change (_catalogue_fingerprint)."""
    from .dendrogram import average_linkage_tree
    plot_cache = caches['plots']
    ck = 'catalogue_overview:{}:{}'.format(COMPARE_PLOTS_CACHE_VERSION, _catalogue_fingerprint())
    bundle = plot_cache.get(ck)
    if bundle is None:
        proteins = list(Protein.objects.select_related(
            'organism', 'clade', 'surface_accesability', 'electrostatics', 'hole_hopping'
        ).prefetch_related('flexibility', 'substrate', 'regioselectivity').order_by('protein_name'))
        for protein in proteins:
            _set_overview_values(protein)
        matrix = DendrogramHandler().distance_matrix()
        bundle = {
            'overview': catalogue_overview(proteins, _catalytic_masked_seq),
            'tree': average_linkage_tree(matrix.values.tolist(), list(matrix.index)),
        }
        plot_cache.set(ck, bundle, timeout=None)
    return bundle


def compare_overview(request):
    """The compare page's selection overview. GET ids (comma-separated protein codes). Returns
    the catalogue overview of every protein (the charts draw the selection in colour, the rest of
    the catalogue in grey), the codes of the selection found in it, and the browse page's tree
    pruned to them."""
    from .dendrogram import prune_tree
    codes = [c for c in request.GET.get('ids', '').split(',') if c]
    if not codes:
        return JsonResponse({'error': 'missing ids'}, status=400)
    bundle = _catalogue_overview_bundle()
    known = {p['code'] for p in bundle['overview']['proteins']}
    selected = [c for c in dict.fromkeys(codes) if c in known]
    return JsonResponse({'overview': bundle['overview'], 'selected': selected,
                         'tree': prune_tree(bundle['tree'], set(selected))})


def _site_cells(protein):
    """The compare table's active-site columns, each {'text', 'title'} (title = the hover): the
    second-sphere motif with its spacing (its residues), the shielding residue over the copper's
    solvent side (its atom, distance and angle; 'None (open)' without one), and the motifs after
    the N-terminal and the internal His of the brace (the His and the residues after it, and the
    motif residues). Empty where the pipeline found nothing."""
    def residue(r):
        return f"{(r.get('name') or '').capitalize()}{r.get('residue', '')}"

    site = second_sphere(protein)
    cells = {'second_sphere': {'text': site['motif'],
                               'title': ' · '.join(residue(r) for r in site['residues'])}}
    shield = site['distal_residue']
    if shield:
        cells['shielding'] = {'text': residue(shield),
                              'title': f"{shield.get('atom', '')} {shield.get('distance', '')} Å from the copper, "
                                       f"{shield.get('angle', '')}° from the axis opposite the aromatic"}
    elif site['distal'] == OPEN_POSITION:
        cells['shielding'] = {'text': OPEN_POSITION, 'title': 'No residue covers the copper on the solvent side'}
    else:
        cells['shielding'] = {'text': '', 'title': ''}
    for which, motif in brace_motifs(protein).items():
        if motif['motif'] == NOT_DETERMINED:
            cells[which] = {'text': '', 'title': ''}
            continue
        title = motif['residues']
        if motif['motif_residues']:
            title += ' · motif residues ' + ', '.join(residue(r) for r in motif['motif_residues'])
        cells[which] = {'text': motif['motif'], 'title': title}
    return cells


def _binding_aromatics_cell(protein):
    """The compare table's aromatic binding residues: {'text': 'Tyr54, Trp178', 'title': each with
    its loop}; '–' when the surface has none, '' without a binding surface."""
    from .overview import protein_binding_aromatics
    found = protein_binding_aromatics(protein)
    if found is None:
        return {'text': '', 'title': 'No binding surface'}
    if not found:
        return {'text': '–', 'title': 'No Phe, Trp, Tyr or His (besides the brace) on the binding surface'}
    return {'text': ', '.join(a['label'] for a in found),
            'title': ' · '.join(f"{a['label']} ({a['region']})" for a in found)}


def _kda_text(value):
    """Mass for a table cell: '24.6', or '' when unknown."""
    return '' if value is None else f'{value:.1f}'


def _get_section_plots(codes, ptype, scale=None, full=False):
    """Per-protein plot JSON (in `codes` order) for one plot-type section.

    By default the plots are trimmed to the catalytic domain (matching the 3D viewer);
    `full=True` shows the whole length. Marker types (Sequence, Charge, Hole Hopping,
    Hydrophobicity) are assembled by column-subsetting the pre-generated full-alignment
    figures. The 3 line-gradient types (and any protein missing a pre-generated file) are
    regenerated live over the same column set (via makePlots keep_columns), so all plots
    share the same x-axis. Cached per (combination, type, scale, view)."""
    plot_cache = _plots_cache(codes)
    # Slugify type + variant so the cache key has no spaces (avoids memcached-key warnings).
    # `scale` is the variant selector: hydrophobicity scale OR flexibility method — include it
    # in the key whenever present so different variants cache separately.
    slug = ptype.replace(' ', '_')
    if scale:
        slug += '_' + str(scale).replace(' ', '_')
    view = 'full' if full else 'cat'
    suffix = 'sec_' + slug + '_' + view
    ck = _section_cache_key(codes, suffix)
    cached = plot_cache.get(ck)
    if cached is not None:
        return cached

    records = _load_records(codes)
    if full:
        keep = compute_keep_indices([list(r.aligned_sequence or '') for r in records])
    else:
        keep = _catalytic_keep_indices(records)

    plots = None
    # Pre-generated figures can only be column-subset for the full-length view (which
    # removes no residues, so plot-relative encodings are preserved). The catalytic-domain
    # view removes residues, so it is regenerated live (over fewer columns, hence cheap).
    if full and ptype in _PREGEN_MARKER_TYPES:
        plots = assemble_from_pregen([r.protein_code for r in records], ptype, scale, keep)
        if plots is not None and ptype == 'Hole Hopping':
            # Pre-generated one protein at a time: colour each row by its position here.
            plots = [json.dumps(recolor_hole_hopping_figure(json.loads(p), i))
                     for i, p in enumerate(plots)]
    if plots is None:
        # Catalytic view masks each protein beyond its own end (see _catalytic_masked_seq);
        # full-length view shows the real aligned sequence.
        seqs = ([r.aligned_sequence for r in records] if full
                else [_catalytic_masked_seq(r) for r in records])
        msa = {'name': [r.protein_code for r in records], 'seq': seqs}
        pm = makePlots(records, msa, keep_columns=(None if full else keep))
        plots = [pio.to_json(pm.make_plot_for_protein(ptype, i, scale_name=scale),
                             False, False, True, 'orjson') for i in range(len(records))]
    plot_cache.set(ck, plots, timeout=None)
    return plots


def _get_molstar_coloring(codes):
    """3D-viewer coloring + hydrophobicity scale names + export data for a selection,
    cached per combination. This needs the ~makeSequencePlots setup (makeQueryParams
    reads sequenceDict), so it is the one non-trivial build on page load."""
    plot_cache = _plots_cache(codes)
    # 'molstar2': the export holds a catalytic and a full-length version (a new key for this
    # bundle only, so the cached plots stay valid)
    ck = _section_cache_key(codes, 'molstar2')
    bundle = plot_cache.get(ck)
    if bundle is None:
        records = _load_records(codes)
        msa = {'name': [r.protein_code for r in records],
               'seq': [r.aligned_sequence for r in records]}
        pm = makePlots(records, msa)
        pm.makeSequencePlots(make_subplots(rows=len(records), cols=1), row=1, col=1)
        molstar = pm.makeQueryParams()
        molstar['catalyticDomainEnd'] = [
            _domain_boundaries.get(r.selected_structure or '', {}).get('end_residue')
            for r in records]
        bundle = {'molstar': molstar, 'scale_names': list(pm.hydrophobicityScaleNames),
                  # Flexibility methods available for this selection, as (value, label) for the
                  # dropdown; per-method 3D colouring is inside molstar['flexibilityMethods'].
                  'flex_methods': [(m, FLEX_METHOD_META.get(m, {}).get('label', m))
                                   for m in pm.flexibilityMethods],
                  'export': pm.makeExportData()}
        plot_cache.set(ck, bundle, timeout=None)
    return bundle


# 3D colouring fetched on demand: {kind in the request: key in the molstar bundle}
_LAZY_COLORINGS = {'hydrophobicityScales': 'hydrophobicityScales',
                   'flexibilityMethods': 'flexibilityMethods',
                   'flexibilityMethodsFull': 'flexibilityMethodsFull'}


def compare_coloring(request):
    """One 3D colouring the compare page does not carry: a hydrophobicity scale or a flexibility
    method. GET ids, kind (hydrophobicityScales | flexibilityMethods | flexibilityMethodsFull)
    and name (the scale or method). Returns {'coloring': [per protein]} from the cached bundle."""
    codes = [c for c in request.GET.get('ids', '').split(',') if c]
    kind = request.GET.get('kind', '')
    name = request.GET.get('name', '')
    if not codes or kind not in _LAZY_COLORINGS:
        return JsonResponse({'error': 'missing ids or kind'}, status=400)
    if compare_limits.too_many(codes):
        return JsonResponse({'error': compare_limits.limit_message(len(set(codes)))}, status=400)
    colorings = _get_molstar_coloring(codes)['molstar'].get(_LAZY_COLORINGS[kind]) or {}
    if name not in colorings:
        return JsonResponse({'error': f'unknown {kind} name'}, status=404)
    return JsonResponse({'coloring': colorings[name]})


def compare_plots(request):
    """Lazy per-section plot data for the compare page.

    GET params: ids (comma-separated protein codes, URL order), ptype, optional scale
    (hydrophobicity), optional full=1 for the whole length (default: catalytic domain only).
    Returns {'plots_html': [<plotly json str>, ...]} in protein order.
    """
    ids = request.GET.get('ids', '')
    ptype = request.GET.get('ptype', '')
    scale = request.GET.get('scale') or None
    full = request.GET.get('full') == '1'
    codes = [c for c in ids.split(',') if c]
    if not codes or not ptype:
        return JsonResponse({'error': 'missing ids or ptype'}, status=400)
    if compare_limits.too_many(codes):
        return JsonResponse({'error': compare_limits.limit_message(len(set(codes)))}, status=400)
    try:
        return JsonResponse({'plots_html': _get_section_plots(codes, ptype, scale, full=full)})
    except Exception as e:
        return JsonResponse({'error': str(e)}, status=500)


def compare_proteins(request, protein_ids):
    protein_id_list = compare_limits.unique(protein_ids.split(','))
    # at most MAX_COMPARE_PROTEINS, or a whole family -- shown in its fixed order, so that the
    # comparison built ahead for it is the one served
    if compare_limits.too_many(protein_id_list):
        return _back_to_selection(protein_id_list)
    named_groups = protein_groups.parse(request.GET.get('groups'), protein_id_list)
    ordered = (protein_groups.ordered_codes(named_groups, protein_id_list) if named_groups
               else compare_limits.canonical_order(protein_id_list))
    if ordered != protein_id_list:
        url = reverse('compare_proteins', kwargs={'protein_ids': ','.join(ordered)})
        query = request.GET.urlencode()
        return redirect(url + ('?' + query if query else ''))
    # Preserve order from URL parameter - create a mapping and sort
    protein_records_queryset = Protein.objects.select_related('hole_hopping').filter(
        protein_code__in=protein_id_list)
    # Create a dict for O(1) lookup
    protein_dict = {p.protein_code: p for p in protein_records_queryset}
    # Reorder to match protein_id_list order
    protein_records = [protein_dict[pid] for pid in protein_id_list if pid in protein_dict]

    # protein_codes_list = [p.protein_code for p in protein_records]

    # Calculate sequence length for each protein (amino acids, excluding gaps)
    for protein in protein_records:
        try:
            seq_to_use = None
            if hasattr(protein, 'aligned_sequence') and protein.aligned_sequence:
                seq_to_use = protein.aligned_sequence
            elif protein.sequence:
                seq_to_use = protein.sequence

            if seq_to_use:
                # Remove gaps (-), whitespace, and newlines to get actual amino acid sequence length
                sequence_clean = seq_to_use.replace('-', '').replace(' ', '').replace('\n', '').replace('\r',
                                                                                                        '').replace(
                    '\t', '')
                # Count only valid amino acid letters (A-Z, excluding gaps and special chars)
                sequence_clean = ''.join([c for c in sequence_clean if c.isalpha()])
                protein.sequence_length = len(sequence_clean)
            else:
                protein.sequence_length = 0
        except Exception:
            protein.sequence_length = 0

    # Build MSA dict from protein.aligned_sequence stored in DB
    msa_names = []
    msa_seqs = []
    for p in protein_records:
        if p.aligned_sequence:
            msa_names.append(p.protein_code)
            msa_seqs.append(p.aligned_sequence)

    msa_list = {"name": msa_names, "seq": msa_seqs}
    codes = [p.protein_code for p in protein_records]

    # Groups: every protein gets its group's index, name and colour (bands in the table and
    # the plots, 3D colouring); with groups, the proteins outside them are "Ungrouped"
    all_groups = protein_groups.with_ungrouped(named_groups, codes)
    group_of = {c: (i, g) for i, g in enumerate(all_groups) for c in g['codes']}
    for protein in protein_records:
        index, group = group_of.get(protein.protein_code, (None, None))
        protein.group_index = index
        protein.group_name = group['name'] if group else ''
        protein.group_color = group['color'] if group else ''

    # On page load, build only the 3D-viewer coloring (the one piece that needs the
    # makeSequencePlots setup) and the Sequence section. The Sequence section is served
    # from pre-generated figures (instant column-subset); every other plot type — and
    # non-default hydrophobicity scales — is fetched lazily via the compare_plots
    # endpoint when its dropdown opens. This avoids building the full plot set up front.
    _molstar = _get_molstar_coloring(codes)
    molstarColoring = _molstar['molstar']
    hydrophobicity_scale_names = _molstar['scale_names']
    flex_methods = _molstar.get('flex_methods', [])
    export_data = _molstar['export']

    # A large comparison opens with every plot section closed, the Sequence one too: it is then
    # fetched on opening like the others, instead of being part of the page.
    plots_closed = compare_limits.plots_closed(codes)
    sequence_plots = [] if plots_closed else _get_section_plots(codes, 'Sequence')
    plots_data = [
        {'type': t, 'label': _COMPARE_PLOT_LABELS.get(t, t),
         'description': _COMPARE_PLOT_DESCRIPTIONS.get(t, 'No description available.'),
         'plots_html': sequence_plots if t == 'Sequence' else [],
         # 1/0, not True/False: the template also writes this list into the page's script
         'open': int(t == 'Sequence' and not plots_closed)}
        for t in _COMPARE_PLOT_TYPES
    ]
    # The page carries the 3D colouring of the default hydrophobicity scale and flexibility
    # method only; the other scales and methods (two thirds of the colouring data) are fetched
    # from compare_coloring when chosen.
    molstarColoring = {**molstarColoring, **{key: {} for key in _LAZY_COLORINGS}}

    def clean_value(value):
        """Clean value to handle NaN and None properly"""
        if value is None:
            return ''
        value_str = str(value).strip()
        if value_str.lower() == 'nan' or value_str == '':
            return ''
        return value_str
    
    table_rows = []
    for idx, protein in enumerate(protein_records):
        # Get uniprot value, handling cases where field might not exist yet
        uniprot_value = ''
        uniprot_raw = None
        try:
            # Directly access the uniprot field
            uniprot_raw = protein.uniprot
            if uniprot_raw is not None and str(uniprot_raw).strip():
                uniprot_value = clean_value(str(uniprot_raw))
        except AttributeError:
            # Field doesn't exist in database yet (migration not run)
            uniprot_value = ''
        except Exception:
            # Any other error, just use empty string
            uniprot_value = ''
        
        _code = protein.protein_code or ''
        _alias_parts = [p.strip() for p in (protein.other_names or '').split(';') if p.strip() and p.strip() != _code]
        _display_name = '; '.join([_code] + _alias_parts) if _alias_parts else _code

        # Substrates, oxidised positions and EC numbers as the papers report them (step 2j),
        # numbered by the row's References column
        record = activity_data(protein)
        references_html, ref_index = _format_references(protein, activity_references(record))
        substrate_names = [clean_value(s.name) for s in protein.substrate.all() if clean_value(s.name)]
        activity_cells = table_cells(record, ref_index.get, substrate_names)

        table_rows.append({
            'index': idx,
            'group_name': protein.group_name,
            'group_color': protein.group_color,
            'name': _display_name,
            'clade': clean_value(protein.clade.name if protein.clade else ''),
            'catalytic_length': _catalytic_domain_length(protein) or '',
            'sequence_length': protein.sequence_length if protein.sequence_length is not None else 0,
            'catalytic_mass': _kda_text(_protein_masses(protein)[0]),
            'sequence_mass': _kda_text(_protein_masses(protein)[1]),
            'structure_display': {'AlphaFold DB': 'AlphaFold', 'own model': 'AlphaFold 3',
                                  'CAZyme3D': 'CAZyme3D', 'PDB': 'PDB'}[
                                      structure_source(protein.selected_structure)],
            'selected_structure': clean_value(protein.selected_structure or ''),
            'pdb_codes': clean_value(protein.pdb_code),
            'uniprot': uniprot_value,
            'uniprot_notes': _uniprot_notes(protein),
            'uniparc': clean_value(getattr(protein, 'uniparc', None) or ''),
            'organism': clean_value(protein.organism.name if protein.organism else ''),
            # NCBI lineage above the species; CAZy's kingdom without one
            'taxonomy': clean_value((lineage_text(_lineage(protein.organism)) or protein.organism.kingdom)
                                    if protein.organism else ''),
            'ec_numbers': clean_value(protein.ec_numbers),
            'substrates': clean_value(', '.join(substrate_names)),
            'regioselectivity': regioselectivity_label(r.name for r in protein.regioselectivity.all()),
            # second-sphere motif, shielding residue, brace histidine motifs (pipeline 8c-8e)
            **_site_cells(protein),
            'binding_aromatics': _binding_aromatics_cell(protein),
            'associated_CBMs': clean_value(protein.associated_CBMs),
            'references': references_html,
            # step 2j: what the papers report, citing the row's reference numbers
            'ec_cells': activity_cells['ec'],
            'ec_differs': activity_cells['ec_differs'],
            'reported_substrates': activity_cells['substrates'],
            'regio_refs': activity_cells['regio'],
        })

    # the count tiles above the table
    summary = catalogue_summary(protein_records)
    summary['catalogue'] = Protein.objects.count()
    summary['predicted_sources'] = ' · '.join(
        f'{label} {summary[key]}' for key, label in
        (('alphafold', 'AlphaFold DB'), ('own_models', 'AlphaFold 3'), ('cazyme3d', 'CAZyme3D')) if summary[key])

    from .similarity import available_loops
    context: dict = {'protein_records': protein_records,
                     'summary': summary,
                     'activity_export': [EXPORT_HEADER] + export_rows(protein_records),
                     'protCount': len(protein_records),
                     "plots_data": plots_data,
                     'molstarColoring': molstarColoring,
                     # the binding-surface plane of each protein (step 7c), for the viewer's
                     # Top / Side / Active site views; not part of the cached colouring
                     'binding_planes': [_viewer_plane(p) for p in protein_records],
                     # residue numbers of the aromatic binding residues, for "Show as sticks"
                     'binding_aromatic_numbers': [[a['residue'] for a in (protein_binding_aromatics(p) or [])]
                                                  for p in protein_records],
                     # the groups ({name, color, codes, ungrouped?}) and each protein's group
                     # index, for the page's script; groups_param for the links that keep them
                     'protein_groups': {'groups': all_groups,
                                        'of': [p.group_index for p in protein_records]},
                     'groups_param': protein_groups.to_param(named_groups) if named_groups else '',
                     'plots_closed': plots_closed,
                     'plots_closed_above': compare_limits.plots_closed_above(),
                     'table_rows':table_rows,
                     'protein_codes_list': protein_id_list,
                     'hydrophobicity_scale_names': hydrophobicity_scale_names,
                     'flex_methods': flex_methods,
                     'analysis_loops': available_loops(protein_records),
                     # the loops present, in sequence order, for the viewer's "Show as sticks"
                     'stick_loops': [name for name in ('L2', 'L3', 'LS', 'L8', 'LC')
                                     if name in available_loops(protein_records)],
                     'export_data': export_data}

    return render(request, 'catalogue/compare_protein.html', context)


def compare_analysis(request):
    """Descriptive similarity / conservation for the selected set (compare page's
    sequence logo and similarity heatmap). GET: ids, region, property, optional flex. Returns the
    JSON from catalogue.similarity.compute(), cached per (combination, region, property)."""
    from .similarity import compute
    ids = request.GET.get('ids', '')
    region = request.GET.get('region', 'loops')
    prop = request.GET.get('property', 'charge')
    flex = request.GET.get('flex', 'cabsflex')
    exposed = request.GET.get('exposed', '') in ('1', 'true', 'True')
    # min_occ=0 (the sequence logo) keeps every column of the region, so it spans the same
    # columns as the compare page's plots; the panel's own requests keep the default.
    try:
        min_occ = min(1.0, max(0.0, float(request.GET.get('min_occ', 0.3))))
    except ValueError:
        min_occ = 0.3
    codes = [c for c in ids.split(',') if c]
    if len(codes) < 2:
        return JsonResponse({'ok': False, 'error': 'Select at least two proteins.'}, status=400)
    if compare_limits.too_many(codes):
        return JsonResponse({'ok': False, 'error': compare_limits.limit_message(len(set(codes)))}, status=400)

    ck = _section_cache_key(codes, f'analysis:{region}:{prop}:{flex}:{int(exposed)}:{min_occ}:v15')
    plot_cache = _plots_cache(codes)
    cached = plot_cache.get(ck)
    if cached is not None:
        return JsonResponse(cached)

    records = _load_records(codes)
    records = [r for r in records if r.aligned_sequence]
    if len(records) < 2:
        return JsonResponse({'ok': False, 'error': 'Need at least two aligned proteins.'}, status=400)
    msa = {'name': [r.protein_code for r in records],
           'seq': [r.aligned_sequence for r in records]}
    result = compute(records, msa, region, prop, flex, min_occ=min_occ, exposed=exposed)
    plot_cache.set(ck, result, timeout=None)
    return JsonResponse(result)


def sequence_search(request):
    """Browse-page search bar in Sequence mode. POST JSON {query}: a pasted protein sequence
    (plain or FASTA). Returns catalogue hits ranked by local alignment, or exact motif
    matches for short queries (see catalogue.sequence_search)."""
    from .sequence_search import QueryError, parse_query, search, SIGNIFICANT_EVALUE
    if request.method != 'POST':
        return JsonResponse({'ok': False, 'error': 'POST required.'}, status=405)
    try:
        query = json.loads(request.body or b'{}').get('query', '')
    except (ValueError, AttributeError):
        return JsonResponse({'ok': False, 'error': 'Malformed request.'}, status=400)
    try:
        header, seq = parse_query(str(query))
        targets = list(Protein.objects.exclude(sequence='').values_list('protein_code', 'sequence'))
        result = search(seq, targets)
    except QueryError as e:
        return JsonResponse({'ok': False, 'error': str(e)}, status=400)

    info = {code: (clade, other_names) for code, clade, other_names in
            Protein.objects.values_list('protein_code', 'clade__name', 'other_names')}
    for hit in result['hits']:
        clade, other_names = info.get(hit['code'], ('', ''))
        hit['clade'] = clade or ''
        hit['other_names'] = other_names or ''
    return JsonResponse({'ok': True, 'header': header, 'query_length': len(seq),
                         'motif': seq if result['mode'] == 'motif' else '',
                         'significant_evalue': SIGNIFICANT_EVALUE, **result})


def _pqr_to_pdb_polar_only(pqr_text, polar_charge_threshold=0.18):
    """Convert PQR text to PDB format, keeping only polar hydrogens (charge > threshold).
    Polar H: N-H, O-H, S-H — charge in AMBER ff consistently > 0.18.
    Non-polar H: C-H — charge typically 0.02–0.16.
    The PQR charge column maps to PDB occupancy column position; Mol* tolerates this.
    """
    lines = []
    for line in pqr_text.splitlines():
        if not line.startswith(('ATOM', 'HETATM')):
            lines.append(line)
            continue
        parts = line.split()
        if len(parts) < 10:
            lines.append(line)
            continue
        atom_name = parts[2].strip()
        # Identify hydrogen atoms: name starts with H
        if atom_name.startswith('H'):
            try:
                charge = float(parts[-2])  # second-to-last column in PQR
            except (ValueError, IndexError):
                lines.append(line)
                continue
            if abs(charge) < polar_charge_threshold:
                continue  # skip non-polar H
        lines.append(line)
    return '\n'.join(lines)


def _with_disulfide_conect(pdb_text, disulfides):
    """pdb_text with CONECT records for its disulfide bonds (Protein.disulfides, residue
    numbers), merged with its own CONECT records into one block sorted by atom serial, each
    pair listed both ways -- the form in which Mol* (pdbe-molstar 3.5) takes bonds from CONECT.
    By distance alone Mol* misses the longer S-S bonds (7ZJB Cys64-Cys175 at 2.72 Å,
    CsAA10 Cys91-Cys205 at 2.32 Å), and a click on a Cys dot should show the bond."""
    try:
        bonds = json.loads(disulfides) if disulfides else []
    except (TypeError, ValueError):
        bonds = []
    if not pdb_text or not bonds:
        return pdb_text
    lines = pdb_text.split('\n')
    sulfur, partners = {}, defaultdict(set)
    for line in lines:
        try:
            if line.startswith('ATOM') and line[12:16].strip() == 'SG':
                sulfur.setdefault(int(line[22:26]), int(line[6:11]))
            elif line.startswith('CONECT'):
                serials = [int(line[c:c + 5]) for c in range(6, len(line.rstrip()), 5)
                           if line[c:c + 5].strip()]
                partners[serials[0]].update(serials[1:])
        except (ValueError, IndexError):
            continue
    pairs = [(sulfur.get(bond[0]), sulfur.get(bond[1])) for bond in bonds]
    pairs = [(a, b) for a, b in pairs if a and b]
    if not pairs:
        return pdb_text
    for a, b in pairs:
        partners[a].add(b)
        partners[b].add(a)
    conect = []
    for a in sorted(partners):
        bonded = sorted(partners[a])
        conect += [f'CONECT{a:5d}' + ''.join(f'{b:5d}' for b in bonded[i:i + 4])
                   for i in range(0, len(bonded), 4)]
    body = [line for line in lines if not line.startswith('CONECT')]
    end = next((i for i, line in enumerate(body)
                if line.startswith('END') and not line.startswith('ENDMDL')), len(body))
    return '\n'.join(body[:end] + conect + body[end:])


_HIS_RING_ATOMS = {'CG', 'ND1', 'CD2', 'CE1', 'NE2'}


def _brace_axis(protein):
    """The vector from the N-His ring centre to the Internal-His ring centre in the protein's
    structure (the aligned coordinates the viewer shows), or None."""
    try:
        his1, his2 = (int(n) for n in protein.histidine_brace.strip('()').split(','))
    except (AttributeError, ValueError):
        return None
    atoms = {his1: [], his2: []}
    for line in (protein.pdb_structure or '').splitlines():
        if not line.startswith('ATOM') or line[12:16].strip() not in _HIS_RING_ATOMS:
            continue
        try:
            number = int(line[22:26])
        except ValueError:
            continue
        if number in atoms and line[26] == ' ':
            atoms[number].append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if not atoms[his1] or not atoms[his2]:
        return None
    centre = {n: [sum(c[i] for c in xyz) / len(xyz) for i in range(3)] for n, xyz in atoms.items()}
    return [round(centre[his2][i] - centre[his1][i], 3) for i in range(3)]


def _viewer_plane(protein):
    """The binding-surface plane of a protein for the 3D viewer's views ({normal, copper, depth}
    of step 7c), with `brace` (N-His -> Internal-His), which Top and Side put horizontal."""
    if not protein.binding_plane:
        return None
    plane = json.loads(protein.binding_plane)
    plane['brace'] = _brace_axis(protein)
    return plane


def _annotate_pdb_for_molstar(pdb_text, protein_code, structure_code):
    """Rewrite TITLE and COMPND MOLECULE in a PDB string so Mol* shows
    protein_code as the entity name and structure_code in the title."""
    if not pdb_text:
        return pdb_text

    lines = pdb_text.split('\n')
    result = []
    compnd_mol_replaced = False
    title_replaced = False

    for line in lines:
        record = line[:6].rstrip()

        if record == 'TITLE':
            if not title_replaced:
                label = f'{protein_code} ({structure_code})'
                result.append(f'TITLE     {label[:70]}')
                title_replaced = True
            # drop continuation TITLE lines — we wrote a single replacement
            continue

        if record == 'COMPND' and not compnd_mol_replaced:
            if 'MOLECULE:' in line.upper():
                mol_pos = line.upper().index('MOLECULE:')
                result.append(f'{line[:mol_pos]}MOLECULE: {protein_code};')
                compnd_mol_replaced = True
                continue

        result.append(line)

    # If neither record existed, insert minimal versions before the first ATOM
    if not title_replaced:
        atom_idx = next((i for i, l in enumerate(result) if l.startswith('ATOM')), 0)
        result.insert(atom_idx, f'TITLE     {protein_code} ({structure_code})')
    if not compnd_mol_replaced:
        atom_idx = next((i for i, l in enumerate(result) if l.startswith('ATOM')), 0)
        result.insert(atom_idx, f'COMPND   2 MOLECULE: {protein_code};')
        result.insert(atom_idx, f'COMPND    MOL_ID: 1;')

    return '\n'.join(result)


def get_pdb_structures(request):

    fetch_start = time.time()

    if request.method != 'GET':
        return JsonResponse({"error": None}, status=405)

    protein_ids = request.GET.get('ids', '')
    if not protein_ids:
        return JsonResponse({"error": None}, status=400)

    protonated = request.GET.get('protonated', '0') == '1'

    protein_id_list = [pid.strip() for pid in protein_ids.split(',') if pid.strip()]
    if not protein_id_list:
        return JsonResponse({"error": None}, status=400)
    if compare_limits.too_many(protein_id_list):
        return JsonResponse({"error": compare_limits.limit_message(len(set(protein_id_list)))}, status=400)

    try:
        db_start = time.time()
        structure_codes = {}  # protein_code -> selected_structure (pdb/AF code) for Mol* labels
        if protonated:
            # Build dict with PQR-as-PDB (polar H only), fall back to pdb_structure
            from django.db.models import Prefetch
            protein_objs = Protein.objects.filter(
                protein_code__in=protein_id_list
            ).prefetch_related('electrostatics')
            structures_dict = {}
            for prot in protein_objs:
                pqr = None
                if hasattr(prot, 'electrostatics') and prot.electrostatics:
                    pqr = getattr(prot.electrostatics, 'pqr_file', None)
                raw = _pqr_to_pdb_polar_only(pqr) if pqr else prot.pdb_structure
                structures_dict[prot.protein_code] = _annotate_pdb_for_molstar(
                    _with_disulfide_conect(raw, prot.disulfides),
                    prot.protein_code, prot.selected_structure or ''
                )
                structure_codes[prot.protein_code] = prot.selected_structure or ''
        else:
            structures = Protein.objects.filter(
                protein_code__in=protein_id_list
            ).values('protein_code', 'pdb_structure', 'selected_structure', 'disulfides')
            structures_dict = {}
            for p in structures:
                structures_dict[p['protein_code']] = _annotate_pdb_for_molstar(
                    _with_disulfide_conect(p['pdb_structure'], p['disulfides']),
                    p['protein_code'], p['selected_structure'] or ''
                )
                structure_codes[p['protein_code']] = p['selected_structure'] or ''
        db_time = time.time() - db_start

    except Exception as e:
        return JsonResponse({"error": f"Database error: {str(e)}"}, status=500)

    missing = set(protein_id_list) - set(structures_dict.keys())
    fetch_time = time.time() - fetch_start
    response_data = {
        "structures": structures_dict,
        "structure_codes": structure_codes,
        "count": len(structures_dict),
        "fetch_time": round(fetch_time, 3),
    }
    if missing:
        response_data["error"] = f"Some proteins not found: {', '.join(missing)}"
        response = JsonResponse(response_data, status=206)
    else:
        response = JsonResponse(response_data)

    if not protonated:
        response['Cache-Control'] = 'public, max-age=3600'
        response['ETag'] = f'"pdb-{hash(tuple(sorted(protein_id_list)))}"'

    return response


def _strip_protein_prefix(code):
    """Remove leading numeric prefix (e.g. '016_') from protein codes for user-facing filenames."""
    import re
    return re.sub(r'^\d+_', '', code)


def download_alignment(request):
    """The structural alignment of the whole catalogue as FASTA: ?view=catalytic (default) with
    every residue past a protein's catalytic-domain end left out, or ?view=full. Columns no
    protein has a residue in are dropped. Proteins by family, then code."""
    full = request.GET.get('view') == 'full'
    view = 'full_length' if full else 'catalytic_domain'
    response = HttpResponse(alignment_fasta(full), content_type='text/plain')
    response['Content-Disposition'] = f'attachment; filename="lpmocat_alignment_{view}.fasta"'
    return response


def alignment_fasta(full=False):
    """download_alignment's FASTA text (the whole catalogue)."""
    proteins = [p for p in Protein.objects.select_related('clade', 'organism')
                if p.aligned_sequence]
    family = lambda p: int(re.sub(r'\D', '', p.clade.name) or 0)
    proteins.sort(key=lambda p: (family(p), p.protein_code))
    seqs = [p.aligned_sequence if full else _catalytic_masked_seq(p) for p in proteins]
    width = max((len(s) for s in seqs), default=0)
    seqs = [s.ljust(width, '-') for s in seqs]
    columns = [c for c in range(width) if any(s[c] != '-' for s in seqs)]
    lines = []
    for p, s in zip(proteins, seqs):
        lines.append(f'>{p.protein_code} {p.selected_structure} {p.clade.name} {p.organism.name}')
        lines.append(''.join(s[c] for c in columns))
    return '\n'.join(lines) + '\n'


def _pqr_to_pdb(pqr_text):
    """A PDB2PQR file as a plain PDB file with every hydrogen: the atom records up to the
    coordinates (fixed PDB columns in PDB2PQR's output), occupancy 1.00, B-factor 0.00 and the
    element, in place of the charge and radius."""
    lines = []
    for line in pqr_text.splitlines():
        if line.startswith(('ATOM', 'HETATM')):
            name = line[12:16].strip()
            element = next((c for c in name if c.isalpha()), '')
            lines.append(f'{line[:54]:<54}{1.0:6.2f}{0.0:6.2f}          {element:>2}')
        elif line.startswith(('TER', 'END')):
            lines.append(line)
    return '\n'.join(lines) + '\n'


def _trim_to_catalytic(text, cat_end):
    """A PDB or PQR text without the residues past the catalytic-domain end (auth numbering):
    their ATOM records, the SSBOND records naming them and the CONECT records of their atoms.
    HETATM records (the active-site metal) stay."""
    if not text or cat_end is None:
        return text
    kept, removed = [], set()
    for line in text.splitlines():
        try:
            if line.startswith('ATOM') and int(line[22:26]) > cat_end:
                removed.add(int(line[6:11]))
                continue
            if line.startswith('SSBOND') and max(int(line[17:21]), int(line[31:35])) > cat_end:
                continue
            if line.startswith('CONECT') and any(
                    int(line[c:c + 5]) in removed for c in range(6, len(line.rstrip()), 5)
                    if line[c:c + 5].strip()):
                continue
        except ValueError:
            pass
        kept.append(line)
    return '\n'.join(kept) + '\n'


def structure_files(proteins, catalytic=False, protonated=False):
    """[(file name, text)] of the proteins' structures as the viewer shows them: catalytic=True
    leaves out the residues past each protein's catalytic-domain end (names end in
    _catalytic_domain); protonated=True gives, per protein, a PDB file with every hydrogen and the
    PDB2PQR file itself (charges and radii)."""
    files = []
    for p in proteins:
        name = _strip_protein_prefix(p.protein_code)
        cat_end = (_domain_boundaries.get(p.selected_structure or '', {}).get('end_residue')
                   if catalytic else None)
        if cat_end is not None:
            name += '_catalytic_domain'
        if protonated:
            electrostatics = getattr(p, 'electrostatics', None)
            pqr = getattr(electrostatics, 'pqr_file', None) if electrostatics else None
            if pqr:
                pqr = _trim_to_catalytic(pqr, cat_end)
                files.append((f'{name}_protonated.pdb', _pqr_to_pdb(pqr)))
                files.append((f'{name}.pqr', pqr))
        elif p.pdb_structure:
            pdb = p.pdb_structure
            pdb = pdb.decode('utf-8') if isinstance(pdb, bytes) else pdb
            files.append((f'{name}.pdb', _trim_to_catalytic(pdb, cat_end)))
    return files


def download_pdb_structures(request):
    """The structures of the given proteins as PDB files (one file, or a zip for several), or
    with ?protonated=1 the protonated ones: for each protein a PDB file with every hydrogen and
    the PDB2PQR file itself (charges and radii), always as a zip. ?view=catalytic leaves out the
    residues past each protein's catalytic-domain end (file names end in _catalytic_domain)."""
    protein_ids = request.GET.get('ids', '')
    protein_id_list = [pid.strip() for pid in protein_ids.split(',') if pid.strip()]
    protonated = request.GET.get('protonated') == '1'
    catalytic = request.GET.get('view') == 'catalytic'
    proteins = Protein.objects.filter(protein_code__in=protein_id_list).select_related('electrostatics')
    files = structure_files(proteins, catalytic=catalytic, protonated=protonated)
    as_bytes = lambda t: t if isinstance(t, bytes) else str(t).encode('utf-8')

    if len(files) == 1:
        response = HttpResponse(as_bytes(files[0][1]), content_type='chemical/x-pdb')
        response['Content-Disposition'] = f'attachment; filename="{files[0][0]}"'
        return response

    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, 'w', zipfile.ZIP_DEFLATED) as zf:
        for filename, text in files:
            zf.writestr(filename, as_bytes(text))
    zip_filename = 'protonated_structures.zip' if protonated else 'pdb_structures.zip'
    response = HttpResponse(buffer.getvalue(), content_type='application/zip')
    response['Content-Disposition'] = f'attachment; filename="{zip_filename}"'
    return response


def download_comparison(request):
    """The zip of one comparison (?ids= in the page's order): README, the protein table, the
    substrate evidence, the per-residue data and the alignment in both views, and the
    structures (downloads.comparison_zip)."""
    codes = compare_limits.unique(request.GET.get('ids', '').split(','))
    if not codes:
        return HttpResponse('No proteins given', status=400)
    if compare_limits.too_many(codes):
        return HttpResponse('Too many proteins for one comparison', status=400)
    groups = protein_groups.with_ungrouped(
        protein_groups.parse(request.GET.get('groups'), codes), codes)
    response = HttpResponse(downloads.comparison_zip(codes, groups), content_type='application/zip')
    response['Content-Disposition'] = 'attachment; filename="lpmocat_comparison.zip"'
    return response


# The Downloads section (About page): the files build_downloads writes, in its order, with what each holds
DOWNLOAD_FILES = [
    ('lpmocat_catalogue.zip', 'Everything below in one archive, with a README',
     'all'),
    ('lpmocat_proteins.csv', 'One row per protein: family, organism and taxonomy, structure and '
     'its source, UniProt, EC numbers, substrates, regioselectivity, active-site metal, brace '
     'and second-sphere motifs, shielding residue, aromatic binding residues, CBMs, lengths, '
     'masses, disulfide bonds, '
     'references', 'table'),
    ('lpmocat_substrate_evidence.csv', 'Every activity claim with its source (CAZy, Swiss-Prot, '
     'BRENDA, literature) and papers, inactive results included', 'table'),
    ('lpmocat_residues_catalytic_domain.csv', 'One row per residue of the catalytic domains: '
     'secondary structure, core strand, loop, binding surface, disulfide partner, relative SASA, '
     'charge, hydrophobicity, flexibility', 'table'),
    ('lpmocat_residues_full_length.csv', 'The same for the whole structures', 'table'),
    ('lpmocat_alignment_catalytic_domain.fasta', 'Structural alignment (US-align) of the '
     'catalytic domains', 'sequence'),
    ('lpmocat_alignment_full_length.fasta', 'Structural alignment of the whole structures',
     'sequence'),
    ('lpmocat_distance_matrix.csv', 'BLOSUM62 distances between the catalytic domains (the '
     'similarity tree)', 'sequence'),
    ('lpmocat_hole_hopping_paths.json', 'Hole-hopping paths to the copper (EHPath) per protein',
     'table'),
    ('lpmocat_structures.zip', 'All structures as PDB files (aligned coordinates): full length '
     'and catalytic domain', 'structure'),
    ('lpmocat_structures_protonated_full_length.zip', 'All whole structures with hydrogens '
     '(PDB2PQR, AMBER): a PDB file and the PQR file (charges, radii) per protein', 'structure'),
]


def _download_files():
    """The Downloads section's files: name, what it holds, kind and size (None if not built)."""
    files = []
    for name, description, kind in DOWNLOAD_FILES:
        path = downloads.DOWNLOADS_DIR / name
        size = path.stat().st_size if path.exists() else None
        files.append({'name': name, 'description': description, 'kind': kind,
                      'size': _size_text(size) if size is not None else None})
    return files


def downloads_page(request):
    """The downloads are a section of the About page."""
    return redirect('/#downloads')


def _size_text(size):
    for unit in ('B', 'kB', 'MB', 'GB'):
        if size < 1000 or unit == 'GB':
            return f'{size:.0f} {unit}' if unit == 'B' else f'{size:.1f} {unit}'
        size /= 1000


def download_file(request, name):
    """One file of the Downloads page (only the names it lists)."""
    if name not in {n for n, _, _ in DOWNLOAD_FILES}:
        raise Http404('No such download')
    path = downloads.DOWNLOADS_DIR / name
    if not path.exists():
        raise Http404('Not built yet')
    return FileResponse(open(path, 'rb'), as_attachment=True, filename=name)
