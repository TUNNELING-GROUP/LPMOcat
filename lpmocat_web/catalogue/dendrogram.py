import pandas as pd
from .models import Protein
from collections import Counter
from html import escape
from pathlib import Path
from .protein_mass import polypeptide_kda

# Pairwise distances between the catalytic domains (BLOSUM62 over the aligned residues up to
# each catalytic domain end), labelled by struct_id: pipeline step 9
# (distance_matrix). The tree is built from them in the browser (circular_dendrogram.js).
DISTANCE_MATRIX_CSV = Path(__file__).resolve().parents[2] / 'data' / '09_distance_matrix' / 'distance_matrix.csv'

# LPMO family colours of the dendrogram and the browse page's overview charts: the eight
# categorical colours of the dataviz reference palette, assigned so that families which touch
# (neighbours around the tree, neighbouring segments of a stacked bar in FAMILY_STACK_ORDER)
# stay apart under protanopia/deuteranopia (OKLab dE >= 9.2) and normal vision (>= 16.3).
CLADE_COLORS = {
    'AA9': '#4a3aa7',   # violet
    'AA10': '#2a78d6',  # blue
    'AA11': '#eb6834',  # orange
    'AA13': '#e87ba4',  # magenta
    'AA14': '#e34948',  # red
    'AA15': '#eda100',  # yellow
    'AA16': '#1baf7a',  # aqua
    'AA17': '#008300',  # green
}

# Families in the order they follow each other around the dendrogram (2026-09-16 data): the
# overview charts stack in this order, so touching segments are tree neighbours, whose
# colours were chosen to differ. Recheck both if a new family or tree layout changes it.
FAMILY_STACK_ORDER = ['AA9', 'AA11', 'AA16', 'AA10', 'AA17', 'AA13', 'AA15', 'AA14']

def average_linkage_tree(matrix, labels):
    """The browse page's tree as nested dicts, leaves {'name'} and internal nodes {'distance',
    'children'}: a port of hierarchical_clustering.js (average linkage, then optimizeTree) that
    merges in the same order, so a pruned copy (prune_tree) keeps the browse tree's topology and
    leaf order. The average distance of two clusters is summed in the JS's order (the earlier
    cluster's indices outside), so ties break the same way; a pair's average never changes, so
    each is computed once."""
    clusters = [{'name': label, 'id': i, 'distance': 0.0, 'children': None, 'indices': [i], 'size': 1}
                for i, label in enumerate(labels)]
    averages = {}   # (id of the earlier cluster, id of the later one) -> average distance

    def average(first, second):
        key = (first['id'], second['id'])
        if key not in averages:
            total, count = 0.0, 0
            for i in first['indices']:
                row = matrix[i]
                for j in second['indices']:
                    total += row[j]
                    count += 1
            averages[key] = total / count if count else 0
        return averages[key]

    next_id = len(labels)
    while len(clusters) > 1:
        best, best_i, best_j = float('inf'), 0, 1
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                d = average(clusters[i], clusters[j])
                if d < best:
                    best, best_i, best_j = d, i, j
        first, second = clusters[best_i], clusters[best_j]
        merged = {'name': f'internal_{next_id}', 'id': next_id, 'distance': best,
                  'children': [first, second], 'indices': first['indices'] + second['indices'],
                  'size': first['size'] + second['size']}
        next_id += 1
        clusters = [c for k, c in enumerate(clusters) if k not in (best_i, best_j)] + [merged]

    def balance(node):   # optimizeTree: the smaller subtree first ...
        if not node['children']:
            return node['size'] or 1
        left, right = balance(node['children'][0]), balance(node['children'][1])
        if left > right:
            node['children'].reverse()
        node['size'] = left + right
        return node['size']

    def reorder(node):   # ... then the lower merge first
        if not node['children']:
            return
        reorder(node['children'][0])
        reorder(node['children'][1])
        if node['children'][0]['distance'] > node['children'][1]['distance']:
            node['children'].reverse()

    def plain(node):
        if not node['children']:
            return {'name': node['name']}
        return {'distance': node['distance'], 'children': [plain(c) for c in node['children']]}

    if not clusters:
        return None
    root = clusters[0]
    balance(root)
    reorder(root)
    return plain(root)


def prune_tree(node, keep):
    """`node` (average_linkage_tree) with only the leaves named in `keep`: an internal node
    left with one child is replaced by it, the others keep their merge distance."""
    if node is None:
        return None
    if 'children' not in node:
        return dict(node) if node['name'] in keep else None
    children = [c for c in (prune_tree(child, keep) for child in node['children']) if c]
    if len(children) < 2:
        return children[0] if children else None
    return {'distance': node['distance'], 'children': children}


class DendrogramHandler:
    def __init__(self, seq_id=[]):
        """

        Args:
            seq_id: List of chosen IDs, default is empty.
        """
        self.proteins = Protein.objects.all()
        self.seq_id = seq_id 

        self.colors = {
            'primary': '#569976',
            'secondary': '#417459', 
            'accent': '#407754',
            'background': '#f9f9f9',
            'light_green': '#eef6f0',
            'border': '#ddd',
            'text': '#333',
            'white': '#ffffff'
        }
        self.clade_colors = dict(CLADE_COLORS)
        
    def distance_matrix(self) -> pd.DataFrame:
        """The pipeline's distance matrix for the proteins in the DB, labelled by protein code
        (a structure shared by several proteins gives each of them a row)."""
        if not DISTANCE_MATRIX_CSV.exists():
            print(f"Warning: {DISTANCE_MATRIX_CSV} not found - run pipeline step 9 (distance_matrix)")
            return pd.DataFrame()
        df = pd.read_csv(DISTANCE_MATRIX_CSV, index_col=0)
        codes, structs = [], []
        for code, struct_id in self.proteins.values_list('protein_code', 'selected_structure'):
            if struct_id in df.index:
                codes.append(code)
                structs.append(struct_id)
        matrix = df.loc[structs, structs]
        matrix.index = codes
        matrix.columns = codes
        return matrix

    def get_protein_info(self, protein_code):
        try:
            protein = Protein.objects.get(protein_code=protein_code)
            substrates = ", ".join([s.name for s in protein.substrate.all()]) or "Unknown"
            mass_kda = polypeptide_kda(protein.sequence)   # polypeptide mass of the structure sequence
            from .views import _catalytic_domain_length   # same catalytic length as the results table
            from .overview import _lineage, lineage_text, regioselectivity_label  # lineage (step 2i), C1/C4
            catalytic_length = _catalytic_domain_length(protein)
            
            # Calculate sequence length (amino acids, excluding gaps)
            # Use aligned_sequence if available, otherwise use sequence
            # Remove gaps and whitespace to get actual amino acid count
            sequence_length = 0
            seq_to_use = None
            if hasattr(protein, 'aligned_sequence') and protein.aligned_sequence:
                seq_to_use = protein.aligned_sequence
            elif protein.sequence:
                seq_to_use = protein.sequence
            
            if seq_to_use:
                # Remove gaps (-), whitespace, and newlines to get actual amino acid sequence length
                sequence_clean = seq_to_use.replace('-', '').replace(' ', '').replace('\n', '').replace('\r', '').replace('\t', '')
                # Count only valid amino acid letters (A-Z, excluding gaps and special chars)
                sequence_clean = ''.join([c for c in sequence_clean if c.isalpha()])
                sequence_length = len(sequence_clean)
            
            # Extract total charge from PQR file
            total_charge = None
            try:
                if hasattr(protein, 'electrostatics') and protein.electrostatics:
                    pqr_content = protein.electrostatics.pqr_file
                    if pqr_content:
                        pqr_lines = pqr_content.splitlines()
                        # Extract total charge from REMARK line
                        import re
                        for line in pqr_lines:
                            if 'REMARK' in line and 'Total charge' in line:
                                try:
                                    # Extract charge value using regex: "Total charge on this protein: -8.0000 e"
                                    # Look for number (with optional sign and decimal) followed by 'e'
                                    match = re.search(r'([-+]?\d+\.?\d*)\s+e', line)
                                    if match:
                                        total_charge = float(match.group(1))
                                    else:
                                        # Fallback: try splitting by colon
                                        if ':' in line:
                                            charge_part = line.split(':')[1].strip()
                                            charge_str = charge_part.split('e')[0].strip()
                                            total_charge = float(charge_str)
                                except (ValueError, IndexError, AttributeError):
                                    pass
                                break
            except Exception as e:
                print(f"Warning: Could not extract charge for {protein_code}: {e}")
            
            return {
                'code': protein.protein_code,
                'other_names': protein.other_names or '',
                'clade': protein.clade.name,
                'organism': protein.organism.name,
                'kingdom': protein.organism.kingdom,
                'taxonomy': lineage_text(_lineage(protein.organism)),
                'substrates': substrates,
                'ec_numbers': protein.ec_numbers or '',
                'regioselectivity': regioselectivity_label(r.name for r in protein.regioselectivity.all()),
                'selected_structure': protein.selected_structure,
                'sequence_length': sequence_length,
                'catalytic_length': catalytic_length,
                'mass_kda': round(mass_kda, 1) if mass_kda else None,
                'associated_CBMs': protein.associated_CBMs or '',
                'uniprot': protein.uniprot or '',
                'uniparc': getattr(protein, 'uniparc', None) or '',
                'total_charge': total_charge
            }
        except Protein.DoesNotExist:
            return {
                'code': protein_code,
                'name': 'Unknown',
                'clade': 'Unknown',
                'organism': 'Unknown', 
                'substrates': 'Unknown',
                'regioselectivity': 'Unknown',
                'selected_structure': 'Unknown',
                'associated_CBMs': '',
                'sequence_length': 0,
                'total_charge': None
            }

    def static_dendrogram_circular(self):
        df = self.distance_matrix()
        full_matrix = df.values.tolist()

        if len(self.seq_id) == 0:
            labels = list(df.index)
            title = f"LPMO Protein Phylogenetic Tree ({len(labels)} proteins)"
        else:
            selected_indices = [i for i, label in enumerate(df.index) if label in self.seq_id]
            labels = [df.index[i] for i in selected_indices]
            #filtered matrix
            filtered_matrix = []
            for i in selected_indices:
                row = []
                for j in selected_indices:
                    row.append(full_matrix[i][j])
                filtered_matrix.append(row)
            full_matrix = filtered_matrix
            title = f"LPMO Protein Tree - {len(self.seq_id)} selected"

        protein_info = {}
        for label in labels:
            protein_info[label] = self.get_protein_info(label)

        # clade legend with the number of proteins in the tree; circular_dendrogram.js
        # (bindLegend) makes each clade highlight and select its proteins
        clade_counts = Counter(info.get('clade') for info in protein_info.values())
        legend_items = ''.join(
            f'<button type="button" class="tree-legend-item" data-family="{escape(clade)}">'
            f'<span class="tree-legend-swatch" style="background: {color};"></span>'
            f'<span class="tree-legend-name">{escape(clade)}</span>'
            f'<span class="tree-legend-count">{clade_counts[clade]}</span></button>'
            for clade, color in self.clade_colors.items() if clade_counts.get(clade))

        import json
        
        html_content = f"""
        <div id="circular-dendrogram-container" style="width: 100%; height: 100%; min-height: 700px; position: relative;">
            <div id="circular-dendrogram" style="width: 100%; height: 100%;"></div>
            
            <!-- Branch lengths and Export SVG: in the browse page's tree toolbar (index.html) -->
            <!-- Legend -->
            <div class="tree-legend" id="tree-legend" role="group" aria-label="LPMO clades">
                <div class="tree-legend-items">{legend_items}</div>
            </div>
        </div>

        <script src="https://d3js.org/d3.v7.min.js"></script>
        <script src="/static/catalogue/js/hierarchical_clustering.js"></script>
        <script src="/static/catalogue/js/circular_dendrogram.js"></script>
        
        <style>
        .dendrogram-tooltip {{
            position: absolute;
            visibility: hidden;
            background: rgba(0, 0, 0, 0.8);
            color: white;
            padding: 10px;
            border-radius: 5px;
            font-size: 12px;
            max-width: 300px;
            z-index: 1000;
            pointer-events: none;
        }}
        
        .zoom-control:hover rect {{
            fill: #e0e0e0 !important;
        }}
        
        .link {{
            cursor: pointer;
        }}
        
        .link:hover {{
            stroke-width: 2.5px !important;
        }}
        
        .label:hover {{
            font-weight: bold;
        }}
        
        .node:hover {{
            stroke-width: 2px !important;
            r: 4px !important;
        }}
        </style>

        <script>
        (function(initFn) {{
            if (document.readyState === 'loading') {{
                document.addEventListener('DOMContentLoaded', initFn);
            }} else {{
                initFn();
            }}
        }})(function() {{
            const distanceMatrix = {json.dumps(full_matrix)};
            const labels = {json.dumps(labels)};
            const proteinInfo = {json.dumps(protein_info)};
            
            const getSelectedProteins = () => {{
                if (typeof window.getEffectiveSelectedProteins === 'function') {{
                    try {{
                        return window.getEffectiveSelectedProteins();
                    }} catch (e) {{}}
                }}
                const selected = new Set();
                document.querySelectorAll('input[name="protein_ids"]:checked').forEach(cb => selected.add(cb.value));
                return Array.from(selected);
            }};
            const selectedProteins = getSelectedProteins();
            
            const container = document.getElementById('circular-dendrogram');
            if (!container) {{
                console.error('Circular dendrogram container not found');
                return;
            }}
            
            const containerRect = container.getBoundingClientRect();
            const size = Math.min(Math.max(containerRect.width, 600), Math.max(containerRect.height, 700));
            
            const branchLengths = document.getElementById('toggle-branch-lengths');
            const treeLabels = document.getElementById('toggle-tree-labels');
            const dendrogram = new CircularDendrogram('circular-dendrogram', {{
                width: size,
                height: size,
                showBranchLengths: !!(branchLengths && branchLengths.checked),
                showLabels: !treeLabels || treeLabels.checked,
                cladeColors: {json.dumps(self.clade_colors)}
            }});
            
            window.currentDendrogram = dendrogram;
            
            dendrogram.setData(distanceMatrix, labels);
            dendrogram.render(distanceMatrix, labels, proteinInfo, selectedProteins);
            dendrogram.bindLegend(document.getElementById('tree-legend'));

            console.log('Circular dendrogram initialized successfully');
        }});
        
        if (typeof window !== 'undefined' && typeof CircularDendrogram !== 'undefined') {{
            window.initializeDendrogramFallback = function() {{
                console.log('Using fallback dendrogram initialization...');
                const container = document.getElementById('circular-dendrogram');
                if (container && !window.currentDendrogram) {{
                    try {{
                        const dendrogram = new CircularDendrogram('circular-dendrogram', {{
                            width: 800,
                            height: 800
                        }});
                        window.currentDendrogram = dendrogram;
                        console.log('Fallback dendrogram created');
                    }} catch (e) {{
                        console.error('Fallback dendrogram creation failed:', e);
                    }}
                }}
            }};
            
            // Try to initialize after a delay
            setTimeout(window.initializeDendrogramFallback, 500);
        }}
        </script>
        """
        
        return html_content

    def get_clade_legend(self):
        clades_in_db = set(self.proteins.values_list('clade__name', flat=True))
        
        legend_html = f"<div style='margin-top: 20px; padding: 15px; background: {self.colors['background']}; border-radius: 8px; border: 1px solid {self.colors['border']};'>"
        legend_html += f"<h4 style='margin: 0 0 10px 0; color: {self.colors['text']};'>Clade Color Legend:</h4>"
        legend_html += "<div style='display: flex; flex-wrap: wrap; gap: 10px;'>"
        
        for clade in sorted(clades_in_db):
            color = self.clade_colors.get(clade, self.colors['secondary'])
            legend_html += f"""
            <div style='display: flex; align-items: center; gap: 5px;'>
                <div style='width: 16px; height: 16px; background: {color}; border-radius: 3px;'></div>
                <span style='font-size: 12px; color: {self.colors['text']};'>{clade}</span>
            </div>
            """
        
        legend_html += "</div></div>"
        return legend_html
