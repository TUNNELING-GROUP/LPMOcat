/* Compare page, selection overview: the browse page's tree pruned to the compared proteins
   (views.compare_overview, dendrogram.prune_tree), drawn left to right with the merge distance
   on the x axis, and a few features of each protein in columns beside its label.
   catalogue_overview.js draws it with the overview data (SelectionTree.draw) and, while a
   chart's bar or dot is hovered, highlights its proteins here (setOverviewHighlight). Hovering
   a row or a branch shows its proteins; clicking a row picks that protein and clicking a branch
   its whole subtree (again: unpicks), "Which features discriminate the picked proteins from …"
   (grouping_features.js; picked(), onPickChange()). Export SVG saves the drawing. */
(function () {
    'use strict';

    const SVG_NS = 'http://www.w3.org/2000/svg';
    const ROW = 20;            // px per protein
    const TREE_WIDTH = 170;    // px of the branches
    const HEAD = 34;           // px of the column heads
    const AXIS = 46;           // px of the distance axis under the tree
    const GAP = 14;            // px between columns
    const FONT = 'Arial, sans-serif';
    const INK = '#333333';
    const MUTED = '#8a8a8a';
    const BRANCH = '#9aa5a0';
    const MAX_CODES = 14;

    // ---- cutting the tree into groups and shading them ----

    // Maximal subtrees whose merge distance is at or below `cut`, left to right. A leaf that
    // joins the tree only above the cut is its own group.
    function clustersAt(node, cut) {
        if (!node.children || node.distance <= cut) return [node];
        return node.children.flatMap(child => clustersAt(child, cut));
    }

    function hexToHsl(hex) {
        const m = /^#?([0-9a-f]{6})$/i.exec(hex || '');
        if (!m) return null;
        const n = parseInt(m[1], 16);
        const r = ((n >> 16) & 255) / 255, g = ((n >> 8) & 255) / 255, b = (n & 255) / 255;
        const max = Math.max(r, g, b), min = Math.min(r, g, b), l = (max + min) / 2;
        if (max === min) return { h: 0, s: 0, l: l * 100 };
        const d = max - min;
        const s = d / (1 - Math.abs(2 * l - 1));
        let h;
        if (max === r) h = ((g - b) / d) % 6;
        else if (max === g) h = (b - r) / d + 2;
        else h = (r - g) / d + 4;
        return { h: (h * 60 + 360) % 360, s: s * 100, l: l * 100 };
    }

    // Shades of one family colour for the groups that belong to it: the hue stays, the
    // lightness spreads around the family's own, so the family is still readable.
    function shadesOf(hex, count) {
        const hsl = hexToHsl(hex);
        if (!hsl || count <= 1) return Array(Math.max(count, 1)).fill(hex);
        const span = Math.min(19, 46 / count) * (count - 1) / 2;
        return Array.from({ length: count }, (_, i) => {
            const t = count === 1 ? 0 : (i / (count - 1)) * 2 - 1;      // -1 .. 1
            const l = Math.min(74, Math.max(24, hsl.l + t * span));
            return `hsl(${hsl.h.toFixed(0)} ${Math.max(18, hsl.s).toFixed(0)}% ${l.toFixed(0)}%)`;
        });
    }

    // text width at the drawing's font size (a canvas measures without layout)
    const context = document.createElement('canvas').getContext('2d');
    function textWidth(text, { size = 11, bold = false, italic = false } = {}) {
        context.font = `${italic ? 'italic ' : ''}${bold ? 'bold ' : ''}${size}px Arial, sans-serif`;
        return context.measureText(text).width;
    }

    function svg(tag, attrs = {}, parent = null) {
        const node = document.createElementNS(SVG_NS, tag);
        Object.entries(attrs).forEach(([key, value]) => node.setAttribute(key, value));
        if (parent) parent.appendChild(node);
        return node;
    }

    function el(tag, className, text) {
        const node = document.createElement(tag);
        if (className) node.className = className;
        if (text !== undefined) node.textContent = text;
        return node;
    }

    const leavesOf = node => (node.children ? node.children.flatMap(leavesOf) : [node.name]);

    // Nice tick step for a distance axis of `span`
    function tickStep(span) {
        const raw = span / 4;
        const power = 10 ** Math.floor(Math.log10(raw));
        return [1, 2, 5, 10].map(m => m * power).find(step => step >= raw);
    }

    const STRUCTURE_NAMES = { 'AlphaFold DB': 'AlphaFold', 'own model': 'AlphaFold 3', CAZyme3D: 'CAZyme3D' };

    // The columns beside the labels: text(record) and, optionally, muted(record) (grey text:
    // a model, a regioselectivity from the EC numbers only, nothing found) and italic
    function columns(data) {
        const notDetermined = value => !value || value === data.not_determined;
        return [
            {
                head: 'Structure',
                text: r => (r.source === 'PDB' ? r.structure : STRUCTURE_NAMES[r.source] || r.source),
                muted: r => r.source !== 'PDB',
            },
            {
                head: 'Regioselectivity',
                text: r => (r.regio[0] === 'Unknown' ? '–' : r.regio[0]),
                muted: r => r.regio[1] !== 'papers',
            },
            {
                head: 'Second sphere',
                text: r => (notDetermined(r.motif[0]) ? '–' : r.motif[0]),
                muted: r => notDetermined(r.motif[0]),
            },
            {
                head: 'Shielding',
                text: r => (notDetermined(r.distal) ? '–' : r.distal === data.open_position ? 'open' : r.distal),
                muted: r => notDetermined(r.distal) || r.distal === data.open_position,
            },
            {
                head: 'N-His',
                text: r => (r.brace.n_his && r.brace.n_his[0]) || '–',
                muted: r => !(r.brace.n_his && r.brace.n_his[0]) || r.brace.n_his[0] === data.other_motif,
            },
            {
                head: 'Internal-His',
                text: r => (r.brace.internal_his && r.brace.internal_his[0]) || '–',
                muted: r => !(r.brace.internal_his && r.brace.internal_his[0])
                    || r.brace.internal_his[0] === data.other_motif,
            },
            { head: 'Species', text: r => r.species, italic: true },
        ];
    }

    const REGIO_SOURCES = { papers: 'reported in papers', ec: 'from the EC numbers', none: '' };

    // tooltip lines of one protein
    function proteinContent(r, data, group = null) {
        const regio = r.regio[0] === 'Unknown' ? 'unknown' : r.regio[0];
        const regioSource = REGIO_SOURCES[r.regio[1]] || '';
        const reported = new Set(r.substrate_papers || []);
        const substrates = r.substrates.length
            ? r.substrates.map(s => (reported.has(s) ? s : `${s} (EC/name)`)).join(', ') : 'none assigned';
        const brace = ['n_his', 'internal_his']
            .map(which => r.brace[which] && r.brace[which][1])
            .filter(Boolean).join(' · ');
        return {
            heading: `${r.family} · ${r.species}`,
            value: r.code,
            lines: [
                group ? `Group: ${group.name}` : '',
                `Structure: ${r.structure}${r.source === 'PDB' ? '' : ` (${r.source === 'own model' ? 'AlphaFold 3 model' : r.source + ' model'})`}`,
                `Regioselectivity: ${regio}${regioSource ? ` (${regioSource})` : ''}`,
                `Substrates: ${substrates}`,
                r.motif[1] ? `Second-sphere motif: ${r.motif[1]}` : '',
                r.distal && r.distal !== data.not_determined ? `Shielding residue: ${r.distal}` : '',
                brace ? `Brace: ${brace}` : '',
                r.cbms.length ? `CBMs: ${r.cbms.join(', ')}` : '',
            ].filter(Boolean),
        };
    }

    // groups: {code: {name, color}} of a comparison with groups (browse page, "Add as group"):
    // a first column with the group's colour bar and name
    function draw(container, { tree, data, familyColor, missing = [], groups = null }) {
        container.replaceChildren();
        const records = new Map(data.proteins.map(p => [p.code, p]));
        const groupOf = code => (groups && groups[code]) || null;
        if (!tree) {
            container.appendChild(el('p', 'ov-hint', 'None of the compared proteins is in the tree.'));
            return null;
        }
        const leaves = leavesOf(tree);
        const cols = columns(data);
        const BAR = 8;   // the group column's colour bar and the gap after it
        if (groups) {
            cols.unshift({
                head: 'Group',
                text: r => (groupOf(r.code) ? groupOf(r.code).name : ''),
                color: r => (groupOf(r.code) ? groupOf(r.code).color : null),
                muted: r => !!(groupOf(r.code) && groupOf(r.code).ungrouped),
            });
        }

        // x of a merge distance (the root at the left, the proteins at TREE_WIDTH); row y of a leaf
        const maxDistance = tree.children ? tree.distance : 0;
        const x = d => (maxDistance ? TREE_WIDTH * (1 - d / maxDistance) : TREE_WIDTH);
        const rowOf = new Map(leaves.map((code, i) => [code, i]));
        const y = i => HEAD + i * ROW + ROW / 2;

        // column positions from the widest entry (or head)
        const labelX = TREE_WIDTH + 16;
        const labelWidth = Math.max(...leaves.map(code => textWidth(code, { bold: true })));
        let left = labelX + labelWidth + GAP;
        cols.forEach(col => {
            col.x = left;
            col.width = Math.max(textWidth(col.head, { size: 10.5, bold: true }),
                ...leaves.map(code => textWidth(col.text(records.get(code)), { italic: col.italic })
                    + (col.color ? BAR : 0)));
            left += col.width + GAP;
        });
        const width = Math.ceil(left - GAP + 4);
        const height = HEAD + leaves.length * ROW + AXIS;

        const root = svg('svg', {
            width, height, viewBox: `0 0 ${width} ${height}`, 'font-family': FONT, 'font-size': 11,
            role: 'img', 'aria-label': `Tree of the ${leaves.length} compared proteins`,
        });
        root.classList.add('sel-tree-svg');
        svg('rect', { x: 0, y: 0, width, height, fill: 'white' }, root);

        // row bands (hover and highlight) under everything else
        const bands = svg('g', {}, root);
        const bandOf = new Map();
        leaves.forEach((code, i) => {
            const band = svg('rect', {
                x: labelX - 14, y: HEAD + i * ROW, width: width - labelX + 14, height: ROW, fill: 'transparent',
            }, bands);
            band.classList.add('sel-tree-band');
            band.dataset.code = code;
            bandOf.set(code, band);
        });

        // column heads
        const heads = svg('g', { 'font-size': 10.5, 'font-weight': 'bold', fill: '#555' }, root);
        svg('text', { x: labelX, y: HEAD - 10 }, heads).textContent = 'Protein';
        cols.forEach(col => {
            svg('text', { x: col.x, y: HEAD - 10 }, heads).textContent = col.head;
        });
        svg('line', {
            x1: 0, x2: width, y1: HEAD - 4.5, y2: HEAD - 4.5, stroke: '#e3e3e3', 'stroke-width': 1,
        }, root);

        // branches: each internal node is drawn as its vertical line and the horizontal lines to
        // its children; a node's hit area covers them, so hovering a branch shows its subtree
        const branches = svg('g', { stroke: BRANCH, 'stroke-width': 1.3, fill: 'none' }, root);
        const hits = svg('g', { stroke: 'transparent', 'stroke-width': 9, fill: 'none' }, root);
        function place(node) {
            if (!node.children) {
                node.y = y(rowOf.get(node.name));
                node.x = TREE_WIDTH;
                return;
            }
            node.children.forEach(place);
            node.x = x(node.distance);
            node.y = (node.children[0].y + node.children[node.children.length - 1].y) / 2;
            const ys = node.children.map(c => c.y);
            const d = `M${node.x},${Math.min(...ys)}V${Math.max(...ys)}`
                + node.children.map(c => `M${node.x},${c.y}H${c.x}`).join('');
            node.__path = svg('path', { d }, branches);
            const hit = svg('path', { d }, hits);
            hit.classList.add('sel-tree-branch');
            hit.__node = node;
        }
        place(tree);
        if (tree.children) {   // the root's stub
            svg('path', { d: `M${tree.x - 8},${tree.y}H${tree.x}` }, branches);
        }

        // proteins: family dot, code, columns
        const dotOf = new Map();
        const labels = svg('g', {}, root);
        leaves.forEach((code, i) => {
            const r = records.get(code);
            const row = svg('g', { transform: `translate(0,${y(i)})` }, labels);
            svg('line', {
                x1: TREE_WIDTH, x2: labelX - 12, y1: 0, y2: 0, stroke: '#dde3e0', 'stroke-dasharray': '2 2',
            }, row);
            dotOf.set(code, svg('circle', {
                cx: labelX - 7, cy: 0, r: 4.5, fill: familyColor[r.family] || '#569976', stroke: 'white', 'stroke-width': 1,
            }, row));
            const label = svg('text', { x: labelX, y: 0, dy: '0.35em', fill: INK, 'font-weight': 'bold' }, row);
            label.textContent = code;
            cols.forEach(col => {
                const color = col.color && col.color(r);
                if (color) {
                    svg('rect', { x: col.x, y: -ROW / 2 + 3, width: 4, height: ROW - 6, rx: 1, fill: color }, row);
                }
                const text = svg('text', {
                    x: col.x + (col.color ? BAR : 0), y: 0, dy: '0.35em', fill: col.muted && col.muted(r) ? MUTED : INK,
                    ...(col.italic ? { 'font-style': 'italic' } : {}),
                }, row);
                text.textContent = col.text(r);
            });
        });

        // distance axis under the branches
        if (maxDistance) {
            const axisY = HEAD + leaves.length * ROW + 8;
            const axis = svg('g', { fill: MUTED, 'font-size': 10 }, root);
            svg('line', { x1: 0, x2: TREE_WIDTH, y1: axisY, y2: axisY, stroke: '#cccccc' }, axis);
            const step = tickStep(maxDistance);
            const decimals = Math.max(0, -Math.floor(Math.log10(step)));
            for (let k = 0; k * step <= maxDistance + 1e-9; k++) {
                const tx = x(k * step);
                svg('line', { x1: tx, x2: tx, y1: axisY, y2: axisY + 4, stroke: '#cccccc' }, axis);
                svg('text', { x: tx, y: axisY + 14, 'text-anchor': 'middle' }, axis)
                    .textContent = (k * step).toFixed(decimals);
            }
            svg('text', { x: TREE_WIDTH / 2, y: axisY + 27, 'text-anchor': 'middle' }, axis)
                .textContent = 'Sequence distance (BLOSUM62)';
        }

        container.appendChild(root);
        if (missing.length) {
            container.appendChild(el('p', 'ov-hint',
                `Not in the tree (no distance to the others): ${missing.join(', ')}`));
        }

        // ---- hover: a row or a branch ----
        const tooltip = el('div', 'ov-tooltip');
        tooltip.hidden = true;
        tooltip.setAttribute('aria-hidden', 'true');
        document.body.appendChild(tooltip);

        function showTooltip(event, content) {
            tooltip.replaceChildren();
            if (content.heading) tooltip.appendChild(el('div', 'ov-tip-heading', content.heading));
            tooltip.appendChild(el('div', 'ov-tip-value', content.value));
            (content.rows || []).forEach(row => {
                const line = el('div', 'ov-tip-row');
                const key = el('span', 'ov-tip-key');
                key.style.background = row.color;
                line.append(key, el('span', 'ov-tip-row-label', row.label));
                if (row.value !== undefined) line.appendChild(el('span', 'ov-tip-row-value', String(row.value)));
                tooltip.appendChild(line);
            });
            (content.lines || []).forEach(text => tooltip.appendChild(el('div', 'sel-tip-line', text)));
            if (content.codes && content.codes.length) {
                const more = content.codes.length - MAX_CODES;
                const listed = content.codes.slice(0, MAX_CODES).join(', ');
                tooltip.appendChild(el('div', 'ov-tip-codes', more > 0 ? `${listed} … and ${more} more` : listed));
            }
            tooltip.hidden = false;
            moveTooltip(event);
        }
        function moveTooltip(event) {
            if (tooltip.hidden) return;
            const pad = 14;
            const { offsetWidth: w, offsetHeight: h } = tooltip;
            let tx = event.clientX + pad;
            let ty = event.clientY + pad;
            if (tx + w > window.innerWidth - 8) tx = event.clientX - pad - w;
            if (ty + h > window.innerHeight - 8) ty = event.clientY - pad - h;
            tooltip.style.left = `${Math.max(8, tx)}px`;
            tooltip.style.top = `${Math.max(8, ty)}px`;
        }
        const hideTooltip = () => { tooltip.hidden = true; };

        let overview = new Set();   // highlighted from a chart
        let hovered = new Set();    // highlighted by hovering here
        const picked = new Set();   // picked by clicking a row or a branch
        const pickListeners = [];
        function paint() {
            leaves.forEach(code => {
                const on = overview.has(code) || hovered.has(code);
                bandOf.get(code).classList.toggle('sel-tree-hl', on);
                bandOf.get(code).classList.toggle('sel-tree-picked', picked.has(code));
            });
        }
        // a click picks the proteins, or unpicks them when all of them are picked already
        function togglePick(codes) {
            const all = codes.every(code => picked.has(code));
            codes.forEach(code => (all ? picked.delete(code) : picked.add(code)));
            paint();
            const list = leaves.filter(code => picked.has(code));
            pickListeners.forEach(fn => fn(list));
        }

        function subtreeContent(node) {
            const codes = leavesOf(node);
            const byFamily = new Map();
            codes.forEach(code => {
                const f = records.get(code).family;
                byFamily.set(f, (byFamily.get(f) || 0) + 1);
            });
            return {
                heading: `Merged at distance ${node.distance.toFixed(3)}`,
                value: `${codes.length} proteins`,
                rows: [...byFamily].map(([f, n]) => ({ color: familyColor[f], label: f, value: n })),
                codes,
            };
        }

        // a branch shows its subtree, a row (its band: the row's text lets the pointer through)
        // its protein; anything else ends the hover
        root.addEventListener('mouseover', event => {
            const target = event.target;
            const band = target.closest('.sel-tree-band');
            if (target.classList.contains('sel-tree-branch')) {
                hovered = new Set(leavesOf(target.__node));
                showTooltip(event, subtreeContent(target.__node));
            } else if (band) {
                hovered = new Set([band.dataset.code]);
                showTooltip(event, proteinContent(records.get(band.dataset.code), data, groupOf(band.dataset.code)));
            } else {
                hovered = new Set();
                hideTooltip();
            }
            paint();
        });
        root.addEventListener('mousemove', moveTooltip);
        root.addEventListener('click', event => {
            const target = event.target;
            const band = target.closest('.sel-tree-band');
            if (target.classList.contains('sel-tree-branch')) togglePick(leavesOf(target.__node));
            else if (band) togglePick([band.dataset.code]);
        });
        root.addEventListener('mouseleave', () => {
            hovered = new Set();
            hideTooltip();
            paint();
        });
        labels.style.pointerEvents = 'none';

        // ---- groups from a distance cut ----
        // Drawn above the branches so the dashed line stays visible over them.
        const cutLine = svg('line', {
            x1: 0, x2: 0, y1: HEAD - 4, y2: HEAD + leaves.length * ROW + 4,
            stroke: '#d33682', 'stroke-width': 1.4, 'stroke-dasharray': '4 3', display: 'none',
        }, root);
        cutLine.classList.add('sel-tree-cut');

        // The groups of the last cut, [{codes, color}] in tree order, or null without one
        // (read by grouping_features.js when the user asks what separates them).
        let currentGroups = null, currentCut = null;

        function setCutoff(cut) {
            // cut >= the root's distance (or null) means no grouping: back to family colours.
            if (cut == null || !maxDistance || cut >= maxDistance) {
                cutLine.setAttribute('display', 'none');
                const reset = node => {
                    if (node.__path) node.__path.removeAttribute('stroke');
                    (node.children || []).forEach(reset);
                };
                reset(tree);
                leaves.forEach(code => {
                    const r = records.get(code);
                    dotOf.get(code).setAttribute('fill', familyColor[r.family] || '#569976');
                });
                currentGroups = null;
                currentCut = null;
                return { groups: 1, cut: null };
            }
            cutLine.setAttribute('x1', x(cut));
            cutLine.setAttribute('x2', x(cut));
            cutLine.setAttribute('display', '');

            const groups = clustersAt(tree, cut);
            // A group takes the family most of its proteins belong to; each family then
            // spreads its groups over shades of its own colour, in tree order.
            const familyOf = groups.map(node => {
                const counts = new Map();
                leavesOf(node).forEach(code => {
                    const f = (records.get(code) || {}).family;
                    counts.set(f, (counts.get(f) || 0) + 1);
                });
                return [...counts.entries()].sort((a, b) => b[1] - a[1])[0][0];
            });
            const perFamily = new Map();
            familyOf.forEach(f => perFamily.set(f, (perFamily.get(f) || 0) + 1));
            const shades = new Map();
            perFamily.forEach((count, f) => shades.set(f, shadesOf(familyColor[f] || '#569976', count)));
            const used = new Map();
            const groupList = [];

            groups.forEach((node, i) => {
                const f = familyOf[i];
                const k = used.get(f) || 0;
                used.set(f, k + 1);
                const color = shades.get(f)[k];
                groupList.push({ codes: leavesOf(node), color });
                const paintNode = n => {
                    if (n.__path) n.__path.setAttribute('stroke', color);
                    (n.children || []).forEach(paintNode);
                };
                paintNode(node);
                leavesOf(node).forEach(code => dotOf.get(code).setAttribute('fill', color));
            });
            // Branches above the cut keep the default grey.
            const clear = node => {
                if (groups.includes(node)) return;
                if (node.__path) node.__path.removeAttribute('stroke');
                (node.children || []).forEach(clear);
            };
            clear(tree);
            currentGroups = groupList;
            currentCut = cut;
            return { groups: groups.length, cut };
        }

        return {
            maxDistance,
            setCutoff,
            groups: () => (currentGroups ? { cut: currentCut, groups: currentGroups } : null),
            // The tree distance of every pair of proteins: the height at which they join (their
            // lowest common branch), in leaf order. grouping_features.js correlates it with the
            // proteins' feature differences ("Which features follow the tree?").
            cophenetic() {
                const index = new Map(leaves.map((code, i) => [code, i]));
                const matrix = leaves.map(() => leaves.map(() => 0));
                const walk = node => {
                    if (!node.children) return [node.name];
                    const parts = node.children.map(walk);
                    for (let a = 0; a < parts.length; a++) {
                        for (let b = a + 1; b < parts.length; b++) {
                            parts[a].forEach(x => parts[b].forEach(y => {
                                matrix[index.get(x)][index.get(y)] = node.distance;
                                matrix[index.get(y)][index.get(x)] = node.distance;
                            }));
                        }
                    }
                    return parts.flat();
                };
                walk(tree);
                return { codes: leaves.slice(), matrix };
            },
            leaves: () => leaves.slice(),   // the proteins in row order, top to bottom
            picked: () => leaves.filter(code => picked.has(code)),
            clearPicked() { togglePick(leaves.filter(code => picked.has(code))); },
            onPickChange(fn) { pickListeners.push(fn); },
            setOverviewHighlight(codes) {
                overview = new Set(codes || []);
                paint();
            },
            exportSvg(filename = 'lpmocat_selection_tree.svg') {
                const copy = root.cloneNode(true);
                copy.querySelectorAll('.sel-tree-branch').forEach(node => node.remove());
                copy.querySelectorAll('.sel-tree-band').forEach(node => node.setAttribute('fill', 'transparent'));
                copy.setAttribute('xmlns', SVG_NS);
                const blob = new Blob([new XMLSerializer().serializeToString(copy)], { type: 'image/svg+xml' });
                const link = document.createElement('a');
                link.href = URL.createObjectURL(blob);
                link.download = filename;
                document.body.appendChild(link);
                link.click();
                link.remove();
                setTimeout(() => URL.revokeObjectURL(link.href), 1000);
            },
        };
    }

    window.SelectionTree = { draw };
})();
