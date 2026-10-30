/* Browse page: the "Catalogue overview" | "Selection criteria" tabs of the left panel, and the
   overview (_catalogue_overview.html): headline counts and small Plotly charts of the
   per-protein records from overview.py. Hovering a bar or dot lists its proteins
   and highlights them in the dendrogram beside it (CircularDendrogram.setOverviewHighlight);
   clicking selects them (or deselects them when all are selected), like a branch click in the
   tree; the category labels of the bar charts act like their bars ("Fungi" stands for all its
   phyla). Every chart has a table view with the same numbers.
   Compare page (_selection_overview.html, data-mode="selection"): the same charts count the
   compared proteins, with the rest of the catalogue in grey behind them (or not, by the card's
   switch); the data come from views.compare_overview after the page has loaded, hovering
   highlights proteins in the selection's tree (selection_tree.js), and clicks select nothing. */
(async function () {
    'use strict';

    const root = document.getElementById('catalogueOverview');
    if (!root) return;
    const selectionMode = root.dataset.mode === 'selection';

    const INK = '#333333';
    const INK_MUTED = '#666666';
    const SELECTED_RING = '#1a1a1a';
    const GRID = '#ececec';
    const AXIS_LINE = '#cccccc';
    const SURFACE = '#ffffff';
    const FONT = { family: 'Arial, sans-serif', size: 11, color: INK_MUTED };
    const PREDICTED_ALPHA = '59';   // hex alpha (35%) of a family colour for predicted models
    // Model sources (overview.structure_source) and how the tiles and tooltips name them
    const MODEL_SOURCES = [['AlphaFold DB', 'AlphaFold DB'], ['own model', 'own AlphaFold 3'],
                           ['CAZyme3D', 'CAZyme3D']];
    const PLOT_HEIGHTS = { normal: 230, tall: 320, wide: 340, xtall: 420 };   // px, as .ov-plot in the CSS
    const AXIS_BAND = 45;           // px of a plot taken by the x axis and its title
    const MAX_BAR = 22;             // px, bar thickness cap
    const MAX_CODES = 14;           // protein codes listed in a tooltip
    const CLICK_SLOP = 6;           // px the pointer may move between press and release of a click
    const STORAGE_KEY = selectionMode ? 'lpmocat.selectionOverview' : 'lpmocat.overview';
    const REST_COLOR = '#d9d9d9';   // the rest of the catalogue (selection overview)
    const REST_DOT = '#c4c4c4';
    const REST_NAME = 'Rest of the catalogue';
    const PLOT_CONFIG = { displayModeBar: false, responsive: true, scrollZoom: false, doubleClick: false };

    // ---- helpers ----
    const plural = (n, word = 'protein') => `${n} ${word}${n === 1 ? '' : 's'}`;
    const signed = (x, decimals) => x.toFixed(decimals).replace('-', '−');
    const codesOf = list => list.map(p => p.code);
    const inOrder = (order, list) => order.filter(f => list.some(p => p.family === f));
    const familiesIn = list => inOrder(families, list);
    const ofFamily = (list, family) => list.filter(p => p.family === family);
    const tidy = x => Math.round(x * 1e6) / 1e6;

    function median(values) {
        const sorted = [...values].sort((a, b) => a - b);
        const mid = sorted.length >> 1;
        return sorted.length % 2 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2;
    }

    function el(tag, className, text) {
        const node = document.createElement(tag);
        if (className) node.className = className;
        if (text !== undefined) node.textContent = text;
        return node;
    }

    function loadState() {
        try {
            return JSON.parse(localStorage.getItem(STORAGE_KEY)) || {};
        } catch (e) {
            return {};
        }
    }

    function saveState(patch) {
        try {
            localStorage.setItem(STORAGE_KEY, JSON.stringify(Object.assign(loadState(), patch)));
        } catch (e) { /* no storage: the page opens in its default state */ }
    }

    // ---- tabs: Catalogue overview | Selection criteria ----
    const tabs = [...document.querySelectorAll('.left-tabs [role="tab"]')];
    const panes = { overview: root, selection: document.getElementById('selectionPane') };

    function showPane(name, { focus = false, initial = false } = {}) {
        tabs.forEach(tab => {
            const active = tab.dataset.pane === name;
            tab.setAttribute('aria-selected', String(active));
            tab.tabIndex = active ? 0 : -1;
            if (active && focus) tab.focus();
        });
        Object.entries(panes).forEach(([pane, node]) => {
            if (node) node.hidden = pane !== name;
        });
        // both tabs have the same width (main_page.css), so the tree stays as it is
        if (name === 'overview') renderCharts();
        if (!initial) saveState({ tab: name });
    }

    tabs.forEach((tab, i) => {
        tab.addEventListener('click', () => showPane(tab.dataset.pane));
        tab.addEventListener('keydown', event => {
            const step = { ArrowRight: 1, ArrowLeft: -1 }[event.key];
            if (!step) return;
            event.preventDefault();
            const next = tabs[(i + step + tabs.length) % tabs.length];
            showPane(next.dataset.pane, { focus: true });
        });
    });

    // Selection-change hook (updateDendrogramHighlighting in index.html calls it); set below
    // once the charts exist
    let renderCharts = () => {};
    window.catalogueOverview = { selectionChanged() {}, showPane };

    // Browse page: the data are in the page. Compare page: fetched (the whole catalogue, the
    // compared codes and their tree), while the card can already be folded.
    const selectionCard = selectionMode ? setupSelectionCard() : null;
    const dataScript = document.getElementById('catalogue-overview-data');
    let payload = null;
    if (typeof Plotly !== 'undefined') {
        try {
            if (selectionMode) {
                const response = await fetch(root.dataset.src);
                if (response.ok) payload = await response.json();
            } else if (dataScript) {
                payload = { overview: JSON.parse(dataScript.textContent) };
            }
        } catch (e) {
            payload = null;
        }
    }
    if (!payload || !payload.overview) {
        // no charts (e.g. the Plotly CDN is unreachable): say so, and open the filters
        root.replaceChildren(el('p', 'ov-hint', 'The overview charts could not be loaded.'));
        if (!selectionMode) showPane('selection', { initial: true });
        return;
    }

    const data = payload.overview;
    const families = data.families.map(f => f.name);
    const familyColor = Object.fromEntries(data.families.map(f => [f.name, f.color]));
    const stackOrder = data.stack_order;
    // The proteins the charts count: all of them on the browse page, the compared ones on the
    // compare page, where `rest` (the rest of the catalogue, while the card's switch shows it)
    // is drawn in grey behind them
    const compared = new Set(selectionMode ? payload.selected : data.proteins.map(p => p.code));
    const proteins = data.proteins.filter(p => compared.has(p.code));
    const notCompared = data.proteins.filter(p => !compared.has(p.code));
    let showRest = selectionMode && loadState().rest !== 'hidden';
    let rest = showRest ? notCompared : [];

    // Thin bars: the gap fraction that keeps n bars in a plot at most MAX_BAR px thick.
    const barGap = (n, height) => Math.max(0.3, 1 - MAX_BAR / ((height - AXIS_BAND) / Math.max(n, 1)));

    const AXIS = {
        gridcolor: GRID, linecolor: AXIS_LINE, showline: true, zeroline: false,
        fixedrange: true, ticks: '', automargin: true, tickfont: FONT,
    };

    function chartLayout({ xaxis = {}, yaxis = {}, ...rest }) {
        const axis = AXIS;
        return {
            margin: { l: 8, r: 28, t: 6, b: 6 },
            paper_bgcolor: SURFACE, plot_bgcolor: SURFACE,
            font: FONT, showlegend: false, dragmode: false,
            hoverdistance: 12, barcornerradius: 4,
            ...rest,
            xaxis: { ...axis, ...xaxis },
            yaxis: { ...axis, ...yaxis },
        };
    }

    const axisTitle = text => ({ text, font: FONT, standoff: 6 });
    // A bar's tooltip shows on the bar (a few px around it for thin ones), not in its whole row
    const BAR_HOVER = { hovermode: 'closest', hoverdistance: 4 };

    // Total at the end of each horizontal bar (y = category label, or index on a two-level axis);
    // with the rest of the catalogue (restTotals): "3 of 79" after its grey part
    const totalLabels = (ys, totals, restTotals = null) => ys.map((y, i) => {
        const all = totals[i] + (restTotals ? restTotals[i] : 0);
        return {
            x: all, y, text: restTotals ? `${totals[i]} of ${all}` : String(totals[i]), showarrow: false,
            xanchor: 'left', xshift: 4, font: { ...FONT, color: INK },
        };
    });
    // x range of horizontal bars ending at `ends`, with room for their total labels
    const barRange = ends => [0, Math.max(...ends, 1) * (rest.length ? 1.3 : 1.14)];

    // The grey part of a horizontal bar chart (the rest of the catalogue in each bar) and its
    // tooltip row
    const restBar = (y, counts) => ({
        type: 'bar', orientation: 'h', name: REST_NAME, y, x: counts,
        marker: { color: REST_COLOR, line: { color: SURFACE, width: 2 } }, hoverinfo: 'none',
    });
    const restRows = list => (rest.length ? [{ color: REST_COLOR, label: REST_NAME, value: list.length }] : []);
    const restCounts = groups => groups.map(g => (g.rest || []).length);

    // ---- selection (the protein checkboxes of the Selection criteria tab) ----
    const checkboxes = new Map(
        [...document.querySelectorAll('input[name="filter_protein_ids"]')].map(cb => [cb.value, cb]));
    const isSelected = code => !!(checkboxes.get(code) && checkboxes.get(code).checked);

    // Select all of `codes`, or deselect them when all are selected already; the page's change
    // handler re-reads every checkbox, so one event is enough (as for a tree branch).
    function toggleSelection(codes) {
        const boxes = codes.map(code => checkboxes.get(code)).filter(Boolean);
        const select = !boxes.every(cb => cb.checked);
        const changed = boxes.filter(cb => cb.checked !== select);
        changed.forEach(cb => { cb.checked = select; });
        if (changed.length) changed[0].dispatchEvent(new Event('change', { bubbles: true }));
    }

    function selectionNote(codes) {
        if (selectionMode) return '';   // the compare page selects nothing
        const n = codes.filter(isSelected).length;
        if (codes.length === 1) return n ? 'Selected · click to deselect' : 'Click to select';
        const selected = n ? `${n} of ${codes.length} selected · ` : '';
        return selected + (n === codes.length ? 'click to deselect all' : 'click to select all');
    }

    // ---- dendrogram link (compare page: the tree of the compared proteins) ----
    let selectionTree = null;
    function highlightInTree(codes) {
        const tree = selectionMode ? selectionTree : window.currentDendrogram;
        if (tree && typeof tree.setOverviewHighlight === 'function') {
            tree.setOverviewHighlight(codes && codes.length ? codes : null);
        }
    }

    // ---- tooltip (one for all charts; the table views carry the same values) ----
    const tooltip = el('div', 'ov-tooltip');
    tooltip.hidden = true;
    tooltip.setAttribute('aria-hidden', 'true');
    document.body.appendChild(tooltip);
    let shown = null;   // {event, content} of the visible tooltip

    // content: {heading, value, detail, rows: [{color, label, value}], codes, members, restCodes}
    // (members = the proteins it stands for, when not listed as codes; restCodes = those of the
    // rest of the catalogue, compare page)
    const codeList = codes => {
        const more = codes.length - MAX_CODES;
        const listed = codes.slice(0, MAX_CODES).join(', ');
        return more > 0 ? `${listed} … and ${more} more` : listed;
    };
    function showTooltip(event, content) {
        shown = { event, content };
        tooltip.replaceChildren();
        if (content.heading) tooltip.appendChild(el('div', 'ov-tip-heading', content.heading));
        tooltip.appendChild(el('div', 'ov-tip-value', content.value));
        if (content.detail) tooltip.appendChild(el('div', 'ov-tip-heading', content.detail));
        (content.rows || []).forEach(row => {
            const line = el('div', 'ov-tip-row');
            const key = el('span', 'ov-tip-key');
            key.style.background = row.color;
            line.append(key, el('span', 'ov-tip-row-label', row.label));
            if (row.value !== undefined) line.appendChild(el('span', 'ov-tip-row-value', String(row.value)));
            tooltip.appendChild(line);
        });
        const codes = content.codes || [];
        if (codes.length) tooltip.appendChild(el('div', 'ov-tip-codes', codeList(codes)));
        const restCodes = content.restCodes || [];
        if (restCodes.length) {
            tooltip.appendChild(el('div', 'ov-tip-codes ov-tip-rest', `Not compared: ${codeList(restCodes)}`));
        }
        const note = selectionNote(membersOf(content));
        if (note) tooltip.appendChild(el('div', 'ov-tip-note', note));
        tooltip.hidden = false;
        moveTooltip(event);
    }

    const membersOf = content => content.members || content.codes || [];

    function moveTooltip(event) {
        if (tooltip.hidden || !event || event.clientX === undefined) return;
        shown.event = event;
        const pad = 14;
        const { offsetWidth: w, offsetHeight: h } = tooltip;
        let x = event.clientX + pad;
        let y = event.clientY + pad;
        if (x + w > window.innerWidth - 8) x = event.clientX - pad - w;
        if (y + h > window.innerHeight - 8) y = event.clientY - pad - h;
        tooltip.style.left = `${Math.max(8, x)}px`;
        tooltip.style.top = `${Math.max(8, y)}px`;
    }

    function hideTooltip() {
        tooltip.hidden = true;
        shown = null;
    }

    // The hovered bars (all their stacked segments) stay, the other bars of their panel fade.
    // focus = {index | indices, curves, curve}: the bar index or a Set of them, the traces of the
    // panel (all when absent), and the one segment kept when only it is hovered (curve = its
    // trace); null restores all.
    function fadeOtherBars(plot, focus) {
        const kept = focus && (focus.indices || new Set([focus.index]));
        plot.querySelectorAll('.barlayer .trace').forEach(traceNode => {
            const calc = traceNode.__data__;
            const curve = calc && calc[0] && calc[0].trace ? calc[0].trace.index : null;
            const inPanel = !!focus && (!focus.curves || curve === null || focus.curves.has(curve));
            traceNode.querySelectorAll('.point').forEach((point, i) => {
                const d = point.__data__;
                const pointIndex = d && typeof d.i === 'number' ? d.i : i;
                const faded = inPanel
                    && (!kept.has(pointIndex) || (focus.curve !== undefined && curve !== focus.curve));
                point.style.opacity = faded ? '0.35' : '';
            });
        });
    }

    // ---- charts ----
    // Each chart: figure(view) -> {traces, layout}; hover(points, view) -> tooltip content or
    // null; barFocus(points, view) -> {index, curves?, curve?} of the hovered bar or null (see
    // fadeOtherBars); labels(plot, card) -> [{node, content, focus}] for the category labels of a
    // bar chart (and legend items in the card's caption); table() -> {head, numeric, rows, total}.

    // Inner category labels of a horizontal bar chart: label i stands for bar i.
    function barLabels(plot, hover) {
        return [...plot.querySelectorAll('.ytick text')]
            .filter(node => node.__data__ && Number.isInteger(node.__data__.x) && node.textContent.trim())
            .map(node => {
                const index = node.__data__.x;
                return { node, content: () => hover([{ pointIndex: index }]), focus: { index } };
            });
    }

    // The solid part of a family's bar (experimental structures) and its light part (predicted
    // models) hover and select on their own; the family label stands for the whole family, and
    // the subtitle's two swatches for all experimental structures or all predicted models.
    function familiesChart(height) {
        const shown = families.filter(f => proteins.some(p => p.family === f) || rest.some(p => p.family === f));
        const groups = shown.map(f => ({ label: f, proteins: ofFamily(proteins, f), rest: ofFamily(rest, f) }));
        const experimental = groups.map(g => g.proteins.filter(p => p.source === 'PDB'));
        const predicted = groups.map(g => g.proteins.filter(p => p.source !== 'PDB'));
        const totals = groups.map(g => g.proteins.length);
        const restTotals = restCounts(groups);
        const bar = (name, lists, colors) => ({
            type: 'bar', orientation: 'h', name, y: shown, x: lists.map(l => l.length),
            marker: { color: colors, line: { color: SURFACE, width: 2 } }, hoverinfo: 'none',
        });
        const PARTS = [   // trace order (the grey rest of the catalogue, if shown, comes last)
            { key: 'experimental', lists: experimental, alpha: '', label: 'Experimental structures (PDB)',
              name: 'experimental structures' },
            { key: 'predicted', lists: predicted, alpha: PREDICTED_ALPHA, label: 'Predicted models',
              name: 'predicted models' },
        ];
        const modelRows = (list, color) => MODEL_SOURCES.map(([source, name]) => ({
            color, label: `${name} model`, value: list.filter(p => p.source === source).length,
        }));
        function familyContent(i) {
            const f = shown[i];
            return {
                heading: `${f} family`,
                value: plural(totals[i]),
                detail: rest.length ? `of ${totals[i] + restTotals[i]} in the catalogue` : '',
                rows: [
                    { color: familyColor[f], label: 'Experimental structure (PDB)', value: experimental[i].length },
                    ...modelRows(predicted[i], familyColor[f] + PREDICTED_ALPHA),
                    ...restRows(groups[i].rest),
                ].filter(row => row.value),
                codes: codesOf(groups[i].proteins),
                restCodes: codesOf(groups[i].rest),
            };
        }
        function partContent(part, i) {
            const f = shown[i];
            const list = part.lists[i];
            if (!list.length) return null;
            return {
                heading: `${f} family: ${part.name}`,
                value: plural(list.length),
                detail: `of ${plural(totals[i])} in the family`,
                rows: part.key === 'predicted' ? modelRows(list, familyColor[f] + PREDICTED_ALPHA).filter(row => row.value) : [],
                codes: codesOf(list),
            };
        }
        // all families' solid (or light) parts
        function legendContent(part) {
            const list = part.lists.flat();
            return {
                heading: part.label,
                value: plural(list.length),
                detail: part.key === 'predicted'
                    ? modelRows(list).filter(row => row.value).map(row => plural(row.value, row.label)).join(' · ')
                    : '',
                rows: shown.map((f, i) => ({
                    color: familyColor[f] + part.alpha, label: f, value: part.lists[i].length,
                })).filter(row => row.value),
                codes: codesOf(list),
            };
        }
        return {
            figure: () => ({
                traces: PARTS.map(part => bar(part.label, part.lists, shown.map(f => familyColor[f] + part.alpha)))
                    .concat(rest.length ? [restBar(shown, restTotals)] : []),
                layout: chartLayout({
                    barmode: 'stack', bargap: barGap(shown.length, height), ...BAR_HOVER,
                    annotations: totalLabels(shown, totals, rest.length ? restTotals : null),
                    xaxis: { title: axisTitle('Proteins'), range: barRange(totals.map((t, i) => t + restTotals[i])) },
                    yaxis: { type: 'category', autorange: 'reversed', showgrid: false },
                }),
            }),
            hover(points) {
                const part = PARTS[points[0].curveNumber];
                return part ? partContent(part, points[0].pointIndex) : familyContent(points[0].pointIndex);
            },
            barFocus: points => ({ index: points[0].pointIndex, curve: points[0].curveNumber }),
            labels(plot, card) {
                const labels = barLabels(plot, points => familyContent(points[0].pointIndex));
                const all = new Set(shown.map((_, i) => i));
                card.querySelectorAll('.ov-legend-item[data-part]').forEach(node => {
                    const curve = PARTS.findIndex(part => part.key === node.dataset.part);
                    if (curve < 0) return;
                    labels.push({ node, content: () => legendContent(PARTS[curve]), focus: { indices: all, curve } });
                });
                return labels;
            },
            table: () => ({
                head: ['Family', 'Experimental', 'Predicted', 'Total', ...(rest.length ? ['Catalogue'] : [])],
                numeric: [false, true, true, true, true],
                rows: shown.map((f, i) => [f, experimental[i].length, predicted[i].length, totals[i],
                                           ...(rest.length ? [totals[i] + restTotals[i]] : [])]),
                total: ['All', experimental.flat().length, predicted.flat().length, proteins.length,
                        ...(rest.length ? [proteins.length + rest.length] : [])],
            }),
        };
    }

    // Bars of a property that papers may report or that may only follow from EC numbers
    // (regioselectivity, substrates), stacked by family. As in the families chart, a family's
    // solid part is what papers report and its light part the rest. show = 'all' | 'papers':
    // with 'papers' only the reported values count, and one last bar holds the proteins with none.
    // groups: [{label, proteins}]; sourceOf(p, label) = 'papers' | 'ec' | 'none'.
    const NOT_REPORTED = 'Not reported in papers';
    function evidenceChart({ groups: allGroups, sourceOf, show, height, groupTitle, derived, byCount = true }) {
        let groups = allGroups;
        if (show === 'papers') {
            const reportedIn = g => p => sourceOf(p, g.label) === 'papers';
            groups = allGroups
                .map(g => ({ label: g.label, proteins: g.proteins.filter(reportedIn(g)), rest: (g.rest || []).filter(reportedIn(g)) }))
                .filter(g => g.proteins.length || g.rest.length);
            // largest first, unless the groups have their own order (C1, C4, C1/C4)
            if (byCount) groups.sort((a, b) => (b.proteins.length - a.proteins.length) || (b.rest.length - a.rest.length));
            const reported = new Set(groups.flatMap(g => [...g.proteins, ...g.rest]));
            const missing = proteins.filter(p => !reported.has(p));
            const restMissing = rest.filter(p => !reported.has(p));
            if (missing.length || restMissing.length) {
                groups.push({ label: NOT_REPORTED, proteins: missing, rest: restMissing, missing: true });
            }
        }
        const labels = groups.map(g => g.label);
        const totals = groups.map(g => g.proteins.length);
        const restTotals = restCounts(groups);
        const sourceIn = (g, p) => (g.missing ? 'ec' : sourceOf(p, g.label));
        const sources = [
            { key: 'papers', label: 'reported in papers', alpha: '' },
            { key: 'ec', label: derived, alpha: PREDICTED_ALPHA },
            { key: 'none', label: 'not known', alpha: PREDICTED_ALPHA },
        ];
        const count = (g, f, key) => g.proteins.filter(p => p.family === f && sourceIn(g, p) === key).length;
        const traces = inOrder(stackOrder, proteins).flatMap(f => sources.map(s => ({
            type: 'bar', orientation: 'h', name: `${f} ${s.key}`, y: labels,
            x: groups.map(g => count(g, f, s.key)),
            marker: { color: familyColor[f] + s.alpha, line: { color: SURFACE, width: 2 } },
            hoverinfo: 'none',
        }))).filter(trace => trace.x.some(Boolean))
            .concat(rest.length ? [restBar(labels, restTotals)] : []);
        return {
            figure: () => ({
                traces,
                layout: chartLayout({
                    barmode: 'stack', bargap: barGap(labels.length, height), ...BAR_HOVER,
                    annotations: totalLabels(labels, totals, rest.length ? restTotals : null),
                    xaxis: { title: axisTitle('Proteins'), range: barRange(totals.map((t, i) => t + restTotals[i])) },
                    yaxis: { type: 'category', autorange: 'reversed', showgrid: false, dtick: 1 },
                }),
            }),
            hover(points) {
                const g = groups[points[0].pointIndex];
                const bySource = sources.map(s => [s, g.proteins.filter(p => sourceIn(g, p) === s.key).length])
                    .filter(([, n]) => n);
                return {
                    heading: g.missing ? NOT_REPORTED : g.label,
                    value: plural(g.proteins.length),
                    detail: g.missing ? `No paper in the sources gives a ${groupTitle.toLowerCase()} for these`
                        : bySource.map(([s, n]) => `${n} ${s.label}`).join(' · '),
                    rows: familiesIn(g.proteins).flatMap(f => sources.map(s => ({
                        color: familyColor[f] + s.alpha,
                        label: g.missing ? f : `${f}, ${s.label}`,
                        value: count(g, f, s.key),
                    }))).filter(row => row.value).concat(restRows(g.rest || [])),
                    codes: codesOf(g.proteins),
                    restCodes: codesOf(g.rest || []),
                };
            },
            barFocus: points => ({ index: points[0].pointIndex }),
            labels(plot) {
                return barLabels(plot, this.hover);
            },
            table: () => ({
                head: [groupTitle, 'Reported in papers', derived[0].toUpperCase() + derived.slice(1), 'Total', 'By family',
                       ...(rest.length ? ['Catalogue'] : [])],
                numeric: [false, true, true, true, false, true],
                rows: groups.map((g, i) => [
                    g.missing ? NOT_REPORTED : g.label,
                    g.proteins.filter(p => sourceIn(g, p) === 'papers').length,
                    g.proteins.filter(p => sourceIn(g, p) !== 'papers').length,
                    g.proteins.length,
                    familiesIn(g.proteins).map(f => `${f} ${ofFamily(g.proteins, f).length}`).join(' · '),
                    ...(rest.length ? [g.proteins.length + restTotals[i]] : []),
                ]),
            }),
        };
    }

    // Regioselectivity (overview.regioselectivity_of: p.regio = [C1 | C4 | C1/C4 | Unknown,
    // papers | ec | none]) and substrates (p.substrate_papers: those a paper reports)
    const REGIO_ORDER = ['C1', 'C4', 'C1/C4', 'Unknown'];
    const SHOW_CONTROL = {
        key: 'show', label: 'Values shown', type: 'buttons',
        options: [{ value: 'all', label: 'All' }, { value: 'papers', label: 'Reported in papers' }],
    };
    const legendHtml = derived => '<span class="ov-swatch" style="background:#666"></span> reported in papers '
        + `&nbsp;<span class="ov-swatch" style="background:#c9c9c9"></span> ${derived}`;
    function regioselectivity({ show }) {
        const counted = [...proteins, ...rest];
        const labels = REGIO_ORDER.filter(r => counted.some(p => p.regio[0] === r))
            .concat([...new Set(counted.map(p => p.regio[0]))].filter(r => !REGIO_ORDER.includes(r)));
        const derived = 'from the EC numbers';
        return {
            title: 'Regioselectivity',
            sub: show === 'papers' ? 'As papers report it' : `Reported in papers or ${derived}`,
            subHtml: show === 'papers' ? 'As papers report it; the last bar has no report'
                : `C1, C4 or both: ${legendHtml(derived)}`,
            chart: evidenceChart({
                groups: labels.map(r => ({
                    label: r, proteins: proteins.filter(p => p.regio[0] === r), rest: rest.filter(p => p.regio[0] === r),
                })),
                sourceOf: p => p.regio[1], show, height: PLOT_HEIGHTS.normal,
                groupTitle: 'Regioselectivity', derived, byCount: false,
            }),
        };
    }
    function substrates({ show }) {
        const derived = 'from the EC numbers or the name';
        return {
            title: 'Substrates',
            sub: show === 'papers' ? 'Substrates a paper reports activity on; a protein can count twice'
                : 'Reported in papers, or from the EC numbers or the protein name; a protein can count twice',
            subHtml: show === 'papers'
                ? 'Substrates a paper reports activity on (a protein can count twice); the last bar has none'
                : `${legendHtml(derived)}; a protein can count twice`,
            chart: evidenceChart({
                groups: groupsBy(p => (p.substrates.length ? p.substrates : ['None assigned']), ['None assigned']),
                sourceOf: (p, label) => ((p.substrate_papers || []).includes(label) ? 'papers' : 'ec'),
                show, height: PLOT_HEIGHTS.normal, groupTitle: 'Substrate', derived,
            }),
        };
    }

    // Horizontal bars, one per group ({label, proteins, rest?, levels?, detail?}), stacked by
    // family, with the group's rest of the catalogue in grey after them. With `levels` ([outer,
    // inner] per group) the y axis has two levels of labels; `column` ({head, value: group ->
    // text}) adds a column to the table view.
    function groupChart(groups, groupTitle, height, column = null) {
        const twoLevel = groups.some(g => g.levels);
        const labels = groups.map(g => g.label);
        const yValues = twoLevel ? [groups.map(g => g.levels[0]), groups.map(g => g.levels[1])] : labels;
        const stacked = inOrder(stackOrder, groups.flatMap(g => g.proteins));
        const totals = groups.map(g => g.proteins.length);
        const restTotals = restCounts(groups);
        const breakdown = i => familiesIn(groups[i].proteins)
            .map(f => ({ f, n: ofFamily(groups[i].proteins, f).length }));
        // several bars together, e.g. all phyla of the fungi
        function groupContent(name, indices) {
            const members = [...new Set(indices.flatMap(i => groups[i].proteins))];
            const others = [...new Set(indices.flatMap(i => groups[i].rest || []))];
            const subgroups = indices.map(i => groups[i].levels[1].trim()).filter(Boolean);
            const species = new Set(members.map(p => p.species)).size;
            return {
                heading: name,
                value: plural(members.length),
                detail: [subgroups.length > 1 ? subgroups.join(', ') : '', `${species} species`]
                    .filter(Boolean).join(' · '),
                rows: familiesIn(members).map(f => ({ color: familyColor[f], label: f, value: ofFamily(members, f).length }))
                    .concat(restRows(others)),
                codes: codesOf(members),
                restCodes: codesOf(others),
            };
        }
        return {
            figure: () => ({
                traces: stacked.map(f => ({
                    type: 'bar', orientation: 'h', name: f, y: yValues,
                    x: groups.map(g => ofFamily(g.proteins, f).length),
                    marker: { color: familyColor[f], line: { color: SURFACE, width: 2 } }, hoverinfo: 'none',
                })).concat(rest.length ? [restBar(yValues, restTotals)] : []),
                layout: chartLayout({
                    barmode: 'stack', bargap: barGap(groups.length, height), ...BAR_HOVER,
                    annotations: totalLabels(twoLevel ? groups.map((_, i) => i) : labels, totals,
                                             rest.length ? restTotals : null),
                    xaxis: { title: axisTitle('Proteins'), range: barRange(totals.map((t, i) => t + restTotals[i])) },
                    yaxis: {
                        type: twoLevel ? 'multicategory' : 'category', autorange: 'reversed',
                        showgrid: false, dividercolor: AXIS_LINE,
                        dtick: 1,   // every bar named (Plotly thins crowded category ticks)
                    },
                }),
            }),
            hover(points) {
                const i = points[0].pointIndex;
                return {
                    heading: groups[i].label,
                    value: plural(totals[i]),
                    detail: groups[i].detail,
                    rows: breakdown(i).map(c => ({ color: familyColor[c.f], label: c.f, value: c.n }))
                        .concat(restRows(groups[i].rest || [])),
                    codes: codesOf(groups[i].proteins),
                    restCodes: codesOf(groups[i].rest || []),
                };
            },
            barFocus: points => ({ index: points[0].pointIndex }),
            labels(plot) {
                const labels = barLabels(plot, this.hover);
                // outer labels of a two-level axis: the whole group ("Fungi")
                plot.querySelectorAll('.ytick2 text').forEach(node => {
                    const name = node.textContent.trim();
                    const indices = groups.map((g, i) => (g.levels && g.levels[0] === name ? i : -1)).filter(i => i >= 0);
                    if (!indices.length) return;
                    labels.push({ node, content: () => groupContent(name, indices), focus: { indices: new Set(indices) } });
                });
                return labels;
            },
            table: () => ({
                head: [groupTitle, 'Proteins', 'By family', ...(column ? [column.head] : []),
                       ...(rest.length ? ['Catalogue'] : [])],
                numeric: [false, true, false, ...(column ? [false] : []), true],
                rows: groups.map((g, i) => [g.label, totals[i], breakdown(i).map(c => `${c.f} ${c.n}`).join(' · '),
                                            ...(column ? [column.value(g)] : []),
                                            ...(rest.length ? [totals[i] + restTotals[i]] : [])]),
            }),
        };
    }

    // groups ({label, proteins, rest}) from a list-valued property of the proteins and of the rest
    // of the catalogue; `last` labels go to the end, the others largest first (by their proteins,
    // then by their rest of the catalogue)
    function groupsBy(valuesOf, last = []) {
        const map = new Map();
        const add = key => p => valuesOf(p).forEach(v => {
            if (!map.has(v)) map.set(v, { proteins: [], rest: [] });
            map.get(v)[key].push(p);
        });
        proteins.forEach(add('proteins'));
        rest.forEach(add('rest'));
        const rank = label => last.indexOf(label);
        return [...map].map(([label, g]) => ({ label, ...g }))
            .sort((a, b) => (rank(a.label) - rank(b.label))
                || (b.proteins.length - a.proteins.length) || (b.rest.length - a.rest.length)
                || a.label.localeCompare(b.label));
    }

    // Taxon groups (overview.taxon_group) in the server's group order, largest subgroup first.
    function taxonomyGroups() {
        const map = new Map();
        const add = key => p => {
            const [group, sub] = p.taxon;
            const id = `${group}\t${sub || ''}`;
            if (!map.has(id)) map.set(id, { group, sub, proteins: [], rest: [] });
            map.get(id)[key].push(p);
        };
        proteins.forEach(add('proteins'));
        rest.forEach(add('rest'));
        const order = group => {
            const i = data.taxon_groups.indexOf(group);
            return i < 0 ? data.taxon_groups.length : i;
        };
        const groupSize = (list, group) => list.filter(p => p.taxon[0] === group).length;
        return [...map.values()]
            .sort((a, b) => (order(a.group) - order(b.group))
                || (groupSize(proteins, b.group) - groupSize(proteins, a.group))
                || (groupSize(rest, b.group) - groupSize(rest, a.group))
                || a.group.localeCompare(b.group) || (b.proteins.length - a.proteins.length)
                || (b.rest.length - a.rest.length))
            .map(g => {
                const members = g.proteins.length ? g.proteins : g.rest;
                const species = new Set(members.map(p => p.species)).size;
                // the lineage down to the phylum, when the group's proteins share it
                const lineages = new Set(members.map(p => p.lineage.split(' › ').slice(0, 3).join(' › ')));
                const lineage = lineages.size === 1 ? [...lineages][0] : '';
                return {
                    label: g.sub ? `${g.group} · ${g.sub}` : g.group,
                    levels: [g.group, g.sub || ' '],
                    detail: [lineage, `${species} species`].filter(Boolean).join(' · '),
                    proteins: g.proteins,
                    rest: g.rest,
                };
            });
    }

    // Bin width for a histogram of `values` with fewer than maxBins bins
    const BIN_WIDTHS = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100];
    function niceBinWidth(values, maxBins = 20) {
        const span = Math.max(...values) - Math.min(...values);
        return BIN_WIDTHS.find(w => span / w < maxBins) || 200;
    }

    // Numeric properties side by side, one panel each (a single panel for a single property):
    // histograms stacked by family ('stacked') or one dot per protein in a row per family with
    // the family median ('family'). The panels share the family rows / count axis. The rest of the
    // catalogue (compare page) stacks in grey on the bins, or lies behind as grey dots with a grey
    // tick at the family's median in the whole catalogue.
    // panels: [{name, value: protein -> number or undefined, binWidth? (or maxBins?),
    //           detail?: protein -> text, zeroline?}]
    function distributionChart({ panels, format, unit, title, height }) {
        const multi = panels.length > 1;
        const label = v => `${format(v)}${unit ? ' ' + unit : ''}`;
        const measuredIn = panel => p => typeof panel.value(p) === 'number';
        const measuredAny = proteins.filter(p => panels.some(panel => measuredIn(panel)(p)));
        const fams = familiesIn(measuredAny.concat(rest.filter(p => panels.some(panel => measuredIn(panel)(p)))));
        const familySize = f => ofFamily(proteins, f).length;
        const medianOf = (list, panel) => (list.length ? median(list.map(panel.value)) : null);

        const prepared = panels.map(panel => {
            const measured = proteins.filter(measuredIn(panel));
            const restMeasured = rest.filter(measuredIn(panel));
            const values = measured.concat(restMeasured).map(panel.value);
            const binWidth = panel.binWidth || (values.length ? niceBinWidth(values, panel.maxBins) : 1);
            const lo = values.length ? tidy(Math.floor(Math.min(...values) / binWidth) * binWidth) : 0;
            const nBins = values.length ? Math.floor((Math.max(...values) - lo) / binWidth + 1e-9) + 1 : 0;
            const bins = Array.from({ length: nBins }, (_, i) => ({ start: tidy(lo + i * binWidth), proteins: [], rest: [] }));
            const binOf = p => bins[Math.min(nBins - 1, Math.floor((panel.value(p) - lo) / binWidth + 1e-9))];
            measured.forEach(p => binOf(p).proteins.push(p));
            restMeasured.forEach(p => binOf(p).rest.push(p));
            const perFamily = fams.map(f => ofFamily(measured, f));
            const restPerFamily = fams.map(f => ofFamily(restMeasured, f));
            const medians = perFamily.map(list => medianOf(list, panel));
            // the family's median in the whole catalogue (compare page)
            const catalogueMedians = perFamily.map((list, k) => medianOf(list.concat(restPerFamily[k]), panel));
            const integer = values.every(Number.isInteger) && Number.isInteger(binWidth);
            return { ...panel, binWidth, integer, measured, restMeasured, bins, perFamily, restPerFamily, medians,
                     catalogueMedians };
        });

        // integer values: the bin's first and last value (210–219); else its limits (−0.30 to −0.25)
        const binLabel = (panel, b) => {
            const last = panel.integer ? b.start + panel.binWidth - 1 : tidy(b.start + panel.binWidth);
            const text = panel.integer && panel.binWidth === 1 ? format(b.start)
                : `${format(b.start)}${panel.integer && b.start >= 0 ? '–' : ' to '}${format(last)}`;
            return unit ? `${text} ${unit}` : text;
        };
        const named = (panel, text) => (multi ? `${panel.name} · ${text}` : text);

        // axes of panel i: x/y, x2/y2, ... side by side; with shareY the y axes follow the first
        const gap = multi ? 0.035 : 0;
        const width = (1 - gap * (panels.length - 1)) / panels.length;
        const axisRef = (i, letter) => (i ? `${letter}${i + 1}` : letter);
        const axisKey = (i, letter) => `${letter}axis${i ? i + 1 : ''}`;
        function panelAxes(x, y, shareY) {
            const axes = {};
            prepared.forEach((panel, i) => {
                const left = i * (width + gap);
                axes[axisKey(i, 'x')] = {
                    ...AXIS, domain: [left, left + width], anchor: axisRef(i, 'y'),
                    ...(multi ? { nticks: 4 } : {}),
                    ...(panel.zeroline ? { zeroline: true, zerolinecolor: '#9a9a9a', zerolinewidth: 1 } : {}),
                    ...x(i),
                };
                axes[axisKey(i, 'y')] = {
                    ...AXIS, anchor: axisRef(i, 'x'),
                    ...(i && shareY ? { matches: 'y', showticklabels: false } : {}), ...y(i),
                };
            });
            return axes;
        }
        // panel names above the panels; the x title under the middle one
        const panelTitles = () => (multi ? prepared.map((panel, i) => ({
            text: panel.name, xref: `${axisRef(i, 'x')} domain`, yref: 'paper', x: 0.5, y: 1,
            yanchor: 'bottom', showarrow: false, font: { ...FONT, size: 12, color: INK },
        })) : []);
        const xTitle = i => ((!multi || i === Math.floor(panels.length / 2)) ? { title: axisTitle(title) } : {});

        // deterministic vertical spread of the dots, from the protein code
        const jitter = code => {
            let h = 7;
            for (const c of code) h = (h * 31 + c.charCodeAt(0)) >>> 0;
            return ((h % 1000) / 1000 - 0.5) * 0.46;
        };

        let traceInfo = [];   // per trace of the last figure: {panel, kind, family?, rows?}

        function stacked() {
            const traces = [];
            traceInfo = [];
            prepared.forEach((panel, i) => {
                // a bin of whole numbers holds start … start + width − 1: its bar is centred on
                // them (a charge of 0 on the zero line), a bin of real numbers on its interval
                const middle = panel.integer ? (panel.binWidth - 1) / 2 : panel.binWidth / 2;
                const centers = panel.bins.map(b => tidy(b.start + middle));
                inOrder(stackOrder, panel.measured).forEach(f => {
                    traces.push({
                        type: 'bar', name: f, x: centers, xaxis: axisRef(i, 'x'), yaxis: axisRef(i, 'y'),
                        y: panel.bins.map(b => ofFamily(b.proteins, f).length),
                        marker: { color: familyColor[f], line: { color: SURFACE, width: multi ? 1 : 1.5 } },
                        hoverinfo: 'none',
                    });
                    traceInfo.push({ panel: i, kind: 'bin' });
                });
                if (rest.length) {
                    traces.push({
                        type: 'bar', name: REST_NAME, x: centers, xaxis: axisRef(i, 'x'), yaxis: axisRef(i, 'y'),
                        y: panel.bins.map(b => b.rest.length),
                        marker: { color: REST_COLOR, line: { color: SURFACE, width: multi ? 1 : 1.5 } },
                        hoverinfo: 'none',
                    });
                    traceInfo.push({ panel: i, kind: 'bin' });
                }
            });
            return {
                traces,
                layout: chartLayout({
                    barmode: 'stack', bargap: 0.08, ...BAR_HOVER,
                    margin: { l: 8, r: 12, t: multi ? 20 : 6, b: 6 },
                    annotations: panelTitles(),
                    ...panelAxes(i => ({ showgrid: false, ...xTitle(i) }),
                                 i => ({ rangemode: 'tozero', nticks: 5, ...(i ? {} : { title: axisTitle('Proteins') }) }),
                                 false),
                }),
            };
        }

        function byFamily() {
            const traces = [];
            traceInfo = [];
            const size = multi ? 6 : 8;
            prepared.forEach((panel, i) => {
                const axes = { xaxis: axisRef(i, 'x'), yaxis: axisRef(i, 'y') };
                if (panel.restMeasured.length) {
                    const list = panel.restPerFamily.flat();
                    traces.push({
                        type: 'scatter', mode: 'markers', name: REST_NAME, ...axes,
                        x: list.map(panel.value), y: list.map(p => fams.indexOf(p.family) + jitter(p.code)),
                        marker: { color: REST_DOT, size: size - 2, line: { color: SURFACE, width: 0.8 } },
                        hoverinfo: 'none',
                    });
                    traceInfo.push({ panel: i, kind: 'restdot', list });
                }
                panel.perFamily.forEach((list, k) => {
                    if (!list.length) return;
                    const selected = list.map(p => isSelected(p.code));
                    traces.push({
                        type: 'scatter', mode: 'markers', name: fams[k], ...axes,
                        x: list.map(panel.value), y: list.map(p => k + jitter(p.code)),
                        marker: {
                            color: familyColor[fams[k]],
                            size: selected.map(on => (on ? size + 2 : size)),
                            line: {
                                color: selected.map(on => (on ? SELECTED_RING : SURFACE)),
                                width: selected.map(on => (on ? 2 : 1.2)),
                            },
                        },
                        hoverinfo: 'none',
                    });
                    traceInfo.push({ panel: i, kind: 'dot', family: k });
                });
                if (rest.length) {
                    const catalogueRows = fams.map((_, k) => k).filter(k => panel.catalogueMedians[k] !== null);
                    traces.push({
                        type: 'scatter', mode: 'markers', name: 'catalogue median', ...axes,
                        x: catalogueRows.map(k => panel.catalogueMedians[k]), y: catalogueRows,
                        marker: { symbol: 'line-ns', size: multi ? 13 : 15, line: { color: '#a6a6a6', width: 2.5 } },
                        hoverinfo: 'none',
                    });
                    traceInfo.push({ panel: i, kind: 'catalogueMedian', rows: catalogueRows });
                }
                const rows = fams.map((_, k) => k).filter(k => panel.perFamily[k].length);
                traces.push({
                    type: 'scatter', mode: 'markers', name: 'median', ...axes,
                    x: rows.map(k => panel.medians[k]), y: rows,
                    marker: { symbol: 'line-ns', size: multi ? 13 : 15, line: { color: INK, width: 2.5 } },
                    hoverinfo: 'none',
                });
                traceInfo.push({ panel: i, kind: 'median', rows });
            });
            return {
                traces,
                layout: chartLayout({
                    hovermode: 'closest',
                    margin: { l: 8, r: 12, t: multi ? 20 : 6, b: 6 },
                    annotations: panelTitles(),
                    ...panelAxes(i => xTitle(i), () => ({
                        tickvals: fams.map((_, k) => k), ticktext: fams, showgrid: false,
                        range: [fams.length - 0.5, -0.5],
                    }), true),
                }),
            };
        }

        return {
            hasViews: true,
            figure: view => (view === 'family' ? byFamily() : stacked()),
            hover(points) {
                const point = points[0];
                const info = traceInfo[point.curveNumber];
                if (!info) return null;
                const panel = prepared[info.panel];
                if (info.kind === 'bin') {
                    const bin = panel.bins[point.pointIndex];
                    if (!bin || !(bin.proteins.length || bin.rest.length)) return null;
                    return {
                        heading: named(panel, binLabel(panel, bin)),
                        value: plural(bin.proteins.length),
                        rows: familiesIn(bin.proteins).map(f => ({
                            color: familyColor[f], label: f, value: ofFamily(bin.proteins, f).length,
                        })).concat(restRows(bin.rest)),
                        codes: codesOf(bin.proteins),
                        restCodes: codesOf(bin.rest),
                    };
                }
                if (info.kind === 'restdot') {
                    const p = info.list[point.pointIndex];
                    return {
                        heading: `${p.code} · not compared`,
                        value: named(panel, label(panel.value(p))),
                        rows: [
                            { color: REST_DOT, label: `${p.family} · ${p.species}` },
                            { color: 'transparent', label: panel.detail ? panel.detail(p) : `Structure ${p.structure}` },
                        ],
                        members: [],
                    };
                }
                if (info.kind === 'catalogueMedian') {
                    const k = info.rows[point.pointIndex];
                    const all = panel.perFamily[k].concat(panel.restPerFamily[k]);
                    const vals = all.map(panel.value);
                    return {
                        heading: `${fams[k]}${multi ? ' ' + panel.name : ''} median of the catalogue`,
                        value: label(panel.catalogueMedians[k]),
                        rows: [{
                            color: REST_DOT,
                            label: `${all.length} proteins, ${format(Math.min(...vals))} to ${label(Math.max(...vals))}`,
                        }],
                        codes: codesOf(panel.perFamily[k]),
                    };
                }
                if (info.kind === 'median') {
                    const k = info.rows[point.pointIndex];
                    const list = panel.perFamily[k];
                    const vals = list.map(panel.value);
                    const of = list.length < familySize(fams[k]) ? ` of ${familySize(fams[k])}` : '';
                    return {
                        heading: `${fams[k]}${multi ? ' ' + panel.name : ''} median`,
                        value: label(panel.medians[k]),
                        rows: [{
                            color: familyColor[fams[k]],
                            label: `${list.length}${of} proteins, ${format(Math.min(...vals))} to ${label(Math.max(...vals))}`,
                        }],
                        codes: codesOf(list),
                    };
                }
                const p = panel.perFamily[info.family][point.pointIndex];
                return {
                    heading: p.code,
                    value: named(panel, label(panel.value(p))),
                    rows: [
                        { color: familyColor[p.family], label: `${p.family} · ${p.species}` },
                        { color: 'transparent', label: panel.detail ? panel.detail(p) : `Structure ${p.structure}` },
                    ],
                    members: [p.code],
                };
            },
            barFocus(points) {
                const info = traceInfo[points[0].curveNumber];
                if (!info || info.kind !== 'bin') return null;
                const curves = new Set(traceInfo.map((t, c) => (t.panel === info.panel ? c : -1)).filter(c => c >= 0));
                return { index: points[0].pointIndex, curves };
            },
            dependsOnSelection: view => view === 'family',
            table: () => {
                const stats = (list, panel) => {
                    const vals = list.map(panel.value);
                    return { n: vals.length, min: Math.min(...vals), med: median(vals), max: Math.max(...vals) };
                };
                if (!multi) {
                    const panel = prepared[0];
                    const row = (name, list, others) => {
                        const st = stats(list, panel);
                        const cells = st.n ? [name, st.n, format(st.min), format(st.med), format(st.max)]
                            : [name, 0, '–', '–', '–'];
                        if (!rest.length) return cells;
                        const all = list.concat(others);
                        return cells.concat(all.length ? format(medianOf(all, panel)) : '–');
                    };
                    return {
                        head: ['Family', 'Proteins', 'Min', 'Median', 'Max', ...(rest.length ? ['Catalogue median'] : [])],
                        numeric: [false, true, true, true, true, true],
                        rows: fams.map((f, k) => row(f, panel.perFamily[k], panel.restPerFamily[k])),
                        total: row('All', panel.measured, panel.restMeasured),
                    };
                }
                // one column per panel: median (min–max), and how many have it when not all do
                const cell = (list, panel, total) => {
                    if (!list.length) return '–';
                    const st = stats(list, panel);
                    const range = `${format(st.min)}${st.min < 0 ? ' to ' : '–'}${format(st.max)}`;
                    return `${format(st.med)} (${range})`
                        + (st.n < total ? ` · ${st.n} of ${total}` : '');
                };
                return {
                    head: ['Family', ...prepared.map(panel => panel.name)],
                    numeric: prepared.map(() => false).concat(false),
                    rows: fams.map((f, k) => [`${f} (${familySize(f)})`,
                        ...prepared.map(panel => cell(panel.perFamily[k], panel, familySize(f)))]),
                    total: [`All (${measuredAny.length})`,
                        ...prepared.map(panel => cell(panel.measured, panel, measuredAny.length))],
                };
            },
        };
    }

    // A card whose chart depends on settings: the scope of the Properties heading (catalytic
    // domain or all domains; not with usesScope: false) and the card's own controls, [{key, label,
    // type: 'buttons' | 'select', options: [{value, label}], when?: values -> shown}], drawn by
    // setupCard. build({scope, ...values}) returns {title, sub, chart}; the choices are remembered.
    function configurable(key, { controls, build, usesScope = true }) {
        const saved = (loadState().measures || {})[key];
        const values = {};
        controls.forEach(control => {
            const stored = saved && typeof saved === 'object' ? saved[control.key] : undefined;
            values[control.key] = control.options.some(o => o.value === stored) ? stored : control.options[0].value;
        });
        const built = new Map();
        const current = () => {
            const settings = usesScope ? { scope, ...values } : { ...values };
            const id = JSON.stringify(settings);
            if (!built.has(id)) built.set(id, build(settings));
            return built.get(id);
        };
        const wrapper = {
            controls,
            values,
            current,
            usesScope,
            hasViews: true,
            set(controlKey, value) {
                values[controlKey] = value;
                saveState({ measures: { ...(loadState().measures || {}), [key]: { ...values } } });
            },
            dependsOnSelection: v => {
                const chart = current().chart;   // bar charts (groupChart) do not
                return !!(chart.dependsOnSelection && chart.dependsOnSelection(v));
            },
        };
        ['figure', 'hover', 'barFocus', 'table'].forEach(name => {
            wrapper[name] = (...args) => current().chart[name](...args);
        });
        // category labels of a bar chart (groupChart); distributions have none
        wrapper.labels = (...args) => {
            const chart = current().chart;
            return chart.labels ? chart.labels(...args) : [];
        };
        return wrapper;
    }

    // ---- assemble ----
    const SCOPES = { cd: 'catalytic domain', all: 'all domains' };
    const charge = v => (v > 0 ? '+' : '') + signed(v, Number.isInteger(v) ? 0 : 1);
    const residues = v => (Number.isInteger(v) ? String(v) : v.toFixed(1));
    // Residue sets a distribution can describe. 'loops' and 'binding' live inside the catalytic
    // domain, so they carry one value each (measures.sets) instead of a catalytic/all pair.
    const RESIDUE_OPTIONS = [
        { value: 'all', label: 'All residues' },
        { value: 'surface', label: 'Surface' },
        { value: 'loops', label: 'Loop regions' },
        { value: 'binding', label: 'Substrate-binding' },
    ];
    // inside a loop, 'loops' would be the loop itself
    const LOOP_RESIDUE_OPTIONS = RESIDUE_OPTIONS.filter(o => o.value !== 'loops');
    const RESIDUE_SET = { loops: 'loops', binding: 'binding' };   // -> measures.sets
    const setOf = (p, which) => ((p.measures.sets || {})[RESIDUE_SET[which]] || {});
    const chargeValue = (p, which, sc) =>
        which === 'all' ? p.measures.charge[sc]
            : which === 'surface' ? p.measures.surface_charge[sc]
                : setOf(p, which).charge;
    const hydrophobicityValue = (p, which, scale, sc) =>
        which === 'all' ? p.measures.hydrophobicity[scale][sc]
            : which === 'surface' ? p.measures.surface_hydrophobicity[scale][sc]
                : (setOf(p, which).hydrophobicity || {})[scale];
    // card titles per residue set: "Net charge (catalytic domain)", "Surface hydrophobicity
    // (all domains)", "Charge of the loop regions", "Substrate-binding hydrophobicity"
    const setTitle = (which, what) => (
        which === 'surface' ? `Surface ${what}`
            : which === 'loops' ? `${capitalise(what)} of the loop regions`
                : which === 'binding' ? `Substrate-binding ${what}`
                    : capitalise(what));
    // the same for a panel of the loop card, which is already one loop
    const loopTitle = (which, what) => `Loop ${what}`
        + (which === 'surface' ? ' (surface residues)'
            : which === 'binding' ? ' (substrate-binding residues)' : '');
    const insideDomain = '; these residues all lie in the catalytic domain';
    const capitalise = text => text.charAt(0).toUpperCase() + text.slice(1);

    const SCALE_OPTIONS = data.hydrophobicity_scales.map((scale, i) => ({
        value: scale.name, label: i === 0 ? `GRAVY (${scale.name})` : scale.name,
    }));
    const isGravy = scale => scale === data.hydrophobicity_scales[0].name;
    // subtitle of a hydrophobicity chart: the scale's meaning, for the residues `of`
    function scaleSub(scale, of) {
        if (isGravy(scale)) return `GRAVY: mean Kyte–Doolittle hydropathy of ${of}; below 0 = hydrophilic`;
        const info = data.hydrophobicity_scales.find(x => x.name === scale);
        return `Mean ${scale} value of ${of}; higher = more ${info && info.hydrophilic ? 'hydrophilic' : 'hydrophobic'}`;
    }
    const surfaceOf = what => `the surface-exposed ${what} (relative SASA ≥ ${data.surface_rsa} %)`;
    // the residues a card describes, named for its subtitle
    const residuesOf = (which, what) => (
        which === 'surface' ? surfaceOf(what)
            : which === 'loops' ? `the ${what} of the loop regions taken together`
                : which === 'binding' ? `the substrate-binding ${what}: those a plane laid on the `
                    + 'copper face touches, the stand-in for a flat substrate'
                    : `the ${what}`);

    // One panel per loop region with a value in any protein (or in the rest of the catalogue). A
    // loop is [residues, 'start-end']; its other values are p.loop_values[loop]
    // (overview.loop_profile). valueOf(p, loop) gives a number or null; withLength adds the loop
    // length to a dot's tooltip.
    const loopPanels = (valueOf, withLength, extra = {}) => data.loop_names
        .map(loop => ({
            name: loop,
            maxBins: 12,   // narrow panels
            ...extra,
            value: p => {
                const v = p.loops[loop] ? valueOf(p, loop) : null;
                return typeof v === 'number' ? v : undefined;
            },
            detail: p => `${loop}: residues ${p.loops[loop][1]}`
                + (withLength ? ` (${plural(p.loops[loop][0], 'residue')})` : ''),
        }))
        .filter(panel => proteins.concat(rest).some(p => panel.value(p) !== undefined));
    const loopValues = (p, loop) => p.loop_values[loop] || {};

    // The loop card: one property of the loop regions (always in the catalytic domain)
    function loopProperty({ property, scale, residues: which }) {
        const setIndex = { all: 0, surface: 1, binding: 2 }[which] || 0;
        const pick = values => (values ? values[setIndex] : null);   // [all, surface, binding]
        const ofLoop = residuesOf(which, 'loop residues');
        if (property === 'flexibility') {
            return {
                title: loopTitle(which, 'flexibility'),
                sub: `Mean CABS-flex RMSF of ${ofLoop} (coarse-grained simulation of each structure)`,
                chart: distributionChart({
                    panels: loopPanels((p, loop) => pick(loopValues(p, loop).rmsf), true),
                    format: v => v.toFixed(2), unit: 'Å', title: 'Mean RMSF (Å)',
                }),
            };
        }
        if (property === 'hydrophobicity') {
            const index = data.hydrophobicity_scales.findIndex(x => x.name === scale);
            return {
                title: loopTitle(which, 'hydrophobicity'),
                sub: scaleSub(scale, ofLoop),
                chart: distributionChart({
                    panels: loopPanels((p, loop) => pick((loopValues(p, loop).hydrophobicity || [])[index]), true),
                    format: v => signed(v, 2), unit: '', title: isGravy(scale) ? 'GRAVY' : `Mean ${scale}`,
                }),
            };
        }
        if (property === 'charge') {
            return {
                title: loopTitle(which, 'charge'),
                sub: which === 'all'
                    ? 'Net charge of the loop residues (atomic charges from pdb2pqr, AMBER)'
                    : `Net charge of ${ofLoop}`,
                chart: distributionChart({
                    panels: loopPanels((p, loop) => pick(loopValues(p, loop).charge), true,
                                       { binWidth: 1, zeroline: true }),
                    format: charge, unit: 'e', title: 'Charge (e)',
                }),
            };
        }
        return {
            title: which === 'all' ? 'Loop lengths'
                : which === 'surface' ? 'Surface-exposed residues per loop'
                    : 'Substrate-binding residues per loop',
            sub: (which === 'all'
                ? 'Residues per loop region, the segments between the core β-strands'
                : `How many residues of each loop are ${which === 'surface'
                    ? `surface-exposed (relative SASA ≥ ${data.surface_rsa} %)`
                    : 'substrate-binding: touched by a plane laid on the copper face'}`)
                + '; a family without dots or bars in a panel lacks that loop',
            chart: distributionChart({
                // binding-surface counts run 0-15, so one bar per residue, the same width in
                // every panel; loop lengths (to 113) and surface counts (to 50) keep the width
                // niceBinWidth picks for their panel
                panels: loopPanels((p, loop) => (which === 'all' ? p.loops[loop][0]
                    : pick(loopValues(p, loop).residues)),
                    which !== 'all', which === 'binding' ? { binWidth: 1 } : {}),
                format: residues, unit: 'residues', title: 'Residues',
            }),
        };
    }
    // The second-sphere card (pipeline step 8c-8d): the motif (H-Q-Y), the motif with its
    // spacing under it (H-Q-Y › H-X8-Q-X-Y), or the shielding residue (on the copper's solvent side)
    const LAST_GROUPS = [data.open_position, data.not_determined];
    const countsOf = values => {
        const counts = new Map();
        values.forEach(v => counts.set(v, (counts.get(v) || 0) + 1));
        return [...counts].sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]));
    };
    function motifGroups(withSpacing) {
        const shortGroups = groupsBy(p => [p.motif[0]], LAST_GROUPS);
        if (!withSpacing) {
            return shortGroups.map(g => {
                const variants = countsOf(g.proteins.map(p => p.motif[1]).filter(Boolean));
                return {
                    ...g,
                    detail: variants.length > 1 || (variants.length && variants[0][0] !== g.label)
                        ? variants.map(([motif, n]) => `${motif} ${n}`).join(' · ') : '',
                };
            });
        }
        return shortGroups.flatMap(g => {
            const variants = countsOf(g.proteins.map(p => p.motif[1])).map(([motif]) => motif);
            const restOnly = countsOf(g.rest.map(p => p.motif[1])).map(([motif]) => motif)
                .filter(motif => !variants.includes(motif));
            return variants.concat(restOnly).map(motif => ({
                label: motif ? `${g.label} · ${motif}` : g.label,
                levels: [g.label, motif || ' '],
                proteins: g.proteins.filter(p => p.motif[1] === motif),
                rest: g.rest.filter(p => p.motif[1] === motif),
            }));
        });
    }
    function secondSphere({ show }) {
        const height = PLOT_HEIGHTS.xtall;   // room for the ~20 rows of the spaced motifs
        if (show === 'distal') {
            return {
                title: 'Shielding residue',
                sub: 'Residue over the copper on its solvent side, opposite the axial aromatic (an atom '
                    + 'within 4.5 Å of the copper and 65° of the axis), e.g. the conserved Ala of AA10s',
                chart: groupChart(groupsBy(p => [p.distal], LAST_GROUPS).map(g => ({
                    ...g,
                    detail: g.label === data.open_position ? 'No residue covers this side of the copper' : '',
                })), 'Residue', height),
            };
        }
        return {
            title: show === 'spacing' ? 'Second-sphere motif with spacing' : 'Second-sphere motif',
            sub: show === 'spacing'
                ? 'Xn: n residues between two motif residues; …: the Gln/Glu lies far before the rest'
                : 'Hydrogen-bonding residues from the one holding the Gln/Glu to the axial aromatic '
                    + '(H-Q-Y = H-Xn-Q-X-Y); E-N-F: the Gln/Glu closest to the copper, then N-X-F',
            chart: groupChart(motifGroups(show === 'spacing'), 'Motif', height),
        };
    }

    // The brace card (pipeline step 8e): the motif after the N-terminal His (N-His) or after
    // the internal His (Internal-His); the tooltip and the table name the residues behind it
    const BRACE_LAST = [data.other_motif, data.not_determined];
    const BRACE_VARIANTS = 8;   // residue variants named in a tooltip
    function braceMotifs({ his }) {
        const groups = groupsBy(p => [p.brace[his][0]], BRACE_LAST).map(g => {
            const variants = countsOf(g.proteins.map(p => p.brace[his][1]).filter(Boolean))
                .map(([residues, n]) => `${residues} ${n}`);
            const more = variants.length - BRACE_VARIANTS;
            return {
                ...g,
                detail: variants.slice(0, BRACE_VARIANTS).join(' · ') + (more > 0 ? ` · ${more} more` : ''),
                residues: variants.join(' · '),
            };
        });
        const column = { head: 'Residues', value: g => g.residues };
        if (his === 'internal_his') {
            return {
                title: 'Motif after the internal His',
                sub: 'The Gly of H-X-G-P (AA9, AA17) and H-X-G lies on β8, too close for a side chain; '
                    + 'tooltip: the His and the residues after it',
                chart: groupChart(groups, 'Motif', PLOT_HEIGHTS.normal, column),
            };
        }
        return {
            title: 'Motif after the N-terminal His',
            sub: 'Xn: n residues of any kind; the Gly lies on β4, the Arg bridges to an Asp of β9; '
                + 'tooltip: the His and the residues after it',
            chart: groupChart(groups, 'Motif', PLOT_HEIGHTS.normal, column),
        };
    }

    // Disulfide bonds with both cysteines inside the catalytic domain (overview.structural_features:
    // p.features.disulfides = their number, p.features.disulfide_pattern = the regions each joins,
    // 'L2–β9 + other–LC', 'none' without a bond): how many a protein has, or the pairs of regions
    // joined -- the commonest as bars of their own (a protein counts in each pair it has), the rest
    // together in one bar.
    const DISULFIDE_PAIRS = 8;      // pairs of regions shown as their own bar
    function disulfides({ show }) {
        const sub = 'Both cysteines in the catalytic domain, sulfur atoms within 3.0 Å';
        if (show === 'count') {
            const groups = groupsBy(p => {
                const n = (p.features || {}).disulfides;
                return [n == null ? data.not_determined : String(n)];
            }).sort((a, b) => (a.label === data.not_determined) - (b.label === data.not_determined)
                || Number(a.label) - Number(b.label));
            return { title: 'Disulfide bonds', sub, chart: groupChart(groups, 'Bonds', PLOT_HEIGHTS.normal) };
        }
        const NONE = 'No bond';
        const OTHER = 'Other pairs';
        const pairsOf = p => {
            const pattern = (p.features || {}).disulfide_pattern;
            if (pattern == null) return [data.not_determined];
            return pattern === 'none' ? [NONE] : [...new Set(pattern.split(' + '))];
        };
        const all = groupsBy(pairsOf, [NONE, data.not_determined]);
        const named = all.filter(g => g.label !== NONE && g.label !== data.not_determined);
        const shown = new Set(named.slice(0, DISULFIDE_PAIRS).map(g => g.label));
        const groups = groupsBy(p => [...new Set(pairsOf(p).map(pair =>
            (shown.has(pair) || pair === NONE || pair === data.not_determined ? pair : OTHER)))],
            [OTHER, NONE, data.not_determined])
            .map(g => (g.label === OTHER
                ? { ...g, detail: named.slice(DISULFIDE_PAIRS).map(o => `${o.label} ${o.proteins.length}`).join(' · ') }
                : g));
        return {
            title: 'Regions joined by disulfide bonds',
            sub: 'A protein counts in every pair it has; other: outside the loops and core strands',
            chart: groupChart(groups, 'Regions', PLOT_HEIGHTS.normal, { head: 'Pairs', value: g => g.detail || '' }),
        };
    }

    // CBM families (overview: p.cbms = {uniprot, structure}) as the protein's UniProt entries
    // annotate them, or only those inside the catalogue's structure (map_domain_ranges); the
    // first bar holds every protein with one, so one click selects them all
    const ANY_CBM = 'Any CBM';
    function cbmModules({ source }) {
        const other = source === 'uniprot' ? 'structure' : 'uniprot';
        const has = (p, which, family) => (family ? p.cbms[which].includes(family) : p.cbms[which].length > 0);
        const withCbm = proteins.filter(p => has(p, source));
        const more = n => `${n} more ${n === 1 ? 'protein has' : 'proteins have'} one outside the structure`;
        // how the other source sees the bar's proteins
        const detail = family => {
            if (source === 'uniprot') {
                const n = proteins.filter(p => has(p, source, family)).length;
                const inside = proteins.filter(p => has(p, source, family) && has(p, other, family)).length;
                const where = 'inside the catalogue\'s structure';
                if (inside === n) return n === 1 ? `Lies ${where}` : `All ${where}`;
                if (!inside) return n === 1 ? `Not ${where}` : `None ${where}`;
                return `${inside} of them ${where}`;
            }
            const outside = proteins.filter(p => !has(p, source, family) && has(p, other, family)).length;
            return outside ? `${more(outside)} (UniProt)` : '';
        };
        const restWithCbm = rest.filter(p => has(p, source));
        const groups = (withCbm.length || restWithCbm.length
            ? [{ label: ANY_CBM, proteins: withCbm, rest: restWithCbm, detail: detail(null) }] : [])
            .concat(groupsBy(p => p.cbms[source]).map(g => ({ ...g, detail: detail(g.label) })));
        const outside = proteins.filter(p => !has(p, source) && has(p, other)).length;
        return {
            title: 'Carbohydrate-binding modules',
            sub: source === 'uniprot'
                ? `${withCbm.length} of ${proteins.length} proteins have one in their UniProt entries, also where `
                    + 'the catalogue structure holds only the catalytic domain'
                : `${withCbm.length} of ${proteins.length} catalogue structures contain one`
                    + (outside ? `; ${more(outside)}` : ''),
            chart: groupChart(groups, 'Module', PLOT_HEIGHTS.normal, { head: 'Note', value: g => g.detail || '' }),
        };
    }

    // The charts; rebuilt when the compare page's switch shows or hides the rest of the catalogue
    const buildCharts = () => ({
        families: familiesChart(PLOT_HEIGHTS.tall),
        taxonomy: groupChart(taxonomyGroups(), 'Taxon', PLOT_HEIGHTS.tall),
        substrates: configurable('substrates', { usesScope: false, controls: [SHOW_CONTROL], build: substrates }),
        regioselectivity: configurable('regioselectivity', {
            usesScope: false, controls: [SHOW_CONTROL], build: regioselectivity,
        }),
        cbms: configurable('cbms', {
            usesScope: false,
            controls: [{
                key: 'source', label: 'CBMs of', type: 'buttons',
                options: [
                    { value: 'uniprot', label: 'Protein (UniProt)',
                      title: 'Modules the protein\'s UniProt entries annotate, whether or not its catalogue structure holds them' },
                    { value: 'structure', label: 'In the structure',
                      title: 'Modules inside the structure the catalogue shows (crystal structures often hold only the catalytic domain)' },
                ],
            }],
            build: cbmModules,
        }),
        disulfides: configurable('disulfides', {
            usesScope: false,
            controls: [{
                key: 'show', label: 'Disulfide bonds', type: 'buttons',
                options: [
                    { value: 'count', label: 'Number' },
                    { value: 'regions', label: 'Regions joined',
                      title: 'The loops and core strands the two cysteines of each bond lie in' },
                ],
            }],
            build: disulfides,
        }),
        motifs: configurable('motifs', {
            usesScope: false,
            controls: [{
                key: 'show', label: 'Second sphere', type: 'buttons',
                options: [
                    { value: 'motif', label: 'Motif' },
                    { value: 'spacing', label: 'With spacing' },
                    { value: 'distal', label: 'Shielding residue' },
                ],
            }],
            build: secondSphere,
        }),
        brace: configurable('brace', {
            usesScope: false,
            controls: [{
                key: 'his', label: 'Histidine', type: 'buttons',
                options: [
                    { value: 'n_his', label: 'N-His' },
                    { value: 'internal_his', label: 'Internal-His' },
                ],
            }],
            build: braceMotifs,
        }),
        size: configurable('size', {
            controls: [
                {
                    key: 'measure', label: 'Measure', type: 'buttons',
                    options: [{ value: 'length', label: 'Length' }, { value: 'mass', label: 'Mass' }],
                },
                { key: 'residues', label: 'Residues', type: 'buttons', options: RESIDUE_OPTIONS },
            ],
            build: ({ scope: sc, measure, residues: which }) => {
                const mass = measure === 'mass';
                const inDomain = Boolean(RESIDUE_SET[which]);
                // loops and binding carry one value (measures.sets), the others a scope pair
                const value = p => (which === 'all' ? p.measures[mass ? 'mass' : 'length'][sc]
                    : which === 'surface' ? (p.measures[mass ? 'surface_mass' : 'surface_length'] || {})[sc]
                        : setOf(p, which)[mass ? 'mass' : 'residues']);
                const title = which === 'all'
                    ? (mass ? 'Mass' : 'Length')
                    : mass ? setTitle(which, 'mass')
                        : which === 'surface' ? 'Surface residues'
                            : which === 'loops' ? 'Residues in the loop regions'
                                : 'Substrate-binding residues';
                const sub = which === 'all'
                    ? {
                        'length cd': 'Residues up to the catalytic domain end (no CBMs, linkers or tails)',
                        'length all': 'Residues of the whole structure, with any CBMs, linkers and tails in it',
                        'mass cd': 'Polypeptide mass up to the catalytic domain end (no glycans or metal)',
                        'mass all': 'Polypeptide mass of the whole structure (no glycans or metal)',
                    }[`${measure} ${sc}`]
                    : (mass ? `Summed residue mass of ${residuesOf(which, 'residues')} (no glycans or metal)`
                        : `Number of ${residuesOf(which, 'residues').replace(/^the /, '')}`)
                        + (inDomain ? insideDomain : '');
                // One bar per residue for the binding surface (about 10-25 residues); tens for the
                // larger sets, so their histograms keep a readable number of bars.
                const binWidth = mass ? undefined
                    : which === 'binding' ? 1
                        : (which === 'all' && sc === 'all') ? undefined : 10;
                return {
                    title: title + (inDomain ? '' : ` (${SCOPES[sc]})`),
                    sub,
                    chart: distributionChart(mass
                        ? { panels: [{ name: 'Mass', value }],
                            format: v => v.toFixed(1), unit: 'kDa', title: 'Mass (kDa)' }
                        : { panels: [{ name: 'Length', value, binWidth }],
                            format: residues, unit: 'residues', title: 'Residues' }),
                };
            },
        }),
        hydrophobicity: configurable('hydrophobicity', {
            controls: [
                { key: 'scale', label: 'Hydrophobicity scale', type: 'select', options: SCALE_OPTIONS },
                { key: 'residues', label: 'Residues', type: 'buttons', options: RESIDUE_OPTIONS },
            ],
            build: ({ scope: sc, scale, residues: which }) => {
                const gravy = isGravy(scale);
                const of = residuesOf(which, 'residues');
                return {
                    title: setTitle(which, 'hydrophobicity')
                        + (RESIDUE_SET[which] ? '' : ` (${SCOPES[sc]})`),
                    sub: scaleSub(scale, of) + (RESIDUE_SET[which] ? insideDomain : ''),
                    chart: distributionChart({
                        panels: [{
                            name: gravy ? 'GRAVY' : scale,
                            value: p => hydrophobicityValue(p, which, scale, sc),
                        }],
                        format: v => signed(v, 2), unit: '', title: gravy ? 'GRAVY' : `Mean ${scale}`,
                    }),
                };
            },
        }),
        charge: configurable('charge', {
            controls: [{ key: 'residues', label: 'Residues', type: 'buttons', options: RESIDUE_OPTIONS }],
            build: ({ scope: sc, residues: which }) => ({
                title: (which === 'all' ? 'Net charge' : setTitle(which, 'charge'))
                    + (RESIDUE_SET[which] ? '' : ` (${SCOPES[sc]})`),
                sub: (which === 'all'
                    ? 'Sum of the atomic charges (pdb2pqr, AMBER), copper excluded'
                    : `Net charge of ${residuesOf(which, 'residues')}`)
                    + (RESIDUE_SET[which] ? insideDomain : ''),
                chart: distributionChart({
                    panels: [{
                        name: 'Charge', zeroline: true,
                        value: p => chargeValue(p, which, sc),
                    }],
                    format: charge, unit: 'e', title: 'Charge (e)',
                }),
            }),
        }),
        loops: proteins.concat(rest).some(p => Object.keys(p.loops).length) ? configurable('loops', {
            usesScope: false,
            controls: [
                {
                    key: 'property', label: 'Loop property', type: 'buttons',
                    options: [
                        { value: 'length', label: 'Length' },
                        { value: 'flexibility', label: 'Flexibility' },
                        { value: 'hydrophobicity', label: 'Hydrophobicity' },
                        { value: 'charge', label: 'Charge' },
                    ],
                },
                {
                    key: 'scale', label: 'Hydrophobicity scale', type: 'select', options: SCALE_OPTIONS,
                    when: v => v.property === 'hydrophobicity',
                },
                { key: 'residues', label: 'Residues', type: 'buttons', options: LOOP_RESIDUE_OPTIONS },
            ],
            build: loopProperty,
        }) : null,
    });
    let charts = buildCharts();

    // A card's chart, looked up in `charts` at each use (setupCard binds it once)
    function liveChart(key) {
        const now = () => charts[key];
        const live = {};
        ['figure', 'hover', 'barFocus', 'table', 'labels', 'dependsOnSelection', 'set', 'current'].forEach(name => {
            live[name] = (...args) => (now()[name] ? now()[name](...args) : undefined);
        });
        ['controls', 'values', 'hasViews', 'usesScope'].forEach(name => {
            Object.defineProperty(live, name, { get: () => now()[name] });
        });
        return live;
    }

    const state = loadState();
    let view = state.view === 'family' ? 'family' : 'stacked';
    let scope = state.scope === 'all' ? 'all' : 'cd';   // of the distribution charts

    function renderStats() {
        const box = document.getElementById('overviewStats');
        if (!box) return;   // compare page: the tiles are above the table (views.compare_proteins)
        const count = source => proteins.filter(p => p.source === source).length;
        const predicted = MODEL_SOURCES.map(([s, name]) => [name, count(s)]).filter(([, n]) => n);
        const tiles = [
            { value: proteins.length, label: 'Proteins',
              note: selectionMode ? `of ${data.proteins.length} in the catalogue` : '' },
            { value: familiesIn(proteins).length, label: 'LPMO families' },
            { value: new Set(proteins.map(p => p.species)).size, label: 'Source species' },
            { value: count('PDB'), label: 'Experimental structures', note: 'PDB' },
            {
                value: predicted.reduce((sum, [, n]) => sum + n, 0), label: 'Predicted models',
                note: predicted.map(([s, n]) => `${s} ${n}`).join(' · '),
            },
        ];
        box.replaceChildren(...tiles.map(t => {
            const tile = el('div', 'ov-stat');
            tile.append(el('div', 'ov-stat-value', String(t.value)), el('div', 'ov-stat-label', t.label));
            if (t.note) tile.appendChild(el('div', 'ov-stat-note', t.note));
            return tile;
        }));
    }

    function fillTable(box, table) {
        const tbl = el('table');
        const headRow = el('tr');
        table.head.forEach((h, i) => headRow.appendChild(el('th', table.numeric[i] ? 'num' : '', h)));
        tbl.appendChild(el('thead')).appendChild(headRow);
        const body = tbl.appendChild(el('tbody'));
        const addRow = (cells, className) => {
            const tr = body.appendChild(el('tr', className || ''));
            cells.forEach((c, i) => tr.appendChild(el('td', table.numeric[i] ? 'num' : '', String(c))));
        };
        table.rows.forEach(r => addRow(r));
        if (table.total) addRow(table.total, 'ov-total');
        box.replaceChildren(tbl);
    }

    function setupCard(card) {
        const chart = liveChart(card.dataset.chart);
        const plot = el('div', 'ov-plot');
        const tableBox = el('div', 'ov-table');
        tableBox.hidden = true;
        const button = el('button', 'ov-view-btn', 'Table');
        button.type = 'button';
        button.setAttribute('aria-pressed', 'false');
        button.title = 'Show the numbers as a table';
        const caption = card.querySelector('figcaption');
        const tools = el('div', 'ov-card-tools');
        caption.appendChild(tools);
        card.append(plot, tableBox);

        // settings of a configurable card: its controls, and its title and subtitle
        let syncCaption = () => {};
        // a control applies only to some settings (control.when); display, as the CSS sets it
        const shownWhen = (control, node, sync) => () => {
            node.style.display = control.when && !control.when(chart.values) ? 'none' : '';
            sync();
        };
        if (chart.controls) {
            const titleNode = caption.querySelector('.ov-card-title');
            const subNode = caption.querySelector('.ov-card-sub');
            const controlsRow = el('div', 'ov-card-controls');   // under the subtitle
            caption.appendChild(controlsRow);
            const syncers = chart.controls.map(control => {
                if (control.type === 'select') {
                    const select = el('select', 'ov-select');
                    select.setAttribute('aria-label', control.label);
                    control.options.forEach(option => {
                        const node = el('option', '', option.label);
                        node.value = option.value;
                        select.appendChild(node);
                    });
                    select.addEventListener('change', () => {
                        chart.set(control.key, select.value);
                        card.refresh();
                    });
                    controlsRow.appendChild(select);
                    return shownWhen(control, select, () => { select.value = chart.values[control.key]; });
                }
                const group = el('div', 'ov-segmented');
                group.setAttribute('role', 'group');
                group.setAttribute('aria-label', control.label);
                const buttons = control.options.map(option => {
                    const b = el('button', '', option.label);
                    b.type = 'button';
                    b.dataset.value = option.value;
                    if (option.title) b.title = option.title;
                    b.addEventListener('click', () => {
                        chart.set(control.key, option.value);
                        card.refresh();
                    });
                    return b;
                });
                group.append(...buttons);
                controlsRow.appendChild(group);
                return shownWhen(control, group, () => buttons.forEach(b => {
                    const active = b.dataset.value === chart.values[control.key];
                    b.classList.toggle('active', active);
                    b.setAttribute('aria-pressed', String(active));
                }));
            });
            syncCaption = () => {
                syncers.forEach(sync => sync());
                const { title, sub, subHtml } = chart.current();
                titleNode.textContent = title;
                if (subHtml) subNode.innerHTML = subHtml;   // our own markup: the colour swatches
                else subNode.textContent = sub;
            };
            syncCaption();
        }
        tools.appendChild(button);
        const svgButton = el('button', 'ov-view-btn', 'SVG');
        svgButton.type = 'button';
        svgButton.title = 'Download the chart as an SVG image';
        tools.appendChild(svgButton);
        svgButton.addEventListener('click', () => {
            if (!plot.dataset.rendered || !window.SvgExport) return;
            const title = caption.querySelector('.ov-card-title').textContent;
            window.SvgExport.plotly(plot, {
                title, sub: caption.querySelector('.ov-card-sub'),
                filename: window.SvgExport.filenameOf(title, selectionMode ? 'lpmocat_selection_' : 'lpmocat_'),
            });
        });

        button.addEventListener('click', () => {
            const showTable = tableBox.hidden;
            if (showTable && !tableBox.firstChild) fillTable(tableBox, chart.table());
            tableBox.hidden = !showTable;
            plot.hidden = showTable;
            button.textContent = showTable ? 'Chart' : 'Table';
            button.title = showTable ? 'Show the chart' : 'Show the numbers as a table';
            button.setAttribute('aria-pressed', String(showTable));
            if (!showTable && plot.dataset.rendered) Plotly.Plots.resize(plot);
        });

        let hovered = null;     // {content, focus} of the mark or label under the pointer
        let press = null;       // {content, focus, x, y} while the button is down on a mark
        let lastClick = null;   // {content, focus, x, y} of the last click while the pointer stays there
        let stickUntil = 0;     // until then a coordinate-less unhover after a click is ignored
        const near = (event, point) => !!point && Math.hypot(event.clientX - point.x, event.clientY - point.y) <= CLICK_SLOP;
        // the plot, or a label of it outside the plot (a legend item in the caption)
        const onChart = node => plot.contains(node)
            || (card.contains(node) && !!(node.closest && node.closest('.ov-axis-label')));
        const clear = () => {
            hovered = null;
            hideTooltip();
            highlightInTree(null);
            fadeOtherBars(plot, null);
        };
        const focusOn = (event, content, focus) => {
            hovered = { content, focus };
            showTooltip(event, content);
            highlightInTree(membersOf(content));
            fadeOtherBars(plot, focus);
        };
        // Unhover / leave ends the hover, except while the button is down, or right after a click
        // with the pointer still on the spot: Plotly unhovers then although the pointer did not
        // move, and the next click there must still find its bar.
        const leave = event => {
            if (press) return;
            const located = event && typeof event.clientX === 'number' && (event.clientX || event.clientY);
            if (lastClick && (located ? near(event, lastClick) : Date.now() < stickUntil)) return;
            clear();
        };

        // A click selects what was under the pointer when the button went down. Plotly's own
        // plotly_click is not used: it is dropped when the pointer moves a pixel or two during
        // the click (and it unhovers), which a hand on a mouse or touchpad often does.
        card.addEventListener('mousedown', event => {
            if (!onChart(event.target)) return;
            const target = hovered || (near(event, lastClick) ? lastClick : null);
            press = event.button === 0 && target
                ? { content: target.content, focus: target.focus, x: event.clientX, y: event.clientY } : null;
        }, true);
        document.addEventListener('mouseup', event => {
            if (press && near(event, press)) {
                toggleSelection(membersOf(press.content));
                lastClick = { ...press, x: event.clientX, y: event.clientY };
                stickUntil = Date.now() + 800;
                if (!hovered) focusOn(event, press.content, press.focus);   // show the new state
            }
            press = null;
        });
        // Moving away forgets the click; outside the chart nothing stays highlighted.
        document.addEventListener('mousemove', event => {
            if (lastClick && !near(event, lastClick)) lastClick = null;
            if (hovered && !press && !onChart(event.target)) clear();
        });

        // category labels hover and click like their bars; after every redraw (plotly_afterplot)
        // each label node gets its listeners once and its current target in labelTargets
        const labelTargets = new WeakMap();
        function bindLabels() {
            (chart.labels(plot, card) || []).forEach(target => {
                const { node } = target;
                node.classList.add('ov-axis-label');
                if (!labelTargets.has(node)) {
                    node.addEventListener('mouseenter', event => {
                        const { content, focus } = labelTargets.get(node);
                        const shownContent = content();
                        if (shownContent) focusOn(event, shownContent, focus);
                    });
                    node.addEventListener('mousemove', moveTooltip);
                    node.addEventListener('mouseleave', leave);
                }
                labelTargets.set(node, target);
            });
        }

        card.render = () => {
            const { traces, layout } = chart.figure(view);
            Plotly.react(plot, traces, layout, PLOT_CONFIG);
            if (plot.dataset.rendered) return;
            plot.dataset.rendered = '1';
            plot.on('plotly_afterplot', bindLabels);
            bindLabels();
            plot.on('plotly_hover', event => {
                const content = chart.hover(event.points, view);
                if (content) focusOn(event.event, content, chart.barFocus(event.points, view));
                else clear();
            });
            plot.on('plotly_unhover', data => leave(data && data.event));
            plot.addEventListener('mousemove', moveTooltip);
            plot.addEventListener('mouseleave', leave);
        };
        card.dependsOnSelection = () => !!(chart.dependsOnSelection && chart.dependsOnSelection(view));
        // after a change of settings: caption, table and chart
        card.refresh = () => {
            syncCaption();
            if (tableBox.firstChild) fillTable(tableBox, chart.table());
            if (plot.dataset.rendered) card.render();
        };
    }

    const cards = [...root.querySelectorAll('.ov-card')].filter(card => {
        card.hidden = !charts[card.dataset.chart];   // e.g. no loop regions loaded
        return !card.hidden;
    });
    cards.forEach(setupCard);

    renderCharts = () => cards.forEach(card => {
        card.render();
        const plot = card.querySelector('.ov-plot');
        if (plot && !plot.hidden) Plotly.Plots.resize(plot);
    });

    renderStats();

    // view and region switches of the distribution charts: beside the Properties heading and,
    // as a copy, in the Length card above it; every copy is kept in step
    const viewButtons = [...root.querySelectorAll('button[data-view]')];
    function syncViewButtons() {
        viewButtons.forEach(b => {
            const active = b.dataset.view === view;
            b.classList.toggle('active', active);
            b.setAttribute('aria-pressed', String(active));
        });
    }
    const scopeButtons = [...root.querySelectorAll('button[data-scope]')];
    function syncScopeButtons() {
        scopeButtons.forEach(b => {
            const active = b.dataset.scope === scope;
            b.classList.toggle('active', active);
            b.setAttribute('aria-pressed', String(active));
        });
    }
    scopeButtons.forEach(b => b.addEventListener('click', () => {
        if (b.dataset.scope === scope) return;
        scope = b.dataset.scope;
        syncScopeButtons();
        saveState({ scope });
        cards.filter(c => (charts[c.dataset.chart] || {}).usesScope).forEach(c => c.refresh());
    }));
    syncScopeButtons();

    viewButtons.forEach(b => b.addEventListener('click', () => {
        if (b.dataset.view === view) return;
        view = b.dataset.view;
        syncViewButtons();
        saveState({ view });
        cards.filter(c => (charts[c.dataset.chart] || {}).hasViews).forEach(c => c.render());
    }));
    syncViewButtons();

    // ---- selection changes (updateDendrogramHighlighting in index.html calls this) ----
    let syncQueued = false;
    window.catalogueOverview = {
        showPane,
        selectionChanged() {
            if (syncQueued) return;
            syncQueued = true;
            requestAnimationFrame(() => {
                syncQueued = false;
                if (!root.hidden) cards.filter(c => c.dependsOnSelection()).forEach(c => c.render());
                if (shown) showTooltip(shown.event, shown.content);
            });
        },
    };

    if (selectionMode) startSelectionOverview();
    else showPane(state.tab === 'selection' ? 'selection' : 'overview', { initial: true });

    // ---- compare page ----
    // The card's fold. The charts are drawn once the card is open, near the viewport (it lies
    // below the plots) and the page has loaded its structures (hideComparisonLoadingOverlay in
    // compare_protein.html fires 'comparison-ready'), so they don't hold up the 3D viewer.
    function setupSelectionCard() {
        const card = document.getElementById('selectionOverviewCard') || root;
        const fold = document.getElementById('selectionOverviewFold');
        const tools = document.getElementById('selectionOverviewTools');
        let open = loadState().open !== false;
        let pageReady = !!window.comparisonReady;
        let seen = !('IntersectionObserver' in window);
        let onShow = null;
        function sync() {
            root.hidden = !open;
            if (tools) tools.hidden = !open;
            if (fold) fold.setAttribute('aria-expanded', String(open));
            if (open && pageReady && seen && onShow) onShow();
        }
        if (!seen) {
            const observer = new IntersectionObserver(entries => {
                if (!entries.some(entry => entry.isIntersecting)) return;
                observer.disconnect();
                seen = true;
                sync();
            }, { rootMargin: '400px' });
            observer.observe(card);
        }
        if (fold) {
            fold.addEventListener('click', () => {
                open = !open;
                saveState({ open });
                sync();
            });
        }
        if (!pageReady) {
            document.addEventListener('comparison-ready', () => {
                pageReady = true;
                sync();
            }, { once: true });
        }
        sync();
        return {
            whenShown(fn) {
                onShow = fn;
                sync();
            },
        };
    }

    // The tree of the compared proteins, the switch for the rest of the catalogue, the charts
    function startSelectionOverview() {
        const loading = root.querySelector('.sel-loading');
        if (loading) loading.remove();

        const treeBox = document.getElementById('selectionTree');
        if (treeBox && window.SelectionTree) {
            const leaves = node => (node ? (node.children ? node.children.flatMap(leaves) : [node.name]) : []);
            const inTree = new Set(leaves(payload.tree));
            // the comparison's groups (browse page, "Add as group"): a Group column in the tree
            const own = window.comparisonGroups ? window.comparisonGroups().groups : [];
            const groupByCode = own.length ? Object.fromEntries(own.flatMap(g => g.codes.map(
                code => [code, {name: g.name, color: g.color, ungrouped: !!g.ungrouped}]))) : null;
            selectionTree = window.SelectionTree.draw(treeBox, {
                tree: payload.tree, data, familyColor, missing: payload.selected.filter(code => !inTree.has(code)),
                groups: groupByCode,
            });
            const exportButton = document.getElementById('selectionTreeExport');
            if (exportButton && selectionTree) exportButton.addEventListener('click', () => selectionTree.exportSvg());
            // The sequence logo and the similarity heatmap highlight their proteins here, and the
            // heatmap follows the tree's row order.
            window.selectionTreeInstance = selectionTree;
            document.dispatchEvent(new Event('selection-tree-ready'));

            // Distance cut: the slider runs 0..max, the far right meaning no grouping; the number
            // field takes an exact distance (empty = no grouping). Each follows the other.
            const cutInput = document.getElementById('selectionTreeCut');
            const cutValue = document.getElementById('selectionTreeCutValue');
            const cutRead = document.getElementById('selectionTreeCutRead');
            // "Which features separate the groups?": scored only on request, marked stale on a new cut
            const grouping = window.GroupingFeatures && window.GroupingFeatures.attach({
                root: document.getElementById('groupingFeatures'),
                button: document.getElementById('groupingFeaturesRun'),
                treeButton: document.getElementById('treeFeaturesRun'),
                pickButton: document.getElementById('pickedFeaturesRun'),
                referenceSelect: document.getElementById('pickedFeaturesReference'),
                // the comparison's own groups (browse page, "Add as group"), when it has them
                sourceSelect: document.getElementById('groupingSource'),
                userGroups: window.comparisonGroups ? window.comparisonGroups().groups : [],
                status: document.getElementById('groupingFeaturesStatus'),
                tree: selectionTree, data, compared: payload.selected,
            });
            const clearPick = document.getElementById('pickedFeaturesClear');
            if (clearPick && selectionTree) clearPick.addEventListener('click', () => selectionTree.clearPicked());
            if (cutInput && selectionTree && selectionTree.maxDistance) {
                const max = selectionTree.maxDistance;
                const decimals = max < 1 ? 2 : 1;
                const steps = +cutInput.max;
                if (cutValue) {
                    cutValue.max = max.toFixed(decimals + 1);
                    cutValue.step = decimals === 2 ? '0.01' : '0.1';
                }
                // `from` is the control the value came from, so it is not rewritten under the cursor
                const setCut = (cut, from) => {
                    const state = selectionTree.setCutoff(cut);
                    if (from !== 'slider') {
                        cutInput.value = state.cut == null ? steps : Math.round(state.cut / max * steps);
                    }
                    if (cutValue && from !== 'field') {
                        cutValue.value = state.cut == null ? '' : state.cut.toFixed(decimals + (state.cut < 1 ? 1 : 0));
                    }
                    cutRead.textContent = state.cut == null ? 'off' : `${state.groups} groups`;
                    if (grouping) grouping.sync(state.cut);
                };
                cutInput.addEventListener('input', () => {
                    const frac = +cutInput.value / steps;
                    setCut(frac >= 1 ? null : frac * max, 'slider');
                });
                if (cutValue) {
                    cutValue.addEventListener('input', () => {
                        const text = cutValue.value.trim();
                        // a half-typed number ("0.") reads as '' but flags badInput: wait for more
                        if (text === '' && cutValue.validity.badInput) return;
                        if (text === '') return setCut(null, 'field');
                        const cut = parseFloat(text);
                        if (!Number.isFinite(cut)) return;          // half-typed ("0.") : wait
                        setCut(Math.max(0, cut), 'field');
                    });
                    // on Enter or leaving the field, show the value the tree actually uses
                    cutValue.addEventListener('change', () => {
                        const cut = parseFloat(cutValue.value);
                        setCut(cutValue.value.trim() === '' || !Number.isFinite(cut) ? null : Math.max(0, cut));
                    });
                }
                setCut(null);
            } else if (cutInput) {
                cutInput.disabled = true;
                if (cutValue) cutValue.disabled = true;
            }
        }

        // With the rest of the catalogue in grey, or the compared proteins only
        const restButtons = [...document.querySelectorAll('#selectionOverviewTools [data-rest]')];
        const syncRest = () => restButtons.forEach(b => {
            const active = (b.dataset.rest === 'shown') === showRest;
            b.classList.toggle('active', active);
            b.setAttribute('aria-pressed', String(active));
        });
        if (!notCompared.length) {   // the whole catalogue is compared
            restButtons.forEach(b => { b.parentElement.hidden = true; });
        }
        restButtons.forEach(b => b.addEventListener('click', () => {
            const show = b.dataset.rest === 'shown';
            if (show === showRest) return;
            showRest = show;
            rest = showRest ? notCompared : [];
            saveState({ rest: showRest ? 'shown' : 'hidden' });
            syncRest();
            hideTooltip();
            highlightInTree(null);
            charts = buildCharts();
            cards.forEach(card => {
                card.hidden = !charts[card.dataset.chart];
                if (!card.hidden) card.refresh();
            });
        }));
        syncRest();

        selectionCard.whenShown(renderCharts);
    }
})();
