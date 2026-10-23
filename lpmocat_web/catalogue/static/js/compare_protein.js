async function checkInstanceReps(instance, reps, interval = 1000) {
    //function checks if instance has visible added representations; works only on initial loading
    return new Promise(resolve => {
        let elapsed = 0;
        const maxWait = 30000; // 30s timeout in case empty reps don't count toward addedReprs
        const intervalId = setInterval(() => {
            elapsed += interval;
            if (Object.keys(instance.addedReprs).length >= reps || elapsed >= maxWait) {
                clearInterval(intervalId);
                resolve(true);
            };
        }, interval);
    });
};

async function waitForRepsReady(instance, reps, interval = 500) {
    // Wait for select() to process: first wait for rep count to drop (reps being rebuilt),
    // then wait for it to reach the expected count again.
    // If count never drops (select was instant), fall back to a stabilization check.
    return new Promise(resolve => {
        let elapsed = 0;
        let dropped = false;
        let stableCount = 0;
        let lastCount = Object.keys(instance.addedReprs).length;
        const maxWait = 30000;
        const intervalId = setInterval(() => {
            elapsed += interval;
            let currentCount = Object.keys(instance.addedReprs).length;

            if (!dropped && currentCount < lastCount) {
                dropped = true;
            }

            if (dropped && currentCount >= reps) {
                // Reps dropped and came back — ready
                clearInterval(intervalId);
                resolve(true);
                return;
            }

            // If count never drops, wait for it to stabilize (same count for 2 consecutive checks)
            if (!dropped && currentCount >= reps) {
                stableCount++;
                if (stableCount >= 2) {
                    clearInterval(intervalId);
                    resolve(true);
                    return;
                }
            } else {
                stableCount = 0;
            }

            lastCount = currentCount;

            if (elapsed >= maxWait) {
                clearInterval(intervalId);
                resolve(true);
            }
        }, interval);
    });
};


// Large comparisons (more than COMPARE_PLOTS_CLOSED_ABOVE proteins, window.ONE_PLOT_SECTION):
// one plot section open at a time. Every protein's row is its own Plotly chart and a hidden
// section kept them all, so with 79 proteins a third open section ran the tab out of memory
// (server test, 2026-10-05). Closing a section frees its charts; reopening redraws them, the
// data coming from the server's cache.
function freePlotSection(section) {
    const scrollable = section.querySelector('.scrollable-plots-column');
    if (scrollable) {
        scrollable.querySelectorAll('.js-plotly-plot').forEach(function (div) {
            try { Plotly.purge(div); } catch (e) {}
        });
        scrollable.innerHTML = '';
    }
    section.setAttribute('data-loaded', '0');
    section.setAttribute('data-initial-loaded', '0');
}
window.freePlotSection = freePlotSection;

function dropdown(id) {
    const section = document.getElementById(id);
    if (!section) return;
    
    const container = section.closest('.dropdownPlots');
    const plotType = container ? container.getAttribute('data-plot-type') : null;
    
    const initialLoaded = section.getAttribute('data-initial-loaded') === '1';
    
    const loaded = section.getAttribute('data-loaded') === '1';
    
    const hasPlots = section.querySelector('.plotWrapper') !== null;
    
    const isCurrentlyVisible = section.style.display !== 'none';
    section.style.display = isCurrentlyVisible ? "none" : "block";
    if (window.ONE_PLOT_SECTION) {
        if (isCurrentlyVisible) {
            freePlotSection(section);
        } else {
            document.querySelectorAll('.dropdownPlots .insideDropdown').forEach(function (other) {
                if (other !== section && other.style.display !== 'none') {
                    other.style.display = 'none';
                    freePlotSection(other);
                }
            });
        }
    }

    if (section.style.display !== 'none' && !initialLoaded && !loaded && !hasPlots && plotType) {
        const scrollable = section.querySelector('.scrollable-plots-column');
        if (!scrollable) return;
        // Mark loaded up-front so a rapid open/close/open can't fire two fetches.
        section.setAttribute('data-loaded', '1');
        loadSectionPlots(section, scrollable, plotType);
    }
    // Opening/closing a section changes how far the plots extend, so re-evaluate the
    // floating horizontal scrollbar (after the display change has laid out).
    if (window.updateFloatingHScroll) {
        [0, 300].forEach(function (t) { setTimeout(updateFloatingHScroll, t); });
    }
};

// Re-fetch + re-render every currently-loaded plot section (used when the "Full length"
// toggle flips, so all plots switch between catalytic-domain and full-length together).
// Sections not yet opened simply pick up the new setting when they load.
function reloadLoadedPlotSections() {
    document.querySelectorAll('.dropdownPlots').forEach(function (container) {
        const section = container.querySelector('.insideDropdown');
        const scrollable = container.querySelector('.scrollable-plots-column');
        const plotType = container.getAttribute('data-plot-type');
        if (!section || !scrollable || !plotType) return;
        if (!scrollable.querySelector('.plotWrapper')) return;   // not loaded yet — skip
        const scale = (plotType === 'Hydrophobicity') ? window.currentHydrophobScale
                    : (plotType === 'Flexibility') ? window.currentFlexMethod : null;
        loadSectionPlots(section, scrollable, plotType, scale);
    });
}
window.reloadLoadedPlotSections = reloadLoadedPlotSections;

// Inject the loading-spinner CSS once (keeps it self-contained in this file).
function ensureLoadingSpinnerCss() {
    if (document.getElementById('lpmo-spinner-css')) return;
    const s = document.createElement('style');
    s.id = 'lpmo-spinner-css';
    s.textContent =
        '@keyframes lpmoSpin{to{transform:rotate(360deg)}}' +
        '.lpmo-loading{display:flex;align-items:center;gap:10px;padding:16px;color:#5a7a6a;font-size:0.9em}' +
        '.lpmo-spinner{width:20px;height:20px;border:3px solid #d7e6dd;border-top-color:#569976;' +
        'border-radius:50%;animation:lpmoSpin .8s linear infinite;flex:0 0 auto}';
    document.head.appendChild(s);
}

// ── Full-page loading lock ─────────────────────────────────────────────────────────
// While plot sections are (re)loading — especially the Full-length toggle over many
// proteins — block interaction with a blurred overlay + spinner, so the user can't click
// something mid-load and tip a heavy render into a crash. Keyed off an in-flight counter so
// it stays up for the whole batch; shown only after a short delay so fast cached loads don't
// flash it.
let _plotLoadCount = 0;
let _plotLockShowTimer = null;
function _ensurePlotLockOverlay() {
    let ov = document.getElementById('plotLoadingLock');
    if (ov) return ov;
    ensureLoadingSpinnerCss();
    ov = document.createElement('div');
    ov.id = 'plotLoadingLock';
    ov.style.cssText = 'position:fixed;inset:0;z-index:2147483000;display:none;align-items:center;'
        + 'justify-content:center;background:rgba(247,251,249,0.55);'
        + 'backdrop-filter:blur(2px);-webkit-backdrop-filter:blur(2px);cursor:progress;';
    ov.innerHTML = '<div class="lpmo-loading" style="background:#fff;border:1px solid #d7e6dd;'
        + 'border-radius:10px;padding:14px 22px;box-shadow:0 4px 18px rgba(0,0,0,0.12);font-size:0.95em;">'
        + '<span class="lpmo-spinner"></span>Updating plots…</div>';
    // Swallow any interaction while the lock is up.
    ['click', 'mousedown', 'wheel', 'keydown', 'touchstart'].forEach(function (ev) {
        ov.addEventListener(ev, function (e) { e.stopPropagation(); e.preventDefault(); }, true);
    });
    document.body.appendChild(ov);
    return ov;
}
function beginPlotLoad() {
    _plotLoadCount++;
    if (_plotLoadCount === 1 && !_plotLockShowTimer) {
        _plotLockShowTimer = setTimeout(function () {
            _plotLockShowTimer = null;
            if (_plotLoadCount > 0) _ensurePlotLockOverlay().style.display = 'flex';
        }, 180);
    }
}
function endPlotLoad() {
    _plotLoadCount = Math.max(0, _plotLoadCount - 1);
    if (_plotLoadCount === 0) {
        if (_plotLockShowTimer) { clearTimeout(_plotLockShowTimer); _plotLockShowTimer = null; }
        const ov = document.getElementById('plotLoadingLock');
        if (ov) ov.style.display = 'none';
    }
}
window.beginPlotLoad = beginPlotLoad;
window.endPlotLoad = endPlotLoad;

// --- Full-length view guard ------------------------------------------------
// The full-length view is cheap in 3D but expensive in the PLOTS: Plotly has to
// draw one trace per protein across the whole (much longer) residue axis, and for
// a large selection with CBM tails this can briefly freeze the main thread while
// it renders. We don't block it — we warn once the summed full-length residue
// count of the *visible* proteins crosses this budget, and let the user proceed
// or back out. Tune this number if the warning fires too eagerly / too rarely.
window.FULL_LENGTH_PLOT_WARN_BUDGET = 9000;

// Total residues (and protein count) the full-length plots would render for the
// currently-visible proteins. Uses the same per-protein residue lists Mol* uses.
function fullLengthPlotLoad() {
    const mc = window.molstarColoring;
    if (!mc || !mc.residueNumbers) return { residues: 0, proteins: 0 };
    const protCount = mc.residueNumbers.length;
    const vis = window.getProteinVisibility ? window.getProteinVisibility(protCount) : null;
    let residues = 0, proteins = 0;
    for (let i = 0; i < protCount; i++) {
        if (vis && !vis[i]) continue;
        residues += (mc.residueNumbers[i] || []).length;
        proteins++;
    }
    return { residues, proteins };
}
window.fullLengthPlotLoad = fullLengthPlotLoad;

// onchange handler for the Full-length / Catalytic-domain switch. Warns (but
// allows) before rendering a heavy full-length plot set, then applies the view
// to the label, flexibility colouring, 3D representation and open plot sections.
function handleFullLengthToggle(checkbox, protCount) {
    if (checkbox.checked) {
        const load = fullLengthPlotLoad();
        if (load.residues > (window.FULL_LENGTH_PLOT_WARN_BUDGET || 9000)) {
            const ok = window.confirm(
                'Full length for ' + load.proteins + ' proteins makes the plots heavy to '
                + 'render — the page may briefly freeze while they load.\n\n'
                + 'Tip: deselect some proteins to make this faster.\n\n'
                + 'Show full length anyway?');
            if (!ok) {
                checkbox.checked = false;   // stay on the catalytic domain
                return;
            }
        }
    }
    if (window.syncFullLenLabel) window.syncFullLenLabel();
    // the flexibility colouring of the other view may have to be fetched first
    Promise.resolve(window.setFlexColoringForView ? window.setFlexColoringForView() : null).then(function () {
        if (window.viewerInstance) changeVisibleRepresentationType(protCount, window.viewerInstance);
    });
    if (window.reloadLoadedPlotSections) window.reloadLoadedPlotSections();
}
window.handleFullLengthToggle = handleFullLengthToggle;

// --- Floating horizontal scrollbar -----------------------------------------
// The plots scroll left/right inside each `.scrollable-plots-column`, but with many
// proteins the column is very tall, so its native scrollbar sits far below the fold.
// This adds a scrollbar pinned to the bottom of the viewport that mirrors the plots'
// width and drives all the (horizontally-synced) columns. It
// only shows when the plots overflow horizontally, one of them is on screen AND their own
// scrollbar is off-screen.
let _hsSyncing = false;          // guards the floating <-> columns scroll-position loop
let _hsRepositionRaf = null;
let _hsMembers = [];             // what the bar scrolls right now: the plot columns or the table

function _ensureFloatingHScroll() {
    let fs = document.getElementById('floatingHScroll');
    if (fs) return fs;
    fs = document.createElement('div');
    fs.id = 'floatingHScroll';
    const inner = document.createElement('div');
    inner.id = 'floatingHScrollInner';
    fs.appendChild(inner);
    // Dragging the floating bar scrolls what it currently stands for: every plot column (they
    // scroll together) or the proteins table.
    fs.addEventListener('scroll', function () {
        if (_hsSyncing) return;
        _hsSyncing = true;
        const sl = fs.scrollLeft;
        _hsMembers.forEach(function (c) {
            if (c.scrollLeft !== sl) c.scrollLeft = sl;
        });
        window.requestAnimationFrame(function () { _hsSyncing = false; });
    });
    document.body.appendChild(fs);
    return fs;
}

// The areas that scroll sideways and can own the floating bar: the plot columns (one bar for
// all of them, they are kept in sync) and the proteins table, which is far taller than the
// window once "Show: All" is picked.
function _hScrollGroups() {
    const groups = [];
    const columns = [].slice.call(document.querySelectorAll('.scrollable-plots-column'))
        .filter(function (c) { return c.offsetParent !== null && (c.scrollWidth - c.clientWidth) > 2; });
    if (columns.length) groups.push(columns);
    const table = document.querySelector('.table-container');
    if (table && table.offsetParent !== null && (table.scrollWidth - table.clientWidth) > 2) {
        groups.push([table]);
    }
    return groups;
}

// Position/size the floating bar over the area being read and decide whether to show it: the
// plot columns or the proteins table, whichever fills more of the window right now.
function updateFloatingHScroll() {
    const fs = _ensureFloatingHScroll();
    const vh = window.innerHeight;
    let best = null;
    _hScrollGroups().forEach(function (members) {
        let maxBottom = -Infinity, ref = null, onScreen = 0;
        members.forEach(function (c) {
            const r = c.getBoundingClientRect();
            if (r.bottom > maxBottom) maxBottom = r.bottom;
            if (!ref && r.top < vh && r.bottom > 0) ref = c;   // first one intersecting the viewport
            onScreen += Math.max(0, Math.min(r.bottom, vh) - Math.max(r.top, 0));
        });
        // Nothing of this group on screen (e.g. reading the selection overview above the plots),
        // or its own scrollbar is already visible: it does not need the floating bar.
        if (!ref || maxBottom <= vh - 2) return;
        if (!best || onScreen > best.onScreen) best = { ref: ref, members: members, onScreen: onScreen };
    });
    if (!best) { fs.style.display = 'none'; _hsMembers = []; return; }
    const ref = best.ref;
    _hsMembers = best.members;
    const r = ref.getBoundingClientRect();
    const wasHidden = (fs.style.display === 'none' || !fs.style.display);
    fs.style.left = r.left + 'px';
    fs.style.width = r.width + 'px';
    const inner = fs.firstChild;
    const w = ref.scrollWidth + 'px';
    if (inner.style.width !== w) inner.style.width = w;
    fs.style.display = 'block';
    // Align to the columns' current offset the moment it appears (without fighting a drag).
    if (wasHidden && fs.scrollLeft !== ref.scrollLeft) {
        _hsSyncing = true;
        fs.scrollLeft = ref.scrollLeft;
        window.requestAnimationFrame(function () { _hsSyncing = false; });
    }
}
window.updateFloatingHScroll = updateFloatingHScroll;

function _scheduleFloatingReposition() {
    if (_hsRepositionRaf) return;
    _hsRepositionRaf = window.requestAnimationFrame(function () {
        _hsRepositionRaf = null; updateFloatingHScroll();
    });
}

// A plot column scrolling horizontally (native scrollbar or syncScrollbars) → keep the
// floating bar's thumb in step. Scroll events don't bubble, so listen in the capture phase.
document.addEventListener('scroll', function (e) {
    const t = e.target;
    if (!(t && _hsMembers.indexOf(t) !== -1)) return;
    if (_hsSyncing) return;
    const fs = document.getElementById('floatingHScroll');
    if (!fs || fs.style.display === 'none') return;
    if (fs.scrollLeft !== t.scrollLeft) {
        _hsSyncing = true;
        fs.scrollLeft = t.scrollLeft;
        window.requestAnimationFrame(function () { _hsSyncing = false; });
    }
}, true);
// Capture phase so this fires for page scroll AND any inner scroll container.
window.addEventListener('scroll', _scheduleFloatingReposition, true);   // -> show/hide + reposition
window.addEventListener('resize', _scheduleFloatingReposition);
window.addEventListener('load', function () {
    [300, 1000, 2000].forEach(function (t) { setTimeout(updateFloatingHScroll, t); });
});

// --- Along the structural alignment (selection overview's right column) ----------------
// The sequence logo and the protein similarity heatmap: computed lazily the first time they
// scroll into view (the selection overview card open), the heatmap again on any control change.
var _analysisLoaded = false;

function _initAnalysisAutoload() {
    var el = document.getElementById('analysis-panel');
    if (!el) return;
    if ('IntersectionObserver' in window) {
        var obs = new IntersectionObserver(function (entries) {
            entries.forEach(function (e) {
                if (e.isIntersecting && !_analysisLoaded) { obs.disconnect(); loadLogo(); loadAnalysis(); }
            });
        }, { rootMargin: '250px' });
        obs.observe(el);
    } else {
        loadLogo(); loadAnalysis();
    }
}
window.addEventListener('load', function () { setTimeout(_initAnalysisAutoload, 400); });

function onAnalysisPropertyChange() {
    var prop = document.getElementById('analysisProperty');
    var flexField = document.getElementById('analysisFlexField');
    if (prop && flexField) flexField.style.display = (prop.value === 'flexibility') ? 'flex' : 'none';
    loadAnalysis();
}
window.onAnalysisPropertyChange = onAnalysisPropertyChange;

function loadAnalysis() {
    if (!document.getElementById('analysis-panel')) return;
    var ids = (window.protIdList || []).join(',');
    var region = (document.getElementById('analysisRegion') || {}).value || 'loops';
    var prop = (document.getElementById('analysisProperty') || {}).value || 'charge';
    var flex = (document.getElementById('analysisFlexMethod') || {}).value || 'cabsflex';
    var exposedEl = document.getElementById('analysisExposed');
    var exposed = exposedEl && exposedEl.checked ? '1' : '0';
    var status = document.getElementById('analysisStatus');
    if (status) status.textContent = 'Computing…';
    _analysisLoaded = true;
    if (window.beginPlotLoad) beginPlotLoad();
    var url = '/catalogue/compare/analysis/?ids=' + encodeURIComponent(ids)
        + '&region=' + encodeURIComponent(region)
        + '&property=' + encodeURIComponent(prop)
        + '&flex=' + encodeURIComponent(flex)
        + '&exposed=' + exposed;
    fetch(url)
        .then(function (r) { return r.json(); })
        .then(function (d) {
            if (!d || !d.ok) {
                if (status) status.textContent = (d && d.error) ? d.error : 'Analysis failed.';
                return;
            }
            renderSimilarity(d);
            if (status) {
                var note = d.n_cols + ' aligned position' + (d.n_cols === 1 ? '' : 's')
                    + ' in ' + (d.region_label || region).toLowerCase() + '.';
                if (d.n_cols < 15) note += '  ⚠ Few positions here — similarity is noisy.';
                status.textContent = note;
            }
        })
        .catch(function () { if (status) status.textContent = 'Analysis request failed.'; })
        .finally(function () { if (window.endPlotLoad) endPlotLoad(); });
}
window.loadAnalysis = loadAnalysis;

// The sequence logo is fixed to the whole catalytic domain — fetched once, independent of the
// Region/Property controls (which drive only the similarity heatmap).
var _logoLoaded = false;
function loadLogo() {
    if (_logoLoaded || !document.getElementById('analysisLogo')) return;
    var ids = (window.protIdList || []).join(',');
    fetch('/catalogue/compare/analysis/?ids=' + encodeURIComponent(ids)
          + '&region=whole&property=sequence&flex=cabsflex&exposed=0&min_occ=0')
        .then(function (r) { return r.json(); })
        .then(function (d) { if (d && d.ok && window.renderLogo) { renderLogo(d); _logoLoaded = true; } })
        .catch(function () {});
}
window.loadLogo = loadLogo;

// Protein x protein similarity heatmap. Rows and columns follow the selection tree on the left
// ("As the tree", the default: a block along the diagonal = the property follows sequence
// similarity) or the clustering by this property; hovering a cell highlights the pair there.
var _lastSimilarity = null;
function _treeOrder(d) {
    var tree = window.selectionTreeInstance;
    if (!tree || !tree.leaves) return null;
    var index = {};
    d.codes.forEach(function (c, i) { index[c] = i; });
    var order = tree.leaves().filter(function (c) { return c in index; })
                             .map(function (c) { return index[c]; });
    var placed = {};
    order.forEach(function (i) { placed[i] = true; });
    d.codes.forEach(function (c, i) { if (!placed[i]) order.push(i); });   // not in the tree
    return order;
}
// The tooltip of one cell: the pair, the value and what it rests on
function _similarityTip(d, i, j) {
    var codes = d.codes;
    if (i === j) return '<b>' + codes[i] + '</b>';
    var head = '<b>' + codes[i] + '  ↔  ' + codes[j] + '</b><br>';
    var v = d.similarity[i][j];
    var n = d.positions ? d.positions[i][j] : null;
    if (v === null || v === undefined) {
        return head + 'too few positions to compare' + (n != null ? ' (' + n + ', at least 5 needed)' : '');
    }
    if (d.categorical) {
        var same = d.matching ? d.matching[i][j] : null;
        var what = d.property === 'aa_groups' ? 'in the same amino-acid group' : 'identical';
        return head + Math.round(v * 100) + '% ' + what
            + (same != null && n != null ? '<br>' + same + ' of ' + n + ' positions' : '');
    }
    var dist = d.distance ? d.distance[i][j] : null;
    return head + 'similarity ' + v.toFixed(2)
        + (dist != null ? '<br>profile difference ' + dist.toFixed(2) + ' (RMS of z-scores)' : '')
        + (n != null ? '<br>over ' + n + ' positions' : '')
        + '<br><span style="color:#888">1 = the same profile, 0 = the most different pair here</span>';
}
// "How is it calculated?" under the controls, for the property shown
function _showSimilarityHow(d) {
    var box = document.getElementById('analysisHow');
    if (!box) return;
    var positions = 'The positions are the aligned columns of the chosen region that at least 30% of the '
        + 'compared proteins occupy (' + d.n_cols + ' here). ';
    var text;
    if (d.categorical) {
        var same = d.property === 'aa_groups'
            ? 'have residues of the same amino-acid group (hydrophobic, aromatic, polar, positive, negative, special)'
            : 'have the same residue';
        text = positions + 'Similarity = the share of these positions where both proteins ' + same
            + '. A residue facing a gap counts as different; a position where both have a gap is left out. '
            + 'So, unlike the tree (BLOSUM62, which credits similar residues and skips gaps), insertions and '
            + 'loop-length differences lower the similarity.';
    } else {
        text = positions + 'At each position the values are standardised across the compared proteins '
            + '(z-scores), so every position weighs the same. Two proteins differ by the root-mean-square of '
            + 'their z-differences over the positions where both have a residue (gaps are left out); '
            + 'similarity = 1 − that difference ÷ the largest difference in this comparison. 1 = the same '
            + 'profile, 0 = the most different pair here: the scale is relative to the compared proteins.';
    }
    box.textContent = text + ' A pair needs at least 5 shared positions.';
}
function renderSimilarity(d) {
    var el = document.getElementById('analysisSimilarity');
    if (!el || !window.Plotly) return;
    _lastSimilarity = d;
    var codes = d.codes;
    var wantTree = ((document.getElementById('analysisOrder') || {}).value || 'tree') === 'tree';
    var order = (wantTree && _treeOrder(d)) || d.order;
    var labels = order.map(function (i) { return codes[i]; });
    var z = [], text = [];
    order.forEach(function (i) {
        var zrow = [], trow = [];
        order.forEach(function (j) {
            var v = d.similarity[i][j];
            zrow.push(v);
            trow.push(_similarityTip(d, i, j));
        });
        z.push(zrow); text.push(trow);
    });
    _showSimilarityHow(d);
    var n = labels.length;
    var room = Math.max(320, (el.clientWidth || 700) - 150);   // the colour bar and labels
    var px = Math.max(320, Math.min(room, 760, 70 + n * 26));
    var data = [{
        type: 'heatmap', z: z, x: labels, y: labels,
        text: text, hovertemplate: '%{text}<extra></extra>',
        colorscale: 'YlGnBu', reversescale: true, zmin: 0, zmax: 1,
        colorbar: { title: { text: 'similarity', side: 'right' }, thickness: 10, len: 0.85 }
    }];
    var layout = {
        title: { text: d.property_label + ' · ' + d.region_label
                       + (order === d.order ? ' · clustered' : ' · tree order'),
                 font: { size: 12 } },
        width: px + 150, height: px,
        hovermode: 'closest',
        margin: { l: 95, r: 20, t: 36, b: 95 },
        xaxis: { tickfont: { size: 9 }, tickangle: -60, automargin: true, constrain: 'domain' },
        yaxis: { tickfont: { size: 9 }, automargin: true, autorange: 'reversed',
                 scaleanchor: 'x', constrain: 'domain' },
        paper_bgcolor: 'white'
    };
    Plotly.react(el, data, layout, { displayModeBar: false, responsive: false });
    if (!el._treeHoverBound && el.on) {
        el._treeHoverBound = true;
        el.on('plotly_hover', function (ev) {
            var p = ev && ev.points && ev.points[0];
            if (p) _logoHighlightTree(p.x === p.y ? [p.x] : [p.x, p.y]);
        });
        el.on('plotly_unhover', function () { _logoHighlightTree([]); });
    }
}
window.renderSimilarity = renderSimilarity;
function rerenderSimilarity() { if (_lastSimilarity) renderSimilarity(_lastSimilarity); }
window.rerenderSimilarity = rerenderSimilarity;
// The tree is drawn when the selection overview opens, maybe after the heatmap: redraw in its order.
document.addEventListener('selection-tree-ready', rerenderSimilarity);

// The sequence composition logo (letter height = the share of the proteins with that residue).
var AA_LOGO_COLOR = {   // WebLogo "chemistry" scheme
    G: '#159a4a', S: '#159a4a', T: '#159a4a', Y: '#159a4a', C: '#159a4a',
    N: '#8a2be2', Q: '#8a2be2',
    K: '#2166d0', R: '#2166d0', H: '#2166d0',
    D: '#d0342b', E: '#d0342b',
    A: '#222', V: '#222', L: '#222', I: '#222', P: '#222', W: '#222', F: '#222', M: '#222'
};
function _svgEl(tag, attrs) {
    var e = document.createElementNS('http://www.w3.org/2000/svg', tag);
    for (var k in attrs) e.setAttribute(k, attrs[k]);
    return e;
}
// Rich hover tooltip for a logo column: full residue composition (incl. gap %). Reuses the
// body-level #plotHoverTip element the plots use, so it can't be clipped by the plot rows.
function _logoTipEl() {
    var el = document.getElementById('plotHoverTip');
    if (!el) { el = document.createElement('div'); el.id = 'plotHoverTip'; document.body.appendChild(el); }
    return el;
}
function _logoColHtml(d, ci, hoveredAa) {
    var comp = (d.comp && d.comp[ci]) || [];
    var members = (d.members && d.members[ci]) || {};
    // The proteins are numbered differently, so a column has no single residue number: the
    // heading is the column, and each protein carries its own number below.
    var head = 'column ' + (ci + 1) + ' of ' + d.n_cols
             + (d.occ && d.occ[ci] != null ? ' · ' + Math.round(d.occ[ci] * 100) + '% present' : '');
    // Hovering a letter marks that row and greys the others, so the tooltip says which
    // residue the tree highlight belongs to while the whole composition stays readable.
    var rows = comp.map(function (s) {
        var sym = s[0], isGap = sym === 'gap';
        var on = !hoveredAa || sym === hoveredAa;
        var color = isGap ? '#999' : (AA_LOGO_COLOR[sym] || '#333');
        var who = members[sym] || [];
        return '<div style="line-height:1.35;margin-bottom:2px;' + (on ? '' : 'opacity:0.45;') + '">'
             + '<div style="display:flex;justify-content:space-between;gap:14px;">'
             + '<span style="color:' + color + ';font-weight:bold;font-family:monospace;">'
             + (isGap ? 'gap' : sym) + (on && hoveredAa ? ' ◄' : '') + '</span>'
             + '<span>' + s[1] + '%</span></div>'
             + (who.length ? '<div style="font-size:11px;color:#5a6b62;padding-left:2px;max-width:300px;">'
                             + who.map(_memberLabel).join(', ') + '</div>' : '')
             + '</div>';
    }).join('');
    var lanes = _logoLaneHtml(d, ci);
    return '<div style="font-weight:600;margin-bottom:3px;">' + head + '</div>' + rows + lanes;
}
// A member entry is [code, residue number]; the number is the one THAT protein uses, and is
// absent where it has a gap at this column.
function _memberCode(m) { return Array.isArray(m) ? m[0] : m; }
function _memberLabel(m) {
    if (!Array.isArray(m)) return m;
    return m[1] == null ? m[0]
        : m[0] + ' <span style="color:#8a9a92;">' + m[1] + '</span>';
}
// The regions this column falls in, as the lanes under the logo show them.
function _logoLaneHtml(d, ci) {
    if (!d.lanes) return '';
    var parts = [];
    (d.lane_names.site || []).forEach(function (name) {
        var m = (d.lanes['site:' + name] || [])[ci] || [];
        if (m.length) parts.push(name + ' (' + m.length + ')');
    });
    (d.lane_names.loops || []).forEach(function (name) {
        var m = (d.lanes['loop:' + name] || [])[ci] || [];
        if (m.length) parts.push('<b style="color:#' + (LOGO_LOOP_COLOR[name] || '888') + '">'
                                 + name + '</b> ' + m.length + '/' + d.codes.length);
    });
    var bind = (d.lanes.binding || [])[ci] || [];
    if (bind.length) parts.push('binding surface ' + bind.length + '/' + d.codes.length);
    if (!parts.length) return '';
    return '<div style="margin-top:4px;padding-top:4px;border-top:1px solid #e3ebe6;'
         + 'font-size:11px;color:#3d5a4c;">' + parts.join(' · ') + '</div>';
}
// Highlight a set of protein codes in the selection tree, when that card has drawn one.
function _logoHighlightTree(codes) {
    var tree = window.selectionTreeInstance;
    if (tree && tree.setOverviewHighlight) tree.setOverviewHighlight(codes || []);
}
function _buildLogoSVG(d) {
    var n = d.n_cols, cols = d.logo || [];
    var COL_W = 16, LOGO_H = 108, PAD_L = 46, PAD_T = 6, PAD_B = 4;
    var MAX = 1.0;                                 // frequency logo: y-axis is 0..100%
    var pph = LOGO_H / MAX;                        // pixels per unit fraction
    var NAT_W = 60, NAT_CAP = 72;                 // monospace-bold glyph metrics at font-size 100
    var width = PAD_L + n * COL_W + 6, height = PAD_T + LOGO_H + PAD_B;
    var svg = _svgEl('svg', { width: width, height: height, viewBox: '0 0 ' + width + ' ' + height });
    var y0 = PAD_T + LOGO_H;                       // baseline (0%)

    [0, 0.25, 0.5, 0.75, 1.0].forEach(function (f) {   // y grid + % ticks
        var y = y0 - f * pph;
        svg.appendChild(_svgEl('rect', { x: PAD_L, y: y, width: n * COL_W, height: 1, fill: f === 0 ? '#bbb' : '#eee' }));
        var t = _svgEl('text', { x: PAD_L - 4, y: y + 3, 'text-anchor': 'end', 'font-size': 8, fill: '#888' });
        t.textContent = Math.round(f * 100); svg.appendChild(t);
    });
    var yl = _svgEl('text', { transform: 'translate(8,' + (PAD_T + LOGO_H / 2) + ') rotate(-90)', 'text-anchor': 'middle', 'font-size': 8, fill: '#666' });
    yl.textContent = '% of proteins'; svg.appendChild(yl);

    // Position cursor: a pink highlight box that follows the hovered column (behind the letters).
    var cursor = _svgEl('rect', { x: 0, y: PAD_T, width: COL_W, height: LOGO_H,
                                  fill: 'rgba(255,20,147,0.15)', stroke: 'rgba(255,20,147,0.75)',
                                  'stroke-width': 1, display: 'none', 'pointer-events': 'none' });
    svg.appendChild(cursor);

    // Letter bands are collected while drawing and turned into their own hit areas below, so
    // hovering one letter names only the proteins carrying THAT residue.
    var letterBands = [];
    for (var i = 0; i < n; i++) {                  // letters, stacked most-frequent on top
        var stack = cols[i] || [], cellX = PAD_L + i * COL_W, cum = 0;
        for (var k = stack.length - 1; k >= 0; k--) {
            var aa = stack[k][0], h = stack[k][1] * pph;
            if (h < 0.7) { cum += h; continue; }
            var cellTopY = y0 - cum - h;
            var g = _svgEl('g', { transform: 'translate(' + cellX + ',' + (cellTopY + h) + ') scale(' + (COL_W / NAT_W) + ',' + (h / NAT_CAP) + ')' });
            var t = _svgEl('text', { x: NAT_W / 2, y: 0, 'text-anchor': 'middle', 'font-size': 100, 'font-weight': 'bold', 'font-family': 'monospace', fill: AA_LOGO_COLOR[aa] || '#333' });
            t.textContent = aa; g.appendChild(t); svg.appendChild(g);
            letterBands.push({ col: i, aa: aa, x: cellX, y: cellTopY, h: h });
            cum += h;
        }
        var hr = _svgEl('rect', { x: cellX, y: PAD_T, width: COL_W, height: LOGO_H,
                                  fill: '#000', opacity: 0, 'data-col': i });
        hr.style.cursor = 'default';
        hr.addEventListener('mousemove', _logoTipMove);
        hr.addEventListener('mouseleave', _logoTipHide);
        svg.appendChild(hr);
    }
    // After the column rects, so a letter takes the pointer and the rest of the column (the
    // gap space above the stack, and letters too short to draw) still falls back to it.
    letterBands.forEach(function (b) {
        var lr = _svgEl('rect', { x: b.x, y: b.y, width: COL_W, height: b.h,
                                  fill: '#000', opacity: 0, 'data-col': b.col, 'data-aa': b.aa });
        lr.style.cursor = 'pointer';
        lr.addEventListener('mousemove', _logoTipMove);
        lr.addEventListener('mouseleave', _logoTipHide);
        svg.appendChild(lr);
    });
    // Handlers read the column index off the rect (data-col) and close over `d`.
    function _logoTipMove(ev) {
        var ci = +this.getAttribute('data-col');
        var aa = this.getAttribute('data-aa');     // set on a letter's own hit area
        var members = (d.members && d.members[ci]) || {};
        // On a letter: that residue's proteins only. Elsewhere in the column: all present.
        // Two highlights: a letter marks only its own band (cursor, behind the letters); the
        // gap space above a stack or any lane cell marks the whole position, one tall outline
        // from the top of the logo down through the last lane (columnCursor, drawn on top).
        var lit, colX = PAD_L + ci * COL_W;
        if (aa && members[aa]) {
            lit = members[aa].map(_memberCode);
            cursor.setAttribute('x', colX);
            cursor.setAttribute('y', this.getAttribute('y'));
            cursor.setAttribute('height', this.getAttribute('height'));
            cursor.setAttribute('display', '');
            columnCursor.setAttribute('display', 'none');
        } else {
            lit = [];
            Object.keys(members).forEach(function (sym) {
                if (sym !== 'gap') lit = lit.concat(members[sym].map(_memberCode));
            });
            cursor.setAttribute('display', 'none');
            columnCursor.setAttribute('x', colX);
            columnCursor.setAttribute('height', +svg.getAttribute('height') - PAD_T - 2);
            columnCursor.setAttribute('display', '');
        }
        _logoHighlightTree(lit);
        var tip = _logoTipEl();
        tip.innerHTML = _logoColHtml(d, ci, aa);
        tip.style.display = 'block';
        tip.style.background = '#ffffff';
        tip.style.color = '#243b30';
        tip.style.border = '1px solid #cddbd3';
        tip.style.boxShadow = '0 4px 14px rgba(0,0,0,0.16)';
        // Measure the (now-rendered) box, then follow the cursor horizontally but flip to the
        // left near the right edge, and anchor y to the hovered column's row so it doesn't
        // jitter vertically — same behaviour as the plot tooltips above.
        var colRect = this.getBoundingClientRect();
        var r = tip.getBoundingClientRect();
        var pad = 10;
        var x = ev.clientX + 14;
        if (x + r.width + pad > window.innerWidth) x = ev.clientX - r.width - 14;
        if (x < pad) x = pad;
        var y = colRect.top + colRect.height / 2 - r.height / 2;
        if (y < pad) y = pad;
        if (y + r.height + pad > window.innerHeight) y = window.innerHeight - r.height - pad;
        tip.style.left = x + 'px';
        tip.style.top = y + 'px';
    }
    function _logoTipHide() {
        cursor.setAttribute('display', 'none');
        columnCursor.setAttribute('display', 'none');
        _logoHighlightTree([]);
        var t = document.getElementById('plotHoverTip');
        if (t) t.style.display = 'none';
    }
    _buildLogoLanes(svg, d, { padLeft: PAD_L, colWidth: COL_W, top: PAD_T + LOGO_H + PAD_B,
                              onHover: _logoTipMove, onLeave: _logoTipHide });
    // The whole-position highlight: appended last, so neither the letters nor the coloured
    // lane cells cover it; its height is set on hover, once the lanes have sized the svg.
    var columnCursor = _svgEl('rect', { x: 0, y: PAD_T, width: COL_W, height: LOGO_H,
                                        fill: 'rgba(255,20,147,0.10)', stroke: 'rgba(255,20,147,0.85)',
                                        'stroke-width': 1.2, display: 'none', 'pointer-events': 'none' });
    svg.appendChild(columnCursor);
    return svg;
}

// Lanes under the logo: how much of the SELECTION shares each region at a column, so the
// conservation of the regions reads against the conservation of the residues above them.
// Darkness = the fraction of proteins; the hover is the logo column's, so the tooltip and the
// tree highlight behave the same whether the pointer is on a letter or a lane.
var LANE_H = 9, LANE_GAP = 2;
// The Sequence plot's colours for the same categories (makePlots: the histidine brace magenta,
// the axial position dark orange, SEQ_MOTIF_COLOR teal, SEQ_DISTAL_COLOR brown) and
// BINDING_SURFACE_COLOR for the binding surface, so a category reads the same in both plots.
var LANE_COLOR = { binding: '00796B' };
var SITE_LANE_COLOR = {
    'N-His': 'BA03BA', 'Internal-His': 'BA03BA', 'Axial aromatic': 'FF8C00',
    'Second sphere': '26A69A', 'Shielding': '8D6E63',
};
// The loop colours of the Secondary Structure and Loops plot (makePlots.MOLSTAR_COLOR_BY_GROUP),
// so a loop is the same colour here, in that plot and in the 3D viewer.
var LOGO_LOOP_COLOR = { L2: 'FF4500', L3: 'FF8C00', LS: '8A2BE2', L8: '1E90FF', LC: '228B22' };
function _laneRows(d) {
    if (!d.lanes || !d.lane_names) return [];
    var rows = [];
    var SITE_LABEL = { 'N-His': 'N-His', 'Internal-His': 'Int-His', 'Axial aromatic': 'axial',
                       'Second sphere': '2nd sph', 'Shielding': 'shield' };
    (d.lane_names.site || []).forEach(function (n) {
        if ((d.lanes['site:' + n] || []).some(function (m) { return m.length; }))
            rows.push({ key: 'site:' + n, label: SITE_LABEL[n] || n, color: SITE_LANE_COLOR[n] || '8E44AD' });
    });
    rows.push({ key: 'binding', label: 'bind', color: LANE_COLOR.binding });
    // One lane for all loops: the colour says which loop, the opacity how much of the
    // selection has it there. Overlapping assignments take the loop most proteins agree on.
    var loops = (d.lane_names.loops || []).filter(function (n) {
        return (d.lanes['loop:' + n] || []).some(function (m) { return m.length; });
    });
    if (loops.length) rows.push({ key: 'loops', label: 'loops', loops: loops });
    return rows.filter(function (r) {
        return r.loops || (d.lanes[r.key] || []).some(function (m) { return m.length; });
    });
}
// Hovering a lane's name marks every column that lane occupies — through the letters as well
// as the lanes — bolds the name, and highlights the proteins that have it in the selection
// tree, so the region can be read off the logo without hunting for its cells.
function _bindLaneLabel(label, svg, d, row, marks, opt) {
    function columnsOf() {
        var out = [], codes = {};
        for (var i = 0; i < d.n_cols; i++) {
            var members = row.loops
                ? (row.loops.reduce(function (acc, name) {
                      return acc.concat((d.lanes['loop:' + name] || [])[i] || []); }, []))
                : ((d.lanes[row.key] || [])[i] || []);
            if (members.length) {
                out.push(i);
                members.forEach(function (m) { codes[_memberCode(m)] = true; });
            }
        }
        return { cols: out, codes: Object.keys(codes) };
    }
    label.addEventListener('mouseenter', function () {
        var hit = columnsOf();
        marks.replaceChildren();
        hit.cols.forEach(function (i) {
            marks.appendChild(_svgEl('rect', {
                x: opt.padLeft + i * opt.colWidth, y: 2, width: opt.colWidth,
                height: (+svg.getAttribute('height')) - 4,
                fill: '#' + (row.color || 'B8860B'), 'fill-opacity': 0.16 }));
        });
        label.setAttribute('font-weight', 'bold');
        label.setAttribute('fill', '#' + (row.color || 'B8860B'));
        _logoHighlightTree(hit.codes);
    });
    label.addEventListener('mouseleave', function () {
        marks.replaceChildren();
        label.removeAttribute('font-weight');
        label.setAttribute('fill', '#666');
        _logoHighlightTree([]);
    });
}
// The loop most of the selection puts at this column, with its share: [name, share] or null.
function _dominantLoop(d, row, ci) {
    var best = null, total = d.codes.length || 1;
    row.loops.forEach(function (name) {
        var k = ((d.lanes['loop:' + name] || [])[ci] || []).length;
        if (k && (!best || k > best[1])) best = [name, k];
    });
    return best ? [best[0], best[1] / total] : null;
}
function _buildLogoLanes(svg, d, opt) {
    var rows = _laneRows(d);
    if (!rows.length) return 0;
    var n = d.n_cols, total = d.codes.length, y = opt.top + 4;
    // Bands drawn behind everything when a lane's name is hovered: they mark, through the
    // letters and the other lanes, every column that lane occupies.
    var marks = _svgEl('g', { 'pointer-events': 'none' });
    svg.insertBefore(marks, svg.firstChild);
    rows.forEach(function (row) {
        var label = _svgEl('text', { x: opt.padLeft - 4, y: y + LANE_H - 1, 'text-anchor': 'end',
                                     'font-size': 8, fill: '#666' });
        label.textContent = row.label;
        label.style.cursor = 'pointer';
        _bindLaneLabel(label, svg, d, row, marks, opt);
        svg.appendChild(label);
        var cells = d.lanes[row.key] || [];
        for (var i = 0; i < n; i++) {
            var share, color;
            if (row.loops) {
                var top = _dominantLoop(d, row, i);
                share = top ? top[1] : 0;
                color = top ? (LOGO_LOOP_COLOR[top[0]] || '888888') : null;
            } else {
                share = (cells[i] || []).length / (total || 1);
                color = row.color;
            }
            if (share > 0) {
                svg.appendChild(_svgEl('rect', {
                    x: opt.padLeft + i * opt.colWidth, y: y, width: opt.colWidth, height: LANE_H,
                    fill: '#' + color, 'fill-opacity': (0.18 + 0.82 * share).toFixed(3) }));
            }
            var hit = _svgEl('rect', { x: opt.padLeft + i * opt.colWidth, y: y,
                                       width: opt.colWidth, height: LANE_H,
                                       fill: '#000', opacity: 0, 'data-col': i });
            hit.addEventListener('mousemove', opt.onHover);
            hit.addEventListener('mouseleave', opt.onLeave);
            svg.appendChild(hit);
        }
        y += LANE_H + LANE_GAP;
    });
    var height = y + 2;
    svg.setAttribute('height', height);
    svg.setAttribute('viewBox', '0 0 ' + svg.getAttribute('width') + ' ' + height);
    return height;
}
// The sequence logo is a FIXED overview of the whole catalytic domain (loaded once, its own
// fetch), so it doesn't change when Region/Property are changed — those drive only the heatmap.
function renderLogo(d) {
    var el = document.getElementById('analysisLogo');
    if (!el) return;
    if (!d || !d.n_cols) { el.innerHTML = ''; return; }
    el.innerHTML = '';
    var wrap = document.createElement('div');
    wrap.className = 'logo-scroll';
    wrap.appendChild(_buildLogoSVG(d));
    el.appendChild(wrap);
    var cap = document.createElement('div');
    cap.className = 'analysis-caption';
    cap.innerHTML = 'Lanes underneath: how many of the compared proteins have each region at that '
        + 'position (darker = more), in the Sequence plot\'s colours: '
        + '<b style="color:#BA03BA">brace histidines</b>, <b style="color:#FF8C00">axial aromatic</b>, '
        + '<b style="color:#26A69A">second-sphere motif</b>, <b style="color:#8D6E63">shielding residue</b>, '
        + '<b style="color:#00796B">binding surface</b>, and the loops '
        + '<b style="color:#FF4500">L2</b>, <b style="color:#FF8C00">L3</b>, '
        + '<b style="color:#8A2BE2">LS</b>, <b style="color:#1E90FF">L8</b>, '
        + '<b style="color:#228B22">LC</b>.';
    el.appendChild(cap);
}
window.renderLogo = renderLogo;

// Fetch one plot-type's plots for the current protein set and render them into
// `scrollable`. The data comes from the server plot-bundle cache, so this is
// cheap once the bundle has been built (by the initial page render).
function loadSectionPlots(section, scrollable, plotType, scale) {
    const ids = (window.protIdList || []).join(',');
    // If the Hydrophobicity section is opened after the user already picked a
    // scale, honour that choice instead of loading the default.
    if (!scale && plotType === 'Hydrophobicity' && window.currentHydrophobScale) {
        scale = window.currentHydrophobScale;
    }
    // Flexibility reuses the `scale` slot to carry the chosen method (cabsflex/gnm/…).
    if (!scale && plotType === 'Flexibility' && window.currentFlexMethod) {
        scale = window.currentFlexMethod;
    }
    let url = `/catalogue/compare/plots/?ptype=${encodeURIComponent(plotType)}&ids=${encodeURIComponent(ids)}`;
    if (scale) url += `&scale=${encodeURIComponent(scale)}`;
    // Plots default to the catalytic domain (like the 3D view); the "Full length"
    // toggle switches every plot to the whole length.
    if (document.getElementById('showFullStructure') && document.getElementById('showFullStructure').checked) {
        url += '&full=1';
    }
    ensureLoadingSpinnerCss();
    scrollable.innerHTML = '<div class="lpmo-loading"><span class="lpmo-spinner"></span>Loading ' + plotType + ' plots…</div>';

    beginPlotLoad();   // lock the page while this section loads (batched across sections)
    fetch(url, { method: 'GET' })
        .then(r => r.json())
        .then(data => {
            if (data.error) throw new Error(data.error);
            // closed (and freed) while its data was on the way: don't draw it
            if (window.ONE_PLOT_SECTION && section.style.display === 'none') {
                freePlotSection(section);
                return;
            }
            return renderPlotsIntoColumn(scrollable, data.plots_html || [], plotType);
        })
        .then(function () {
            // Wire hover/click → 3D viewer for the newly added plots (idempotent).
            try { hoverSequenceSelection((window.protIdList || []).length); } catch (e) {}
        })
        .catch(error => {
            console.error('Error loading plots:', error);
            section.setAttribute('data-loaded', '0');   // allow a retry on next open
            scrollable.innerHTML = '<div style="padding:8px;color:#d32f2f;">Error loading plots. Please refresh the page.</div>';
        })
        .finally(function () { endPlotLoad(); });
}

// Current shared horizontal scroll offset of the open plot sections (all synced by
// syncScrollbars). Collapsed/empty columns report 0, so take the max across the others.
function currentPlotsScrollLeft(exclude) {
    let sl = 0;
    document.querySelectorAll('.scrollable-plots-column').forEach(function (col) {
        if (col !== exclude && col.scrollLeft > sl) sl = col.scrollLeft;
    });
    return sl;
}

// Protein visibility as a per-INDEX array (from the checkbox `value` = protein index), NOT
// DOM order — so sorting the protein table doesn't scramble which proteins Mol*/plots show.
function getProteinVisibility(protCount) {
    const n = protCount || (window.protIdList || []).length;
    const vis = new Array(n).fill(true);
    document.querySelectorAll("input.proteinVisibilityCheckbox").forEach(function (cb) {
        const idx = parseInt(cb.value, 10);
        if (!isNaN(idx) && idx >= 0 && idx < n) vis[idx] = cb.checked;
    });
    return vis;
}
window.getProteinVisibility = getProteinVisibility;

// Re-apply protein visibility to the plot rows: hidden proteins get display:none on their
// plot row + y-axis label. Keyed by the stable `data-protein-index` (NOT the
// protein-plot-container-N class, which reorderPlots rewrites to the visual position on sort),
// so it stays correct after sorting. `:not(.proteinEntry)` excludes the table rows, which also
// carry data-protein-index. Needed because re-rendered plot DOM starts fully visible.
// Visible rows get their CSS display back (the labels are flex boxes centring their text).
function applyPlotVisibility() {
    document.querySelectorAll("input.proteinVisibilityCheckbox").forEach(function (cb) {
        const idx = parseInt(cb.value, 10);
        if (isNaN(idx)) return;
        const disp = cb.checked ? '' : 'none';
        document.querySelectorAll('[data-protein-index="' + idx + '"]:not(.proteinEntry)').forEach(function (c) {
            c.style.display = disp;
        });
    });
}
window.applyPlotVisibility = applyPlotVisibility;

// Render an array of Plotly-JSON strings into `scrollable` as .plotWrapper divs,
// one per protein (in order). Returns a Promise that resolves once every plot is
// drawn. Shared by lazy section loading and the hydrophobicity scale switcher.
function renderPlotsIntoColumn(scrollable, plotsHtml, plotType) {
    // Match the horizontal scroll of the sections already open, so a freshly loaded
    // section lines up with them instead of jumping back to the alignment start. Capture
    // before clearing: own position (scale switch on the only open column) or the siblings'.
    const targetScroll = Math.max(scrollable.scrollLeft || 0, currentPlotsScrollLeft(scrollable));
    scrollable.innerHTML = '';
    const pairs = [];
    plotsHtml.forEach(function (h, idx) {
        const w = document.createElement('div');
        w.className = 'plotWrapper protein-plot-container-' + idx;
        w.setAttribute('data-protein-index', idx);
        w.setAttribute('data-plot-type', plotType);
        scrollable.appendChild(w);
        pairs.push([w, h]);
    });
    const promises = pairs.map(function (pair) {
        let parsed;
        try { parsed = JSON.parse(pair[1]); } catch (e) { return Promise.resolve(); }
        return Plotly.newPlot(pair[0], parsed.data, parsed.layout,
            { displayModeBar: false, responsive: false, scrollZoom: false, doubleClick: false });
    });
    return Promise.all(promises).then(function () {
        // Now that the plots have their full width, snap this column to the shared scroll
        // offset (re-apply on the next frame in case the width isn't flushed yet).
        if (targetScroll > 0) {
            scrollable.scrollLeft = targetScroll;
            window.requestAnimationFrame(function () { scrollable.scrollLeft = targetScroll; });
        }
        // Re-rendered plots come back in the original protein order; if the table was sorted,
        // re-apply that order so plot rows line up with the (persistent, still-sorted) labels.
        if (window.currentSortOrder && window.reorderPlots) {
            try { window.reorderPlots(window.currentSortOrder); } catch (e) {}
        }
        // Freshly-rendered rows default to visible; re-hide de-selected proteins (else the
        // Full-length toggle / a section reload re-shows proteins the user had hidden).
        if (window.applyPlotVisibility) window.applyPlotVisibility();
        scheduleLabelAlign();
        [200, 600, 1500].forEach(function (t) { setTimeout(scheduleLabelAlign, t); });
        setTimeout(function () { if (window.initializeVerticalCursor) window.initializeVerticalCursor(); }, 100);
        if (window.updateFloatingHScroll) {
            updateFloatingHScroll();
            [250, 800].forEach(function (t) { setTimeout(updateFloatingHScroll, t); });
        }
        if (window.decorateRegionLabels) {
            if (plotType === 'Hole Hopping') {
                [50, 400].forEach(function (t) { setTimeout(function () {
                    decorateRegionLabels('Hole Hopping', 'holeHopping', { r: 233, g: 30, b: 140 }, 'hole-hopping residues');
                }, t); });
            } else if (plotType === SS_LOOPS_PLOT_TYPE) {
                [50, 400].forEach(function (t) { setTimeout(function () {
                    decorateRegionLabels(SS_LOOPS_PLOT_TYPE, 'loopRegions', { r: 255, g: 20, b: 147 }, 'loops');
                }, t); });
            }
        }
    });
}

function setupProteinVisibilityDropdown(protCount){
    const dropdown = document.getElementById('proteinVisibilityDropdown');
    const toggleBtn = document.getElementById('proteinVisibilityToggle');
    const panel = document.getElementById('proteinVisibilityPanel');
    const selectAllBtn = document.getElementById('selectAllProteinsBtn');
    const clearAllBtn = document.getElementById('clearAllProteinsBtn');
    const dropdownCheckboxes = panel ? panel.querySelectorAll('.dropdownProteinCheckbox') : [];
    const mainBoxes = document.querySelectorAll('input.proteinVisibilityCheckbox');

    if (!dropdown || !toggleBtn || !panel) return;

    const buildNamesLabel = (names) => {
        if (!names.length) {
            return 'No proteins selected';
        }

        const total = dropdownCheckboxes.length;
        if (names.length === total) {
            return 'All proteins';
        }

        let result = '';
        let count = 0;
        for (const name of names) {
            const delimiter = result.length ? ', ' : '';
            if ((result + delimiter + name).length > 32) {
                result += '…';
                return result;
            }
            result += delimiter + name;
            count++;
        }

        if (count < names.length) {
            result += ', …';
        }
        return result;
    };

    const updateLabel = () => {
        const selectedBoxes = Array.from(dropdownCheckboxes).filter(cb => cb.checked);
        const names = selectedBoxes
            .map(cb => cb.getAttribute('data-name') || '')
            .filter(Boolean);

        const label = buildNamesLabel(names);
        toggleBtn.textContent = label;
        toggleBtn.title = names.length ? names.join(', ') : label;
    };

    const applyFromDropdown = () => {
        dropdownCheckboxes.forEach(cb => {
            const idx = parseInt(cb.getAttribute('data-index'));
            // Match the table checkbox by value (protein index), not DOM position — the table
            // may have been sorted, which would otherwise write to the wrong protein's box.
            const main = document.querySelector(`input.proteinVisibilityCheckbox[value="${idx}"]`);
            if (!isNaN(idx) && main) {
                main.checked = cb.checked;
            }
        });

        try {
            const viewer = window.viewerInstance || viewerInstance;
            const numReps = window.numberOfDifferentRepresentations || numberOfDifferentRepresentations;
            if (viewer && numReps) {
                toggleProteinVisibility(viewer, numReps, protCount);
            }
        } catch (e) {
            console.error('Error in applyFromDropdown:', e);
        }
        updateLabel();
    };

    dropdownCheckboxes.forEach(cb => {
        cb.addEventListener('change', applyFromDropdown);
    });

    if (selectAllBtn) {
        selectAllBtn.addEventListener('click', () => {
            dropdownCheckboxes.forEach(cb => cb.checked = true);
            applyFromDropdown();
        });
    }

    if (clearAllBtn) {
        clearAllBtn.addEventListener('click', () => {
            dropdownCheckboxes.forEach(cb => cb.checked = false);
            applyFromDropdown();
        });
    }

    // "Only: <group>": show just that group's proteins
    panel.querySelectorAll('.groupOnlyBtn').forEach(btn => {
        btn.addEventListener('click', () => {
            const groupIndex = parseInt(btn.getAttribute('data-group-index'), 10);
            const of = comparisonGroups().of;
            dropdownCheckboxes.forEach(cb => {
                cb.checked = of[parseInt(cb.getAttribute('data-index'), 10)] === groupIndex;
            });
            applyFromDropdown();
        });
    });

    mainBoxes.forEach((box, idx) => {
        box.addEventListener('change', () => {
            const target = panel.querySelector(`.dropdownProteinCheckbox[data-index="${idx}"]`);
            if (target) {
                target.checked = box.checked;
                updateLabel();
            }
        });
    });

    toggleBtn.addEventListener('click', (event) => {
        event.stopPropagation();
        panel.classList.toggle('open');
    });

    document.addEventListener('click', (event) => {
        if (!dropdown.contains(event.target)) {
            panel.classList.remove('open');
        }
    });

    updateLabel();
}


async function toggleProteinVisibility(viewerInstance, repCount, protCount){
    // Visibility keyed by protein index (checkbox value), so table sorting can't scramble it.
    const visArr = getProteinVisibility(protCount);

    let visIndexArr = [];

    if (lastAllData.length > 0) {
        for (let i = 0; i < protCount; i++){
            if (visArr[i]){
                viewerInstance.visual.select({data: lastAllData[i], structureNumber: i+1, keepColors: false, keepRepresentations: false});
            } else {
                viewerInstance.visual.select({data: [], structureNumber: i+1, keepColors: false, keepRepresentations: false});
            }
        }
    }

    // (compare_layout.css greys the table rows of hidden proteins)
    for (let i = 0; i < protCount; i++){
        if (visArr[i]) visIndexArr.push(i);
    }

    // Hide/show plot rows by stable data-protein-index (robust to table sorting).
    applyPlotVisibility();

    const plots = alignmentPlotDivs();
    typeArr = groupPlotsByType(plots, protCount);

    for (let i = 0; i < typeArr.length; i++){
        Plotly.relayout(typeArr[i][visIndexArr.length - 1], { 'xaxis.visible': true});
        if (visIndexArr.length - 2 >= 0){
            Plotly.relayout(typeArr[i][visIndexArr.length - 2], { 'xaxis.visible': false});
        }
    }

    hoverSequenceSelection(protCount);
};



// The comparison's groups (json_script "proteinGroups"): {groups: [{name, color, codes,
// ungrouped?}], of: [group index per protein, or null]}. Read once; not window.proteinGroups,
// which the browser maps to the element of that id.
function comparisonGroups(){
    if (!window._comparisonGroups) {
        try { window._comparisonGroups = JSON.parse(document.getElementById('proteinGroups').textContent); }
        catch (e) { window._comparisonGroups = {groups: [], of: []}; }
    }
    return window._comparisonGroups;
}
window.comparisonGroups = comparisonGroups;
// Residue numbers of each protein's aromatic binding residues (json_script "bindingAromaticNumbers")
function bindingAromaticNumbers(){
    if (!window._bindingAromaticNumbers) {
        try { window._bindingAromaticNumbers = JSON.parse(document.getElementById('bindingAromaticNumbers').textContent) || []; }
        catch (e) { window._bindingAromaticNumbers = []; }
    }
    return window._bindingAromaticNumbers;
}
function groupColorOf(proteinIndex){
    const data = comparisonGroups();
    const group = data.groups[data.of[proteinIndex]];
    const hex = group ? group.color : '#a3a3a3';
    return {r: parseInt(hex.slice(1, 3), 16), g: parseInt(hex.slice(3, 5), 16), b: parseInt(hex.slice(5, 7), 16)};
}

async function initializeRepresentations(molstarColoring, instance, coloringOption, protCount){

    const isFirstLoad = (proteinLoadingTime === undefined);

    let molstarContainer = document.getElementById('myViewer');
    let sequenceColoringOptions = Object.keys(molstarColoring['sequence']);
    let allData = [];

    // Read active main representation type from radio buttons
    const mainRepSel = document.getElementById('mainRepSelect');
    const mainRep = mainRepSel ? mainRepSel.value : 'cartoon';

    // Read overlay checkbox state
    const overlayCbs = [...document.querySelectorAll('input.proteinRepresentationCheckbox:checked')];
    const overlayOn = name => overlayCbs.some(cb => cb.dataset.overlay === name);
    const showBrace        = overlayOn('brace');
    const showSecondSphere = overlayOn('secondSphere');
    const showHoleHopping  = overlayOn('holeHopping');
    const showBindingSurface = overlayOn('bindingSurface');
    const showDisulfides = overlayOn('disulfides');
    const showBindingAromatics = overlayOn('bindingAromatics');
    const colorOnlySticks = !!document.getElementById('colorOnlySticks')?.checked;
    const stickLoops = new Set(overlayCbs.filter(cb => cb.dataset.overlay === 'loop').map(cb => cb.dataset.loop));

    // Catalytic domain / full structure toggle
    const showFullStructure = document.getElementById('showFullStructure')?.checked || false;

    molstarContainer.style.filter = 'blur(3px)';
    // Let the browser paint the blur before the rebuild below blocks the page, so the viewer
    // reacts at once instead of changing only when the new representations are done
    if (!isFirstLoad) await new Promise(resolve => requestAnimationFrame(() => setTimeout(resolve, 0)));
    // Clear focus/selection only on first load — on re-coloring we preserve the user's view
    if (isFirstLoad) {
        instance.visual.interactivityFocus({data:[]});
        instance.visual.clearSelection();
    }

    for (let proteinIndex = 0; proteinIndex < protCount; proteinIndex++){

        let mainData = [];

        // Effective residue limit: catalytic domain end or full length. Both this and
        // the residue numbers below are in the structure's own (original) numbering.
        const domainEnd = (molstarColoring['catalyticDomainEnd'] && molstarColoring['catalyticDomainEnd'][proteinIndex] != null)
            ? molstarColoring['catalyticDomainEnd'][proteinIndex]
            : Infinity;
        const resLimit = showFullStructure ? Infinity : domainEnd;

        // Original residue numbers from the structure file — these do not necessarily
        // start at 1 and may contain gaps, so never derive them from a loop index.
        const protResNums = (molstarColoring['residueNumbers'] && molstarColoring['residueNumbers'][proteinIndex]) || [];

        if (mainRep !== 'none') {
            if (coloringOption === 'group'){
                // each protein in its group's colour (groups of the comparison, see groups.py)
                const color = groupColorOf(proteinIndex);
                for (const resNum of protResNums){
                    if (resNum > resLimit) break;
                    mainData.push({residue_number: resNum, representation: mainRep, representationColor: color});
                }
            } else if (coloringOption === 'none'){
                for (const resNum of protResNums){
                    if (resNum > resLimit) break;
                    mainData.push({residue_number: resNum, representation: mainRep, representationColor: {r:255, g:255, b:255}});
                }
            } else if (coloringOption === 'sequence'){
                let rawdata = [];
                for (let j = 0; j < sequenceColoringOptions.length; j++){
                    rawdata.push(...molstarColoring[coloringOption][sequenceColoringOptions[j]][proteinIndex]);
                }
                rawdata.sort((a,b) => a['residue_number'] - b['residue_number']);

                let rawdataIndex = 0;
                for (let pos = 0; pos < protResNums.length; pos++){
                    const resNum = protResNums[pos];
                    if (resNum > resLimit) break;
                    if (rawdata.length > rawdataIndex && resNum == rawdata[rawdataIndex]['residue_number']){
                        mainData.push({residue_number: resNum, representation: mainRep, representationColor: rawdata[rawdataIndex]['color']});
                        rawdataIndex++;
                    } else {
                        mainData.push({residue_number: resNum, representation: mainRep, representationColor: {r:255, g:255, b:255}});
                    }
                }

            } else if (coloringOption === 'hydrophobicity' || coloringOption === 'flexibility' || coloringOption === 'charge' || coloringOption === 'holeHopping' || coloringOption === 'secondaryStructure' || coloringOption === 'rsa'){
                let chosenOptionArr = molstarColoring[coloringOption][proteinIndex];
                for (let k = 0; k < chosenOptionArr.length; k++){
                    const resNum = chosenOptionArr[k]['residue_number'];
                    if (resNum > resLimit) break;
                    mainData.push({residue_number: resNum, representation: mainRep, representationColor: chosenOptionArr[k]['color']});
                }

            } else if (coloringOption === 'bindingSurface'){
                let bindingSurface = molstarColoring['bindingSurface'][proteinIndex];
                let bindingSurfaceIndex = 0;
                for (let pos = 0; pos < protResNums.length; pos++){
                    const residueId = protResNums[pos];
                    if (residueId > resLimit) break;
                    if (bindingSurface.length === 0){
                        mainData.push({residue_number: residueId, representation: mainRep, representationColor: {r:255,g:255,b:255}});
                        continue;
                    }
                    if (residueId > bindingSurface[bindingSurfaceIndex]['end_residue_number'] && bindingSurfaceIndex < bindingSurface.length - 1){
                        bindingSurfaceIndex++;
                    }
                    if (residueId >= bindingSurface[bindingSurfaceIndex]['start_residue_number'] && residueId <= bindingSurface[bindingSurfaceIndex]['end_residue_number']){
                        mainData.push({residue_number: residueId, representation: mainRep, representationColor: bindingSurface[bindingSurfaceIndex]['color']});
                    } else {
                        mainData.push({residue_number: residueId, representation: mainRep, representationColor: {r:255,g:255,b:255}});
                    }
                }

            } else if (coloringOption === 'loopRegions'){
                let loopRegions = molstarColoring['loopRegions'][proteinIndex];
                let loopRegionIndex = 0;
                for (let pos = 0; pos < protResNums.length; pos++){
                    const residueId = protResNums[pos];
                    if (residueId > resLimit) break;
                    if (loopRegions.length === 0){
                        mainData.push({residue_number: residueId, representation: mainRep, representationColor: {r:255,g:255,b:255}});
                        continue;
                    }
                    if (residueId > loopRegions[loopRegionIndex]['end_residue_number'] && loopRegionIndex < loopRegions.length - 1){
                        loopRegionIndex++;
                    }
                    if (residueId >= loopRegions[loopRegionIndex]['start_residue_number'] && residueId <= loopRegions[loopRegionIndex]['end_residue_number']){
                        mainData.push({residue_number: residueId, representation: mainRep, representationColor: loopRegions[loopRegionIndex]['color']});
                    } else {
                        mainData.push({residue_number: residueId, representation: mainRep, representationColor: {r:255,g:255,b:255}});
                    }
                }
            }
        } // end mainRep !== 'none'
        const mainCount = mainData.length;   // the main representation's entries come first

        // Active site metal ion shown as ball-and-stick — but hidden too when the
        // representation is "None" (so "None" leaves the viewer empty of that protein).
        const metalResNum = molstarColoring['metalResidueNumbers'] && molstarColoring['metalResidueNumbers'][proteinIndex];
        if (metalResNum != null && mainRep !== 'none') {
            mainData.push({residue_number: metalResNum, representation: 'ball-and-stick'});
        }

        // Histidine brace overlay
        if (showBrace) {
            let braceData = molstarColoring['sequence']['brace'][proteinIndex];
            for (let i = 0; i < braceData.length; i++){
                mainData.push({residue_number: braceData[i]['residue_number'], representation: 'ball-and-stick'});
            }
        }

        // Second-sphere overlay
        if (showSecondSphere) {
            const ssCategories = ['Axial position (to Cu)', 'Second-sphere motif', 'Shielding residue',
                                  'Pi-pi stacking', 'Cation-pi', 'H-bond'];
            for (let cat of ssCategories){
                let catData = molstarColoring['sequence'][cat] && molstarColoring['sequence'][cat][proteinIndex]
                    ? molstarColoring['sequence'][cat][proteinIndex] : [];
                for (let i = 0; i < catData.length; i++){
                    mainData.push({residue_number: catData[i]['residue_number'], representation: 'ball-and-stick'});
                }
            }
        }

        // Hole-hopping overlay
        if (showHoleHopping && molstarColoring['holeHopping'] && molstarColoring['holeHopping'][proteinIndex]) {
            // The residues on a path are the ones the hole-hopping colouring does not leave
            // white; as sticks they take the chosen colouring, like the other sets
            let hhData = molstarColoring['holeHopping'][proteinIndex];
            for (let i = 0; i < hhData.length; i++){
                let c = hhData[i]['color'];
                if (!(c.r === 255 && c.g === 255 && c.b === 255)){
                    mainData.push({residue_number: hhData[i]['residue_number'], representation: 'ball-and-stick'});
                }
            }
        }

        // Binding-surface and loop overlays: runs of residue numbers (start..end in the
        // structure's numbering), so take the structure's residues inside each run
        const sticksInRuns = (runs) => {
            for (const run of runs || []){
                for (const resNum of protResNums){
                    if (resNum >= run['start_residue_number'] && resNum <= run['end_residue_number']){
                        mainData.push({residue_number: resNum, representation: 'ball-and-stick'});
                    }
                }
            }
        };
        if (showBindingSurface && molstarColoring['bindingSurface']) {
            sticksInRuns(molstarColoring['bindingSurface'][proteinIndex]);
        }
        if (showBindingAromatics) {
            (bindingAromaticNumbers()[proteinIndex] || []).forEach(resNum => {
                mainData.push({residue_number: resNum, representation: 'ball-and-stick'});
            });
        }
        if (showDisulfides && molstarColoring['disulfides']) {
            sticksInRuns(molstarColoring['disulfides'][proteinIndex]);
        }
        if (stickLoops.size && molstarColoring['loopRegions']) {
            sticksInRuns((molstarColoring['loopRegions'][proteinIndex] || []).filter(run => stickLoops.has(run['group'])));
        }

        // "Colour only the sticks": the main representation white, the sticks in the colour
        // their residue has in the chosen colouring (the metal keeps its own)
        if (colorOnlySticks && mainCount) {
            const colorOf = new Map();
            for (let k = 0; k < mainCount; k++){
                const entry = mainData[k];
                if (entry.representationColor) colorOf.set(entry.residue_number, entry.representationColor);
                mainData[k] = Object.assign({}, entry, {representationColor: {r:255, g:255, b:255}});
            }
            for (let k = mainCount; k < mainData.length; k++){
                const entry = mainData[k];
                if (!entry.representationColor && entry.residue_number !== metalResNum && colorOf.has(entry.residue_number)){
                    mainData[k] = Object.assign({}, entry, {representationColor: colorOf.get(entry.residue_number)});
                }
            }
        }

        allData.push(mainData);
    }

    // Store allData globally so visibility toggling can re-apply correct data
    lastAllData = allData;

    // Visibility keyed by protein index (checkbox value), robust to table sorting.
    let visArr = getProteinVisibility(protCount);

    if (proteinLoadingTime === undefined) {
        // First load
        for (let i = 0; i < allData.length; i++){
            instance.visual.select({data: allData[i], structureNumber: i+1, keepColors: false, keepRepresentations: false});
        }
        const startTime = performance.now();
        await checkInstanceReps(instance, protCount);
        proteinLoadingTime = performance.now() - startTime;
    } else {
        // A recolour/overlay change must rebuild the representations (keepRepresentations:
        // false) — keeping them skips recolouring entirely. The rebuild does not move the
        // camera: Mol*'s automatic re-fit is switched off once the first view is set.
        for (let i = 0; i < allData.length; i++){
            if (visArr[i]) {
                instance.visual.select({data: allData[i], structureNumber: i+1, keepColors: false, keepRepresentations: false});
            } else {
                instance.visual.select({data: [], structureNumber: i+1, keepColors: false, keepRepresentations: false});
            }
        }
        // The camera stays where it is (Mol*'s automatic re-fit is off, see applyViewMode);
        // only its clipping range grows when the scene did (e.g. Full length), so nothing new
        // is cut off by the far plane
        await waitForMolstarFrames(instance.plugin, 3, 800);
        widenCameraClipping(instance.plugin);
    }

    molstarContainer.style.filter = 'blur(0px)';
};

// Structures load with Mol*'s cheap default representation; per-residue colouring is
// applied only on demand through this wrapper. The first time any scheme is applied we
// delete the default polymer component (once) so initializeRepresentations can rebuild it
// via visual.select(). Deferring this is what keeps loading many structures fast — the
// per-residue rebuild ×N structures up front was the source of the many-protein crash.
window._molstarPolymerDeleted = false;
async function applyMolstarColoring(coloringOption, protCount){
    const instance = window.viewerInstance;
    if (!instance || !window.molstarColoring) return;
    if (!window._molstarPolymerDeleted) {
        try {
            const refs = instance.getStructures();
            const reps = window.numberOfDifferentRepresentations || 1;
            for (let i = 0; i < protCount * reps && i < refs.length; i++){
                const structureRef = refs[i]['structureRef']['cell']['sourceRef'];
                instance.deleteStructureComponents(structureRef, ['polymer']);
            }
        } catch (e) { console.warn('applyMolstarColoring: polymer delete failed', e); }
        window._molstarPolymerDeleted = true;
    }
    await initializeRepresentations(window.molstarColoring, instance, coloringOption, protCount);
}
window.applyMolstarColoring = applyMolstarColoring;

// Domain bar traces (catalytic domain / CBMs on the Sequence plot) are tagged with
// meta.domainBar by makePlots.py. Returns the meta when `point` belongs to one, else null.
function domainBarMeta(point){
    const meta = point && point.data && point.data.meta;
    return (meta && meta.domainBar) ? meta : null;
}
window.domainBarMeta = domainBarMeta;

// Highlight a hovered domain / loop band by editing its SVG fill directly — no Plotly
// redraw, so the tooltip that just opened is left untouched. The band paths are the plot's
// only filled traces, rendered as `.js-fill` in the same order as their segIndex. We use
// the same molstar-pink as the residue cursor / 3D highlight so a hovered bar and its 3D
// highlight read as the same colour.
const BAR_HOVER_PINK = 'rgba(255, 20, 147, 0.55)';
let _hoveredDomainBand = null;
// meta.dots (disulfide bonds): the hovered point's trace holds both Cys dots of the bond, so
// both turn pink; its trace is found by point.curveNumber (every trace of these plots is a
// scatter trace).
function domainBandHoverOn(plotDiv, meta, point){
    domainBandHoverOff();
    try {
        let els = [];
        if (meta.dots) {
            const traces = plotDiv.querySelectorAll('.scatterlayer .trace');
            // a whole set (binding surface, disulfides): every trace of it in this plot
            const curves = WHOLE_SET_REGIONS.has(meta.regionKey)
                ? (plotDiv.data || []).map((trace, i) => trace.meta && trace.meta.regionKey === meta.regionKey ? i : -1)
                                      .filter(i => i >= 0)
                : (point ? [point.curveNumber] : []);
            curves.forEach(i => { if (traces[i]) els.push(...traces[i].querySelectorAll('path.point')); });
        } else {
            const el = plotDiv.querySelectorAll('.js-fill')[meta.segIndex];
            if (el) els = [el];
        }
        // Colours are solid; force full opacity so no lingering fill-opacity dilutes it. A tick
        // (line marker) is drawn by its outline, so it turns pink by its stroke; the inline
        // styles are kept to be put back.
        const saved = els.map(el => ({fill: el.style.fill, fillOpacity: el.style.fillOpacity, stroke: el.style.stroke}));
        els.forEach(function (el) {
            el.style.fillOpacity = '1';
            el.style.fill = BAR_HOVER_PINK;
            if (/^M[^Z]*$/i.test(el.getAttribute('d') || '') ) el.style.stroke = BAR_HOVER_PINK;
        });
        if (els.length) _hoveredDomainBand = {els, saved};
    } catch (e) { /* non-fatal */ }
}
function domainBandHoverOff(){
    if (_hoveredDomainBand) {
        try {
            _hoveredDomainBand.els.forEach(function (el, i) {
                const saved = _hoveredDomainBand.saved[i];
                el.style.fill = saved.fill;
                el.style.fillOpacity = saved.fillOpacity;
                el.style.stroke = saved.stroke;
            });
        } catch (e) {}
        _hoveredDomainBand = null;
    }
}

// ── Custom hover tooltip ────────────────────────────────────────────────────
// Plotly draws its tooltip inside the plot's SVG, which is clipped by the short plot
// rows and the horizontally-scrolling column, so tall (multi-line) tooltips get cut off
// or vanish. We hide Plotly's native label (CSS, body.compare-custom-hover) and render
// our own <div> at document-body level — nothing can clip it, and it also lets us style
// it (rounded, shadow, no pointer) to match the site. We reuse the text Plotly already
// computed (reading the hidden native label), so there's no need to replicate each plot's
// hovertemplate.
function _hoverTipEl(){
    let el = document.getElementById('plotHoverTip');
    if (!el){
        el = document.createElement('div');
        el.id = 'plotHoverTip';
        document.body.appendChild(el);
    }
    return el;
}
function _escapeHtml(s){
    return (s || '').replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
}
// Strip any alpha so the tooltip background is fully opaque (readable over busy plots).
function _opaqueColor(c){
    const m = /rgba?\(([^)]+)\)/.exec(c || '');
    if (!m) return c;
    const p = m[1].split(',').map(s => s.trim());
    return `rgb(${p[0]}, ${p[1]}, ${p[2]})`;
}
function _tipTextColor(bg){
    // Choose dark or white text based on the background's luminance.
    const m = (bg || '').match(/rgba?\(([^)]+)\)/);
    if (!m) return '#243b30';
    const [r, g, b] = m[1].split(',').map(Number);
    const lum = (0.299 * r + 0.587 * g + 0.114 * b) / 255;
    return lum > 0.62 ? '#243b30' : '#ffffff';
}
// Reproduce a plot's hover text from the point data by interpolating its hovertemplate.
// This is what the tooltip shows — reading it from the data (not the rendered SVG label)
// avoids the alignment-position line getting dropped, and needs no rAF. Handles the
// specifiers our templates use: %{customdata}, %{x}, %{y} with optional `.Nf` format.
function _interpolateHover(tmpl, pt){
    let s = (tmpl || '').replace(/<extra>[\s\S]*?<\/extra>/g, '');
    s = s.replace(/%\{(customdata|x|y|text)(?::([^}]+))?\}/g, (m, key, fmt) => {
        const v = key === 'customdata' ? pt.customdata
                : key === 'text' ? pt.text
                : pt[key];
        if (v == null) return '';
        if (key !== 'customdata' && key !== 'text' && fmt){
            const dm = /\.(\d+)f/.exec(fmt);
            if (dm) return (+v).toFixed(+dm[1]);
        }
        return String(v);
    });
    return s.trim();
}

function showPlotHoverTip(plotDiv, point, mouseEvent, isDomainBar){
    if (!mouseEvent || !point) return;
    try {
        const tip = _hoverTipEl();
        const tmpl = point.data && point.data.hovertemplate;
        // Content is our own hovertemplate + customdata (only <br>/<b> markup) — safe HTML.
        let html = tmpl ? _interpolateHover(tmpl, point)
                        : String(point.customdata != null ? point.customdata
                                 : (point.text != null ? point.text : '')).trim();
        if (!html){ tip.style.display = 'none'; return; }
        // Background = the trace's per-point hover colour (encodes the property, e.g.
        // charge/hydrophobicity), so the tooltip stays informative. Residue traces carry a
        // bgcolor array; the domain bar carries a single colour. Force it fully opaque so
        // the text is readable over the busy plots underneath.
        const hl = point.data && point.data.hoverlabel;
        const bgc = hl && hl.bgcolor;
        let bg = Array.isArray(bgc) ? (bgc[point.pointNumber] || bgc[0]) : bgc;
        bg = _opaqueColor(bg) || 'rgb(60, 90, 75)';
        tip.innerHTML = html;
        tip.style.background = bg;
        tip.style.color = _tipTextColor(bg);
        tip.style.display = 'block';
        // Vertically centre the tooltip on the hovered plot row, so it sits level with the
        // row (consistent for every residue, not tracking the cursor's y or the curve's
        // height); follow the cursor horizontally, beside it. Clamp to viewport.
        const wrap = plotDiv.closest('.plotWrapper') || plotDiv;
        const wr = wrap.getBoundingClientRect();
        const r = tip.getBoundingClientRect();
        const pad = 10;
        let x = mouseEvent.clientX + 14;
        let y = wr.top + wr.height / 2 - r.height / 2;
        if (x + r.width + pad > window.innerWidth) x = mouseEvent.clientX - r.width - 14;
        if (x < pad) x = pad;
        if (y < pad) y = pad;
        if (y + r.height + pad > window.innerHeight) y = window.innerHeight - r.height - pad;
        tip.style.left = x + 'px';
        tip.style.top = y + 'px';
    } catch (e) { /* non-fatal */ }
}
function hidePlotHoverTip(){
    const el = document.getElementById('plotHoverTip');
    if (el) el.style.display = 'none';
}
window.showPlotHoverTip = showPlotHoverTip;
window.hidePlotHoverTip = hidePlotHoverTip;
// Hide Plotly's native (clip-prone) label; our div replaces it.
document.addEventListener('DOMContentLoaded', () => document.body.classList.add('compare-custom-hover'));
if (document.body) document.body.classList.add('compare-custom-hover');

function hoverSequenceSelection(protCount){
    var plotArr = alignmentPlotDivs();
    plotArr = groupPlotsByType(plotArr, protCount);

    // With 1 structure per protein, structure number for protein i is i+1
    var representationIndexes = [];
    for (let protRep = 0; protRep < protCount; protRep++){
        representationIndexes.push([protRep + 1]);
    }

    for (let typeIndex = 0; typeIndex < plotArr.length; typeIndex++){
        for (plotIndex = 0; plotIndex < plotArr[typeIndex].length; plotIndex++){

            if (typeIndex == 0){
                if (!plotArr[typeIndex][plotIndex].__clickAnnotationBound) {
                    plotArr[typeIndex][plotIndex].__clickAnnotationBound = true;
                    plotArr[typeIndex][plotIndex].on('plotly_clickannotation', async function (data){
                        let alignmentIndex = +data.annotation.x;
                        let residueNum = molstarColoring.residueIndex[structureNum - 1][alignmentIndex - 1];
                        if (!isNaN(residueNum)){
                            let structNum = representationIndexes[structureNum - 1][0];
                            await viewerInstance.visual.interactivityFocus({data:[{residue_number: residueNum}], structureNumber: structNum});
                            await viewerInstance.visual.focus([{residue_number: residueNum}], structNum);
                        }
                    });
                }
            }

            let structureNum = plotIndex + 1;
            // Capture the graph div: inside a Plotly event handler `this` is the event
            // emitter, not the graph div, so passing `this` to Plotly.restyle() throws.
            const plotDiv = plotArr[typeIndex][plotIndex];

            if (!plotDiv.__hoverBound) {
                plotDiv.__hoverBound = true;
                plotDiv.on('plotly_hover', async function (data){
                    // Domain bar: darken the band (direct SVG, no redraw), let Plotly show
                    // its tooltip, and highlight the whole domain in the 3D viewer.
                    const barMeta = domainBarMeta(data.points[0]);
                    if (barMeta) {
                        domainBandHoverOn(plotDiv, barMeta);
                        viewerInstance.visual.clearHighlight();
                        if (barMeta.startRes != null && barMeta.endRes != null){
                            let structNum = representationIndexes[structureNum - 1][0];
                            viewerInstance.visual.highlight({
                                data: [{start_residue_number: barMeta.startRes, end_residue_number: barMeta.endRes}],
                                color:{r:255, g:20, b:147}, structureNumber: structNum});
                        }
                        return;
                    }
                    // Loop / core-strand bar: pink the hovered band (like the domain bar) and
                    // highlight the whole region in 3D (tooltip shown by Plotly).
                    const loopMeta = loopBarMeta(data.points[0]);
                    if (loopMeta) {
                        domainBandHoverOn(plotDiv, loopMeta, data.points[0]);
                        try { viewerInstance.visual.clearHighlight(); } catch (e) {}
                        if (window.highlightLoopGroup) highlightLoopGroup(structureNum - 1, regionGroupOf(loopMeta), loopMeta.regionKey);
                        return;
                    }
                    domainBandHoverOff();
                    viewerInstance.visual.clearHighlight();
                    var point = data.points[0];
                    var alignmentIndex = +point.x;
                    var residueNum = molstarColoring.residueIndex[structureNum - 1][alignmentIndex - 1];
                    if (!isNaN(residueNum)){
                        let structNum = representationIndexes[structureNum - 1][0];
                        viewerInstance.visual.highlight({data: [{residue_number: residueNum}], color:{r:255, g:20, b:147}, structureNumber: structNum});
                    }
                });

                // Leaving the plot: drop the band darkening AND clear the 3D highlight, so
                // it doesn't linger once the cursor is off the residues/domain.
                plotDiv.on('plotly_unhover', function (){
                    domainBandHoverOff();
                    try { viewerInstance.visual.clearHighlight(); } catch(e) {}
                });
            }

            plotArr[typeIndex][plotIndex].on('plotly_click', async function (data){
                // Clicking the domain bar shouldn't zoom the viewer to a residue
                if (domainBarMeta(data.points[0])) return;
                // Clicking a loop / core-strand bar toggles that whole region as sticks (not a
                // residue zoom).
                const loopMeta = loopBarMeta(data.points[0]);
                if (loopMeta) { if (window.toggleLoopGroupSticks) toggleLoopGroupSticks(structureNum - 1, regionGroupOf(loopMeta), loopMeta.regionKey); return; }
                viewerInstance.visual.clearHighlight();
                let alignmentIndex = +data.points[0].x;
                let residueNum = molstarColoring.residueIndex[structureNum - 1][alignmentIndex - 1];
                if (!isNaN(residueNum)){
                    let structNum = representationIndexes[structureNum - 1][0];
                    await viewerInstance.visual.interactivityFocus({data:[{residue_number: residueNum}], structureNumber: structNum});
                    await viewerInstance.visual.focus([{residue_number: residueNum}], structNum);
                };


            });
        };
    };

};

async function changeVisibleRepresentationType(protCount, instance){
    if (!window.molstarColoring || proteinLoadingTime === undefined) return;
    var coloringOption = document.getElementById('coloringSelect')?.value || 'sequence';
    if (coloringOption === 'default') return;   // nothing coloured yet — leave Mol* default
    await applyMolstarColoring(coloringOption, protCount);
};

// --- Per-protein region interaction (Hole Hopping + Secondary Structure and Loops) -------
// Each protein label in these sections is interactive: hover highlights that protein's
// residues (hole-hopping / all loops) in the 3D viewer, click toggles them on as sticks.
// molstarColoring['holeHopping'][pi] is per-residue; ['loopRegions'][pi] is per-loop-segment
// ({group,start_residue_number,end_residue_number,color}); regionResiduesFor handles both.
function regionResiduesFor(pi, key) {
    const mc = window.molstarColoring;
    if (!mc || !mc[key] || !mc[key][pi]) return [];
    const out = [];
    mc[key][pi].forEach(function (h) {
        if (h.start_residue_number != null && h.end_residue_number != null) {
            for (let r = h.start_residue_number; r <= h.end_residue_number; r++) {
                out.push({ residue_number: r, color: h.color });
            }
        } else {
            const c = h.color || {};
            if (!(c.r === 255 && c.g === 255 && c.b === 255)) {   // white = not in the region
                out.push({ residue_number: h.residue_number, color: h.color });
            }
        }
    });
    return out;
}
function highlightRegion(pi, key, color) {
    const inst = window.viewerInstance;
    if (!inst) return;
    const res = regionResiduesFor(pi, key);
    if (!res.length) return;
    try {
        inst.visual.highlight({
            data: res.map(function (h) { return { residue_number: h.residue_number }; }),
            color: color, structureNumber: pi + 1
        });
    } catch (e) {}
}
function clearRegionHighlight() {
    const inst = window.viewerInstance;
    if (!inst) return;
    try { inst.visual.clearHighlight(); } catch (e) {}
}
window._regionSticks = window._regionSticks || {};   // key -> { pi: bool }
function toggleRegionSticks(pi, key, dataPlotType) {
    const inst = window.viewerInstance;
    if (!inst) return;
    const res = regionResiduesFor(pi, key);
    if (!res.length) return;
    const base = (typeof lastAllData !== 'undefined' && lastAllData && lastAllData[pi])
        ? lastAllData[pi].slice() : null;
    if (base === null) { highlightRegion(pi, key, { r: 233, g: 30, b: 140 }); return; }
    const st = window._regionSticks[key] = window._regionSticks[key] || {};
    const on = !st[pi];
    st[pi] = on;
    let data = base;
    if (on) {
        data = base.concat(res.map(function (h) {
            return { residue_number: h.residue_number, representation: 'ball-and-stick',
                     representationColor: h.color };
        }));
    }
    try {
        inst.visual.select({ data: data, structureNumber: pi + 1, keepColors: false, keepRepresentations: false });
    } catch (e) {}
    document.querySelectorAll('.dropdownPlots[data-plot-type="' + dataPlotType + '"] .yaxis-label[data-protein-index="' + pi + '"]')
        .forEach(function (el) { el.classList.toggle('hh-sticks-on', on); });
}
function _ensureRegionHint(section, what) {
    if (section.querySelector('.region-hint')) return;
    const inside = section.querySelector('.insideDropdown') || section;
    const hint = document.createElement('div');
    hint.className = 'region-hint';
    let msg = '💡 Hover a protein name to highlight its ' + what
        + ' in 3D; click the name to show them as sticks.';
    if (what === 'loops') {
        msg += ' Or hover / click a single loop or core strand in the bars above the plot to'
            + ' highlight / show just that one.';
    }
    hint.innerHTML = msg;
    const grid = inside.querySelector('.protein-plots-grid');
    if (grid && grid.parentNode) grid.parentNode.insertBefore(hint, grid);
    else inside.appendChild(hint);
}
// Wire a section's protein labels (idempotent; re-run after each render).
function decorateRegionLabels(dataPlotType, key, highlightColor, whatLabel) {
    const section = document.querySelector('.dropdownPlots[data-plot-type="' + dataPlotType + '"]');
    if (!section) return;
    _ensureRegionHint(section, whatLabel);
    const st = (window._regionSticks && window._regionSticks[key]) || {};
    section.querySelectorAll('.yaxis-label[data-protein-index]').forEach(function (label) {
        const pi = parseInt(label.getAttribute('data-protein-index'), 10);
        if (isNaN(pi)) return;
        const res = regionResiduesFor(pi, key);
        if (!res.length) { label.classList.add('hh-none'); return; }
        label.classList.remove('hh-none');
        label.classList.toggle('hh-sticks-on', !!st[pi]);
        if (label.getAttribute('data-region-wired') === '1') return;
        label.setAttribute('data-region-wired', '1');
        label.classList.add('hh-interactive');
        label.title = 'Hover to highlight in 3D, click to show as sticks';
        label.addEventListener('mouseenter', function () { highlightRegion(pi, key, highlightColor); });
        label.addEventListener('mouseleave', function () { clearRegionHighlight(); });
        label.addEventListener('click', function () { toggleRegionSticks(pi, key, dataPlotType); });
    });
}
window.decorateRegionLabels = decorateRegionLabels;

// The binding-surface ticks and the disulfide dots stand for their whole set: hovering one
// highlights all of the protein's binding-surface residues (or bonded Cys), a click shows them
// all as sticks. Loops and core strands act one at a time.
const WHOLE_SET_REGIONS = new Set(['bindingSurface', 'disulfides']);
function regionGroupOf(meta) {
    return WHOLE_SET_REGIONS.has(meta.regionKey) ? '*' : meta.group;
}
window.regionGroupOf = regionGroupOf;

// --- Per-region interaction (hover/click a single loop or core strand in the bars of the
// Secondary Structure and Loops plot, or a disulfide-bonded Cys dot below it). key: the
// molstarColoring entry, 'loopRegions', 'betaCore' or 'disulfides' (both Cys of a bond).
function loopResiduesForGroup(pi, group, key) {
    const mc = window.molstarColoring;
    key = key || 'loopRegions';
    if (!mc || !mc[key] || !mc[key][pi]) return [];
    const out = [];
    mc[key][pi].forEach(function (h) {
        if ((group === '*' || h.group === group) && h.start_residue_number != null && h.end_residue_number != null) {
            for (let r = h.start_residue_number; r <= h.end_residue_number; r++) {
                out.push({ residue_number: r, color: h.color });
            }
        }
    });
    return out;
}
let _lastLoopHover = null;
function highlightLoopGroup(pi, group, key) {
    const inst = window.viewerInstance;
    if (!inst) return;
    const res = loopResiduesForGroup(pi, group, key);
    if (!res.length) return;
    try {
        inst.visual.highlight({ data: res.map(function (h) { return { residue_number: h.residue_number }; }),
            color: { r: 255, g: 20, b: 147 }, structureNumber: pi + 1 });
    } catch (e) {}
}
window._loopSticks = window._loopSticks || {};   // "pi|key|group" -> bool
function toggleLoopGroupSticks(pi, group, key) {
    const inst = window.viewerInstance;
    if (!inst) return;
    key = key || 'loopRegions';
    if (!loopResiduesForGroup(pi, group, key).length) return;
    const base = (typeof lastAllData !== 'undefined' && lastAllData && lastAllData[pi]) ? lastAllData[pi].slice() : null;
    if (base === null) { highlightLoopGroup(pi, group, key); return; }
    const id = pi + '|' + key + '|' + group;
    window._loopSticks[id] = !window._loopSticks[id];
    // Rebuild this protein with every region currently toggled on (so toggling one keeps the rest).
    const extra = [];
    Object.keys(window._loopSticks).forEach(function (k) {
        if (!window._loopSticks[k]) return;
        const parts = k.split('|');
        if (parseInt(parts[0], 10) !== pi) return;
        loopResiduesForGroup(pi, parts.slice(2).join('|'), parts[1]).forEach(function (h) {
            extra.push({ residue_number: h.residue_number, representation: 'ball-and-stick', representationColor: h.color });
        });
    });
    try {
        inst.visual.select({ data: base.concat(extra), structureNumber: pi + 1, keepColors: false, keepRepresentations: false });
    } catch (e) {}
}
// Loop / core-strand bars (added above the secondary structure by makePlots._addRegionBar,
// tagged meta.loopBar, with meta.regionKey naming the molstarColoring entry).
// Handled in the shared plot hover/click dispatch alongside the domain bar.
function loopBarMeta(point) {
    const meta = point && point.data && point.data.meta;
    return (meta && meta.loopBar) ? meta : null;
}
window.loopBarMeta = loopBarMeta;

// The plots along the alignment (one per protein and plot type): every Plotly plot of the page
// except the selection overview's charts (#catalogueOverview), which have their own hover and
// must not get the residue cursor, the 3D highlight or a place in groupPlotsByType's rows.
function alignmentPlotDivs() {
    return Array.from(document.querySelectorAll('.plotly-graph-div, .js-plotly-plot'))
        .filter(plotDiv => !plotDiv.closest('#catalogueOverview'));
}
window.alignmentPlotDivs = alignmentPlotDivs;

function groupPlotsByType(inputArr, protCount){
    // this is a helper function that takes arr of plots and groups them by type into a matrix
    let typeArr = [];

    for (let i = 0; i < inputArr.length; i += protCount) {
        let chunk = [];

        for (let j = i; j < i + protCount && j < inputArr.length; j++) {
            chunk.push(inputArr[j]);
        };

        typeArr.push(chunk);
    };

    return typeArr
};



function syncScrollbars() {
    const scrollableElements = document.querySelectorAll('.scrollable-plots-column');
    
    let isSyncing = false;

    scrollableElements.forEach(sourceElement => {
        sourceElement.addEventListener('scroll', () => {
            if (isSyncing) return;

            isSyncing = true;
            
            window.requestAnimationFrame(() => {
                const scrollLeft = sourceElement.scrollLeft;

                scrollableElements.forEach(targetElement => {
                    if (targetElement !== sourceElement) {
                        targetElement.scrollLeft = scrollLeft;
                    }
                });
                
                isSyncing = false;
            });
        });
    });
};


// Vertically align each protein label with the data row of its plot. Two things drift:
//   1. Plotly can render a plot row TALLER than the CSS label height, so the label column
//      and plot column stack out of sync and labels creep toward the top of each plot.
//      -> we set each label's height to its plot's actual rendered height so they stack
//         in lockstep.
//   2. Within the row, the data sits wherever Plotly's axis lands.
//      -> we pad the label so its text centres on the data (measured in screen coords, so
//         it's robust to any offset between the two columns).
// For the plots with bars above the residue row (Sequence, Secondary Structure and Loops)
// the target is the residue row (y=0), below the bars: at ymax / (ymax - ymin) down the
// plotting area, from the plot's y-range.
const SS_LOOPS_PLOT_TYPE = 'Secondary Structure and Loops';
const BAR_PLOT_TYPES = ['Sequence', SS_LOOPS_PLOT_TYPE];
window.BAR_PLOT_TYPES = BAR_PLOT_TYPES;
function residueRowFraction(gd){
    const range = gd.layout && gd.layout.yaxis && gd.layout.yaxis.range;
    if (range && range.length === 2 && range[1] > range[0]) return range[1] / (range[1] - range[0]);
    return 1.08 / (1.08 + 0.5);
}
function alignLabelsToData(){
    document.querySelectorAll('.protein-plots-grid').forEach(grid => {
        const labels = grid.querySelectorAll('.yaxis-label');
        const plots = grid.querySelectorAll('.plotWrapper');
        const n = Math.min(labels.length, plots.length);
        for (let i = 0; i < n; i++){
            const gd = plots[i], label = labels[i];
            try {
                const gr = gd.getBoundingClientRect();
                if (!gr.height) continue;
                // Keep the label row height equal to the plot row height (lockstep stacking).
                const hpx = gr.height + 'px';
                if (label.style.height !== hpx) label.style.height = hpx;

                const area = gd.querySelector('.nsewdrag') || gd.querySelector('.gridlayer')
                          || gd.querySelector('.cartesianlayer');
                if (!area) continue;
                const ar = area.getBoundingClientRect();
                if (!ar.height) continue;
                // Plots with bars above the residue row align their label to that row, not
                // the middle.
                const pt = gd.getAttribute('data-plot-type');
                const target = BAR_PLOT_TYPES.includes(pt) ? (ar.top + residueRowFraction(gd) * ar.height)
                                                           : (ar.top + ar.height / 2);
                // Measure after setting height so lr reflects the corrected stacking.
                const lr = label.getBoundingClientRect();
                if (!lr.height) continue;
                // With align-items:center the text centre is at lr.top + H/2 at pad=0; each
                // px of padding-top shifts it down ½px, so pad = 2·(target − current centre).
                const pad = Math.max(0, Math.round(2 * (target - lr.top - lr.height / 2)));
                const px = pad + 'px';
                if (label.style.paddingTop !== px) label.style.paddingTop = px;
            } catch (e) { /* non-fatal */ }
        }
    });
}
window.alignLabelsToData = alignLabelsToData;

// Coalesce re-alignment requests to one per frame (idempotent, so calling it from
// plotly_afterplot / resize / observers can't cause flicker).
let _alignRaf = null;
function scheduleLabelAlign(){
    if (_alignRaf) return;
    _alignRaf = requestAnimationFrame(() => { _alignRaf = null; alignLabelsToData(); });
}
window.scheduleLabelAlign = scheduleLabelAlign;
window.addEventListener('resize', scheduleLabelAlign);
window.addEventListener('load', () => { scheduleLabelAlign(); setTimeout(scheduleLabelAlign, 800); });

// Keep labels aligned no matter when the plots reach their final size (fonts, scrollbar,
// the Mol* viewer settling can all reflow late). A ResizeObserver on each plot catches
// size changes; the timed passes below cover reflows that don't change the wrapper size.
let _labelResizeObs = null;
function watchLabelAlignment(){
    try {
        if (typeof ResizeObserver === 'undefined') return;
        if (_labelResizeObs) _labelResizeObs.disconnect();
        _labelResizeObs = new ResizeObserver(scheduleLabelAlign);
        document.querySelectorAll('.protein-plots-grid .plotWrapper').forEach(p => _labelResizeObs.observe(p));
    } catch (e) { /* non-fatal */ }
}
window.watchLabelAlignment = watchLabelAlignment;

async function makePlotsFromJSON(jsonArr){
    let plotJSONs = [];
    for (let j = 0; j<jsonArr.length; j++){
        plotJSONs.push(...jsonArr[j].plots_html);
    };
    // A large comparison opens with every plot section closed: nothing to draw yet, and any
    // wrapper present belongs to a section the visitor has opened in the meantime.
    if (!plotJSONs.length) return;
    let plotWrappers = document.querySelectorAll("div.plotWrapper");

    // Create all plots and wait for them to be ready
    const plotPromises = [];
    for (let i = 0; i < plotWrappers.length; i++){
        let parsedJSON = JSON.parse(plotJSONs[i]);
        // displayModeBar: false — hide Plotly's hover toolbar; it isn't useful here and
        // covers part of the (short) plots when it appears.
        const plotPromise = Plotly.newPlot(plotWrappers[i], parsedJSON.data, parsedJSON.layout, {displayModeBar: false, responsive: false, scrollZoom: false, doubleClick: false});
        plotPromises.push(plotPromise);
    }
    
    // Wait for all plots to be created
    await Promise.all(plotPromises);

    // Align the protein labels to the plots' data rows. Run now, again after the layout
    // settles, watch for later reflows (ResizeObserver), and re-align on every redraw
    // (plotly_afterplot) so it can't drift back after the page finishes loading.
    scheduleLabelAlign();
    [200, 600, 1500, 3000].forEach(t => setTimeout(scheduleLabelAlign, t));
    watchLabelAlignment();

    // Set up afterplot listeners for each plot to initialize cursor when ready
    for (let i = 0; i < plotWrappers.length; i++){
        const plotDiv = plotWrappers[i].querySelector('.js-plotly-plot') || plotWrappers[i];
        if (plotDiv && typeof plotDiv.on === 'function') {
            plotDiv.on('plotly_afterplot', function() {
                scheduleLabelAlign();
                // When plot is ready, initialize cursor
                setTimeout(function() {
                    if (window.initializeVerticalCursor) {
                        window.initializeVerticalCursor();
                    }
                }, 50);
            });
        }
    }
    
    // Initialize cursor after all plots are ready (with multiple attempts)
    setTimeout(function() {
        if (window.initializeVerticalCursor) {
            console.log('Calling initializeVerticalCursor after plots created');
            window.initializeVerticalCursor();
        }
    }, 500);
    
    // Also try again after a longer delay to catch any late-loading plots
    setTimeout(function() {
        if (window.initializeVerticalCursor) {
            window.initializeVerticalCursor();
        }
    }, 1500);
};



// Hide (but don't delete) non-protein Molstar components: water, ligands, ions, etc.
// Returns a Promise that resolves once the Mol* state update is committed.
async function hideNonProteinComponents(viewerInstance) {
    try {
        const plugin = viewerInstance.plugin;
        if (!plugin?.managers?.structure) return;

        const structures = plugin.managers.structure.hierarchy.current.structures;
        const componentsToHide = [];

        for (const structure of structures) {
            for (const comp of structure.components) {
                const typeName = comp.cell?.params?.values?.type?.name;
                // polymer was already deleted; skip it and any unrecognised types
                if (!typeName || typeName === 'polymer') continue;
                // Only queue currently-visible components
                if (!comp.cell?.state?.isHidden) {
                    componentsToHide.push(comp);
                }
            }
        }

        if (componentsToHide.length > 0) {
            await plugin.managers.structure.component.toggleVisibility(componentsToHide);
        }
    } catch (e) {
        console.warn('hideNonProteinComponents: could not hide via plugin API', e);
    }
}

// Grow the camera's clipping radius (around its target) to hold the whole visible scene,
// without moving the camera
function widenCameraClipping(plugin) {
    try {
        const canvas3d = plugin && plugin.canvas3d;
        const sphere = canvas3d && canvas3d.boundingSphereVisible;
        if (!sphere || !sphere.radius) return;
        const state = canvas3d.camera.state;
        const need = Math.hypot(state.target[0] - sphere.center[0], state.target[1] - sphere.center[1],
                                state.target[2] - sphere.center[2]) + sphere.radius;
        if (need > state.radius) {
            canvas3d.camera.setState({radius: need, radiusMax: Math.max(need, state.radiusMax || 0)}, 0);
        }
    } catch (e) { /* non-fatal */ }
}

// Wait for n rendered frames from Mol*'s canvas3d, with a ms fallback timeout.
// This is needed before camera reset so the scene bounding box reflects the latest visibility state.
function waitForMolstarFrames(plugin, n, fallbackMs) {
    return new Promise(resolve => {
        let done = false;
        const finish = () => { if (!done) { done = true; resolve(); } };
        try {
            let frames = 0;
            const sub = plugin.canvas3d?.didDraw?.subscribe?.(() => {
                if (++frames >= n) { sub?.unsubscribe?.(); finish(); }
            });
            if (!sub) finish(); // didDraw not available — fall through to timeout
        } catch(e) { finish(); }
        setTimeout(finish, fallbackMs);
    });
}

// --- Views from the binding-surface plane -------------------------------------------------
// _bindingPlanes()[i] = {normal (unit, towards the solvent), copper, depth} of protein i
// (pipeline step 7c), in the aligned coordinates the viewer shows, so one camera fits all
// structures: the copper and the normal are averaged over the visible proteins.
//   top  - looking straight down onto the copper face, the copper in the centre, the brace
//          horizontal with the N-His on the left (default)
//   side - from the side, the copper face up, the brace horizontal
//   site - close up on the copper, from above
const _v3 = {
    add: (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]],
    sub: (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]],
    scale: (a, k) => [a[0] * k, a[1] * k, a[2] * k],
    dot: (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2],
    cross: (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]],
    len: a => Math.hypot(a[0], a[1], a[2]),
    unit: a => { const l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / l, a[1] / l, a[2] / l]; },
};
function _bindingPlanes(){
    if (!Array.isArray(window._bindingPlaneData)) {
        try { window._bindingPlaneData = JSON.parse(document.getElementById('bindingPlanes').textContent) || []; }
        catch (e) { window._bindingPlaneData = []; }
    }
    return window._bindingPlaneData;
}
function viewModeCamera(mode, plugin){
    const planes = _bindingPlanes();
    const visible = getProteinVisibility(planes.length);
    let shown = planes.filter((plane, i) => plane && visible[i]);
    if (!shown.length) shown = planes.filter(Boolean);
    if (!shown.length) return null;
    let copper = [0, 0, 0], normal = [0, 0, 0], brace = [0, 0, 0];
    shown.forEach(plane => {
        copper = _v3.add(copper, plane.copper);
        normal = _v3.add(normal, plane.normal);
        if (plane.brace) brace = _v3.add(brace, _v3.unit(plane.brace));
    });
    copper = _v3.scale(copper, 1 / shown.length);
    if (_v3.len(normal) < 1e-3) return null;
    normal = _v3.unit(normal);
    // The screen's left -> right: the brace (N-His -> Internal-His ring centres) laid onto the
    // face, so the brace lies horizontal with the N-His on the left; without a brace, the x axis
    // of the aligned frame
    const inPlane = v => _v3.sub(v, _v3.scale(normal, _v3.dot(v, normal)));
    let right = inPlane(brace);
    if (_v3.len(right) < 1e-3) right = inPlane(Math.abs(normal[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]);
    right = _v3.unit(right);
    // Mol*'s screen right = dir x up, so looking down the normal (dir = -normal) the up that
    // puts `right` to the right is normal x right
    const across = _v3.unit(_v3.cross(normal, right));

    const canvas3d = plugin.canvas3d;
    const sphere = (canvas3d.boundingSphereVisible || canvas3d.boundingSphere || {});
    const center = sphere.center ? [sphere.center[0], sphere.center[1], sphere.center[2]] : copper;
    const sceneRadius = sphere.radius || 30;
    // fit = the radius to frame; clip = the radius Mol* keeps unclipped around the target;
    // dir = from the camera to the target. Top and Side frame the whole visible scene: the
    // scene's bounding sphere plus how far its centre lies off the line of sight, so the whole
    // protein is in view wherever the target sits.
    const offCentre = _v3.len(_v3.sub(copper, center));
    const clip = sceneRadius + offCentre;
    const wholeScene = (target, dir) => {
        const offset = _v3.sub(center, target);
        const sideways = _v3.sub(offset, _v3.scale(dir, _v3.dot(offset, dir)));
        return sceneRadius + _v3.len(sideways);
    };
    if (mode === 'side') {
        // looking along `across` with the face up: dir x normal = right again, the brace
        // horizontal as in the top view
        const dir = across;
        return {target: center, fit: wholeScene(center, dir), clip: clip, up: normal, dir: dir};
    }
    const down = _v3.scale(normal, -1);
    if (mode === 'site') {
        return {target: copper, fit: 14, clip: clip, up: across, dir: down};
    }
    // top: centred on the copper
    return {target: copper, fit: wholeScene(copper, down), clip: clip, up: across, dir: down};
}
// Apply a view (mode null = the one ticked); without a plane, Mol*'s own reset
// at once on a click; on page load (instance given) once the scene is complete, so its
// bounding sphere holds every structure
async function applyViewMode(mode, instance){
    const onLoad = !!instance;
    instance = instance || window.viewerInstance;
    const plugin = instance && instance.plugin;
    if (!plugin || !plugin.canvas3d) return;
    if (!mode) mode = (document.querySelector('input[name="viewMode"]:checked') || {}).value || 'top';
    if (onLoad) await waitForSceneSettled(plugin);
    // From now on only explicit actions move the camera (these views, Mol*'s reset button):
    // by default Mol* re-fits the camera whenever newly drawn objects fall outside the view,
    // so showing sticks or recolouring zoomed the view in or out, depending on where it was
    try { plugin.canvas3d.setProps({camera: {manualReset: true}}); } catch (e) {}
    const view = viewModeCamera(mode, plugin);
    if (!view) { plugin.canvas3d.requestCameraReset({ durationMs: 400 }); return; }
    try {
        // Camera.getFocus keeps the current viewing direction (it only flips it), so the
        // position is set here: back from the target along dir, at Mol*'s distance for `fit`
        const camera = plugin.canvas3d.camera;
        const distance = camera.getTargetDistance(view.fit);
        const position = _v3.sub(view.target, _v3.scale(view.dir, distance));
        camera.setState(Object.assign({}, camera.getSnapshot(), {
            target: view.target, position: position, up: view.up,
            radius: view.clip, radiusMax: view.clip,
        }), 0);   // instant: an animated move is lost if a recolour saves the camera half-way
    } catch (e) {
        console.warn('applyViewMode failed', e);
        plugin.canvas3d.requestCameraReset({ durationMs: 400 });
    }
    if (onLoad) refitWhileSceneGrows(plugin);
}
window.applyViewMode = applyViewMode;

// After the load view: for up to a minute, if the scene still grows (a very large set can take
// longer than waitForSceneSettled's limit), fit the ticked view again — unless the viewer was
// touched or a view button clicked meanwhile.
function refitWhileSceneGrows(plugin) {
    const viewer = document.getElementById('myViewer');
    let active = true;
    const stop = () => { active = false; };
    if (viewer) {
        viewer.addEventListener('pointerdown', stop, {capture: true, once: true});
        viewer.addEventListener('wheel', stop, {capture: true, once: true});
    }
    document.querySelectorAll('input[name="viewMode"]').forEach(el => el.addEventListener('click', stop, {once: true}));
    let fitted = (plugin.canvas3d.boundingSphereVisible || {}).radius || 0;
    const t0 = performance.now();
    const check = () => {
        if (!active || performance.now() - t0 > 60000) return;
        const radius = (plugin.canvas3d.boundingSphereVisible || {}).radius || 0;
        if (radius > fitted * 1.05) {
            fitted = radius;
            waitForSceneSettled(plugin).then(() => { if (active) applyViewMode(null); });
        }
        setTimeout(check, 500);
    };
    setTimeout(check, 500);
}

// Resolve once the Mol* scene is complete: nothing left in the canvas's commit queue (Mol*
// adds a large scene in batches over many frames), the plugin's state not updating, and the
// visible bounding sphere unchanged over a few checks. With many structures a fixed number of
// frames was too few: the view was fitted to a partial scene and came out zoomed in.
function waitForSceneSettled(plugin, timeoutMs = 20000) {
    return new Promise(resolve => {
        const t0 = performance.now();
        let lastRadius = -1, stable = 0;
        const check = () => {
            try {
                const canvas3d = plugin.canvas3d;
                const queued = canvas3d.commitQueueSize && canvas3d.commitQueueSize.value;
                const updating = plugin.behaviors && plugin.behaviors.state && plugin.behaviors.state.isUpdating
                    && plugin.behaviors.state.isUpdating.value;
                const radius = (canvas3d.boundingSphereVisible || {}).radius || 0;
                if (!queued && !updating && radius > 0 && Math.abs(radius - lastRadius) < 0.01) stable++;
                else stable = 0;
                lastRadius = radius;
            } catch (e) { stable++; }
            if (stable >= 3 || performance.now() - t0 > timeoutMs) resolve();
            else setTimeout(check, 100);
        };
        check();
    });
}

// Reset camera to show all visible content — mirrors what the built-in reset button does.
async function resetMolstarCamera(viewerInstance) {
    const plugin = viewerInstance.plugin;
    if (!plugin?.canvas3d) return;
    // Wait for the scene bounding box to reflect the current visibility state
    await waitForMolstarFrames(plugin, 3, 800);
    plugin.canvas3d.requestCameraReset({ durationMs: 400 });
}

// Mol* names a structure after its mmCIF `_entry.id`, falling back to the label of the
// data cell it was downloaded from — for us a blob: URL, which is what leaks into the
// sequence panel header ("Sequence of blob:https…"). The PDB HEADER id code that would
// supply `_entry.id` is capped at 4 characters, so an AlphaFold code like AF-B6EQB6 can
// never fit there. Instead, give the download cell an explicit label: Mol*'s Download
// transform accepts one, and it flows through to the structure name everywhere.
//
// Must be called before viewerInstance.render(), which kicks off the first load.
function labelStructureDownloads(viewerInstance, urlToLabel) {
    const originalLoad = viewerInstance.load.bind(viewerInstance);

    viewerInstance.load = function (params, fullLoad) {
        try {
            const plugin = viewerInstance.plugin;   // exists by the time load() runs
            const builder = plugin && plugin.builders && plugin.builders.data;
            if (builder && !builder.__labelPatched) {
                builder.__labelPatched = true;
                const originalDownload = builder.download.bind(builder);
                builder.download = function (downloadParams, options) {
                    const url = typeof downloadParams.url === 'string'
                        ? downloadParams.url
                        : downloadParams.url && downloadParams.url.url;
                    const label = urlToLabel[url];
                    return originalDownload(
                        label ? Object.assign({}, downloadParams, {label: label}) : downloadParams,
                        options
                    );
                };
            }
        } catch (e) {
            console.warn('Could not label structure downloads; Mol* will fall back to the blob URL', e);
        }
        return originalLoad(params, fullLoad);
    };
}

// ── Mol* hover tooltip ──────────────────────────────────────────────────────
// Mol*'s stock label prints the entry followed by "Model 1" and "Instance 1",
// which are always 1 for these single-copy structures. Replace the label providers
// with one that shows only what's meaningful here: the protein, the structure it
// came from, and the residue under the pointer.

// A Mol* OrderedSet is either a SortedArray (array-like of positions) or an Interval,
// which packs its [start, end) bounds as two int32s inside a float64.
const _tupleBuf = new ArrayBuffer(8);
const _tupleI32 = new Int32Array(_tupleBuf);
const _tupleF64 = new Float64Array(_tupleBuf);

function forEachOrderedSetIndex(indices, limit, fn) {
    if (indices == null) return;
    if (typeof indices === 'number') {
        _tupleF64[0] = indices;
        const start = _tupleI32[0];
        const end = Math.min(_tupleI32[1], start + limit);
        for (let i = start; i < end; i++) fn(i);
    } else {
        const end = Math.min(indices.length, limit);
        for (let i = 0; i < end; i++) fn(indices[i]);
    }
}

// Index of the loaded structure a loci belongs to, or -1.
function structureIndexOfLoci(plugin, loci) {
    const structures = plugin.managers.structure.hierarchy.current.structures || [];
    for (let i = 0; i < structures.length; i++) {
        const data = structures[i] && structures[i].cell && structures[i].cell.obj
            ? structures[i].cell.obj.data : null;
        if (data && loci.structure && data.root === loci.structure.root) return i;
    }
    return -1;
}

function buildLociLabel(plugin, loci, proteinCodes, structureCodes) {
    if (!loci || loci.kind !== 'element-loci' || !loci.elements || loci.elements.length === 0) {
        return undefined;
    }

    const first = loci.elements[0];
    const hierarchy = first.unit.model.atomicHierarchy;

    // Hover granularity is per-residue, so the loci holds every atom of the residue.
    // Count distinct residues (capped — a select-all can cover the whole structure).
    const residues = new Set();
    for (const element of loci.elements) {
        const h = element.unit.model.atomicHierarchy;
        forEachOrderedSetIndex(element.indices, 200, pos => {
            residues.add(h.residueAtomSegments.index[element.unit.elements[pos]]);
        });
    }

    let firstPos = -1;
    forEachOrderedSetIndex(first.indices, 1, pos => { firstPos = pos; });
    if (firstPos < 0) return undefined;

    const atom = first.unit.elements[firstPos];
    const residueRow = hierarchy.residueAtomSegments.index[atom];
    const chainRow = hierarchy.chainAtomSegments.index[atom];

    const compId = hierarchy.atoms.label_comp_id.value(atom);
    const seqId = hierarchy.residues.auth_seq_id.value(residueRow);
    const chain = hierarchy.chains.auth_asym_id.value(chainRow);

    const structIdx = structureIndexOfLoci(plugin, loci);
    const protein = structIdx >= 0 ? (proteinCodes[structIdx] || '') : '';
    const code = structIdx >= 0 ? (structureCodes[structIdx] || '') : '';

    const heading = [
        protein ? `<b>${protein}</b>` : '',
        code ? `<small>${code}</small>` : ''
    ].filter(Boolean).join(' ');

    const extra = residues.size > 1 ? ` <small>[+ ${residues.size - 1} more]</small>` : '';
    const residue = `${compId} ${seqId}${chain ? ` | Chain ${chain}` : ''}${extra}`;

    return heading ? `${heading}<br/>${residue}` : residue;
}

function setupLociLabels(viewerInstance, proteinCodes, structureCodes) {
    try {
        const plugin = viewerInstance.plugin;
        const lociLabels = plugin && plugin.managers && plugin.managers.lociLabels;
        if (!lociLabels) return;

        lociLabels.clearProviders();
        lociLabels.addProvider({
            priority: 100,
            label: (loci) => {
                try {
                    return buildLociLabel(plugin, loci, proteinCodes, structureCodes);
                } catch (e) {
                    return undefined;   // never break the tooltip
                }
            }
        });
    } catch (e) {
        console.warn('setupLociLabels: keeping default Mol* tooltip', e);
    }
}

// ── Mol* hover -> plot cursor ───────────────────────────────────────────────
// The forward direction (plot hover -> Mol* highlight) lives in hoverSequenceSelection.
// Here we do the reverse: hovering a residue in the 3D viewer moves the vertical
// cursor to that residue's alignment column in every plot.

// residue number -> alignment column (1-based, as used for the plot x axis), per protein.
// molstarColoring.residueIndex[p] lists the residue number of each alignment column,
// with 'NaN' at gaps, so this is just its inverse.
let residueToAlignmentIndex = null;

function buildResidueToAlignmentIndex(molstarColoring) {
    const byProtein = [];
    const residueIndex = molstarColoring['residueIndex'] || [];
    for (let p = 0; p < residueIndex.length; p++) {
        const map = new Map();
        const cols = residueIndex[p] || [];
        for (let col = 0; col < cols.length; col++) {
            const resNum = cols[col];
            if (resNum === 'NaN' || resNum == null) continue;
            map.set(Number(resNum), col + 1);  // plot x is 1-based
        }
        byProtein.push(map);
    }
    return byProtein;
}

// Which loaded structure the pointer is currently over. Mol*'s hover stream carries the
// structure, while PDBe-Molstar's mouseover event carries the residue number — we need both.
let hoveredProteinIndex = null;

function setupMolstarToPlotHover(viewerInstance, molstarColoring) {
    const plugin = viewerInstance.plugin;
    if (!plugin) return;

    residueToAlignmentIndex = buildResidueToAlignmentIndex(molstarColoring);

    // Track which structure is under the pointer. PDBe-Molstar's mouseover event
    // identifies the residue but not which of the loaded structures it belongs to.
    try {
        plugin.behaviors.interaction.hover.subscribe(e => {
            const loci = e?.current?.loci;
            const structure = loci?.structure;
            if (!structure) { hoveredProteinIndex = null; return; }

            const structures = plugin.managers.structure.hierarchy.current.structures || [];
            hoveredProteinIndex = null;
            for (let i = 0; i < structures.length; i++) {
                const data = structures[i]?.cell?.obj?.data;
                if (data && data.root === structure.root) {
                    hoveredProteinIndex = i;
                    break;
                }
            }
        });
    } catch (err) {
        console.warn('Mol* hover subscription failed; plot cursor sync disabled', err);
        return;
    }

    const container = document.getElementById('myViewer');
    if (!container) return;

    container.addEventListener('PDB.molstar.mouseover', (e) => {
        const resNum = e.eventData && e.eventData.auth_seq_id;
        if (resNum == null || hoveredProteinIndex == null) return;
        const map = residueToAlignmentIndex[hoveredProteinIndex];
        if (!map) return;
        const alignmentIndex = map.get(Number(resNum));
        // Ligands, ions and waters have no alignment column — leave the cursor alone
        if (alignmentIndex == null) return;
        if (window.setPlotCursor) window.setPlotCursor(alignmentIndex);
    });

    container.addEventListener('PDB.molstar.mouseout', () => {
        if (window.clearPlotCursor) window.clearPlotCursor();
    });
}

// Map each plot section's type to the matching Mol* coloring option, so a button in the
// plot section can colour the 3D structures by the property being viewed.
const PLOT_TYPE_TO_COLORING = {
    'Sequence': 'sequence',
    'Flexibility': 'flexibility',
    'Hydrophobicity': 'hydrophobicity',
    'Charge': 'charge',
    'Secondary Structure and Loops': 'secondaryStructure',
    'Relative Surface Accessibility': 'rsa',
    'Binding Surface': 'bindingSurface',
    'Hole Hopping': 'holeHopping',
};
function colorStructureByPlot(plotType){
    colorStructureBy(PLOT_TYPE_TO_COLORING[plotType]);
}
window.colorStructureByPlot = colorStructureByPlot;
// Colour the 3D structures by one of the Mol* colouring options (coloringSelect values).
function colorStructureBy(opt){
    if (!opt || !window.viewerInstance || !window.molstarColoring) return;
    const sel = document.getElementById('coloringSelect');
    if (sel) sel.value = opt;   // keep the dropdown in sync
    const protCount = (window.protIdList || []).length;
    applyMolstarColoring(opt, protCount);
}
window.colorStructureBy = colorStructureBy;

console.log('compare_protein.js file loaded');


