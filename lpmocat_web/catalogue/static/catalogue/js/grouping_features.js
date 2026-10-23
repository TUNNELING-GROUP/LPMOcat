/* Compare page, selection overview: "Which features separate the groups?" — for the groups the
   selection tree's distance cut forms (selection_tree.js), ranks the proteins' features by how
   well they separate the groups: eta squared for numbers (the share of the variance that lies
   between the groups), Cramér's V for categories. It says what the groups have in common beyond
   their sequence; the tree is built from sequence distance, so no feature *causes* a grouping.

   Guards against reading too much into few proteins:
   - a group of one protein counts as a group (a feature it stands out in says what makes that
     protein different), but at least one group must hold two or more proteins: with only
     singletons every feature would separate them perfectly
   - chance level: the group labels are shuffled 1000 times (group sizes kept, singletons
     included, so a lone protein's pull is part of the chance level too) and the 95th
     percentile of the best score over ALL features is drawn as a line; a feature below it
     does no better than the best of many features on random groups
   - features that move together across the compared proteins (|r|, Cramér's V or the
     correlation ratio >= REDUNDANT) are listed under the highest-ranked one, since they are
     one observation, not several
   - each feature is tagged: needs the 3D structure / from the sequence alone (composition
     that the sequence distance itself partly measures) / annotation (papers, taxonomy, data
     source)

   - set features (the hole-hopping positions: the alignment columns of the path residues, with
     or without the residue) compare two proteins by the Jaccard distance and score PERMANOVA R²
     for groups, the distance in the Mantel test, and |2 AUC - 1| of the within-pick against the
     pick-to-reference distances for picked proteins

   "Which features follow the tree?" needs no cut: over every pair of proteins it correlates (Spearman)
   their distance in the tree with how different they are in each feature, a Mantel test, with
   the chance level from shuffling the tree's proteins (analyseTree).

   "Which features discriminate the picked proteins from …" compares the proteins picked in the tree (rows or
   branches; none picked = all the compared ones) with a reference set: the rest of the
   comparison, the rest of the catalogue in the same families (the default against the
   catalogue: against everything, the answer is mostly the family) or the whole rest of the
   catalogue (analysePicked). Numbers score |2 AUC - 1|, how cleanly the values of the picked
   proteins sit above or below the reference (1 = every picked value beyond every reference
   value); categories the total variation distance between the two sets' category shares (for a
   single protein: the share of the reference that does not share its value). Chance level as
   for the groups, from shuffling which proteins are picked (sizes kept). A single picked
   protein is a description, not a finding: with one value many features look extreme.

   All three run only when their button is pressed; moving the slider or changing the pick
   afterwards marks a result as stale. catalogue_overview.js attaches it (GroupingFeatures.attach). */
(function () {
    'use strict';

    const PERMUTATIONS = 1000;
    const REDUNDANT = 0.8;
    const MIN_SUBSTRATE = 8;    // proteins reported on a substrate before it becomes a feature
    const SHOWN = 12;           // leading features listed before "Show all"
    const TAGS = {
        structure: { label: 'needs the 3D structure', color: '#2f6f5e' },
        sequence: { label: 'from the sequence alone', color: '#9aa9a1' },
        annotation: { label: 'annotation', color: '#c9a227' },
    };

    function el(tag, className, text) {
        const node = document.createElement(tag);
        if (className) node.className = className;
        if (text != null) node.textContent = text;
        return node;
    }

    // ---- the features (all from the overview payload the page already has) ----

    function featureList(data) {
        const scale = data.hydrophobicity_scales[0].name;       // GRAVY (Kyte–Doolittle)
        const m = p => p.measures || {};
        const set = (p, name) => ((m(p).sets || {})[name] || {});
        const nd = new Set([data.not_determined, '', null, undefined]);
        const cat = v => (nd.has(v) ? undefined : v);
        const list = [
            // catalytic domain, as composition: computable from the sequence alone
            { label: 'Length', part: 'catalytic domain', tag: 'sequence', unit: 'residues', get: p => (m(p).length || {}).cd },
            { label: 'Mass', part: 'catalytic domain', tag: 'sequence', unit: 'kDa', get: p => (m(p).mass || {}).cd },
            { label: 'Net charge', part: 'catalytic domain', tag: 'sequence', unit: 'e', get: p => (m(p).charge || {}).cd },
            { label: 'GRAVY', part: 'catalytic domain', tag: 'sequence', get: p => ((m(p).hydrophobicity || {})[scale] || {}).cd },
            // residue sets the structure defines
            { label: 'Residues', part: 'surface', tag: 'structure', unit: 'residues', get: p => (m(p).surface_length || {}).cd },
            { label: 'Net charge', part: 'surface', tag: 'structure', unit: 'e', get: p => (m(p).surface_charge || {}).cd },
            { label: 'GRAVY', part: 'surface', tag: 'structure', get: p => ((m(p).surface_hydrophobicity || {})[scale] || {}).cd },
            { label: 'Residues', part: 'all loops', tag: 'structure', unit: 'residues', get: p => set(p, 'loops').residues },
            { label: 'Net charge', part: 'all loops', tag: 'structure', unit: 'e', get: p => set(p, 'loops').charge },
            { label: 'GRAVY', part: 'all loops', tag: 'structure', get: p => (set(p, 'loops').hydrophobicity || {})[scale] },
            { label: 'Residues', part: 'substrate-binding', tag: 'structure', unit: 'residues', get: p => set(p, 'binding').residues },
            { label: 'Net charge', part: 'substrate-binding', tag: 'structure', unit: 'e', get: p => set(p, 'binding').charge },
            { label: 'GRAVY', part: 'substrate-binding', tag: 'structure', get: p => (set(p, 'binding').hydrophobicity || {})[scale] },
        ];
        // flexibility of all loops together and of their substrate-binding residues: each loop's
        // mean RMSF weighted by its residue count ([all, surface, binding] triples in loop_values)
        const loopRmsf = (p, part) => {
            let sum = 0, count = 0;
            Object.values(p.loop_values || {}).forEach(lv => {
                const r = (lv.rmsf || [])[part], k = (lv.residues || [])[part];
                if (typeof r === 'number' && k) { sum += r * k; count += k; }
            });
            return count ? sum / count : undefined;
        };
        list.push(
            { label: 'Flexibility (RMSF)', part: 'all loops', tag: 'structure', unit: 'Å', get: p => loopRmsf(p, 0) },
            { label: 'Flexibility (RMSF)', part: 'substrate-binding residues of the loops', tag: 'structure', unit: 'Å', get: p => loopRmsf(p, 2) },
        );
        // per loop: [all, surface, binding] triples in loop_values; hydrophobicity per scale
        (data.loop_names || []).forEach(loop => {
            const lv = p => (p.loop_values || {})[loop] || {};
            const first = arr => (Array.isArray(arr) ? arr[0] : undefined);
            list.push(
                { label: 'Length', part: loop, tag: 'structure', unit: 'residues', get: p => ((p.loops || {})[loop] || [])[0] },
                { label: 'Flexibility (RMSF)', part: loop, tag: 'structure', unit: 'Å', get: p => first(lv(p).rmsf) },
                { label: 'Net charge', part: loop, tag: 'structure', unit: 'e', get: p => first(lv(p).charge) },
                { label: 'GRAVY', part: loop, tag: 'structure', get: p => first(first(lv(p).hydrophobicity)) },
                { label: 'Binding residues', part: loop, tag: 'structure', unit: 'residues', get: p => (lv(p).residues || [])[2] },
                { label: 'Net charge of the binding residues', part: loop, tag: 'structure', unit: 'e', get: p => (lv(p).charge || [])[2] },
            );
        });
        // overview.structural_features (p.features); a missing value is null there
        const fx = (p, key) => { const v = (p.features || {})[key]; return v == null ? undefined : v; };
        list.push(
            { label: 'Disulfide bonds', part: 'catalytic domain', tag: 'structure', unit: 'bonds', get: p => fx(p, 'disulfides') },
            { label: 'Disulfide pattern', part: 'regions joined', tag: 'structure', kind: 'cat', get: p => fx(p, 'disulfide_pattern') },
            { label: 'β2 hairpin', part: 'core strands', tag: 'structure', kind: 'cat', get: p => fx(p, 'b2') },
            { label: 'Isoelectric point', part: 'catalytic domain', tag: 'sequence', get: p => fx(p, 'pi') },
            { label: 'N-glycosylation sequons (N-X-S/T)', part: 'catalytic domain', tag: 'sequence', unit: 'sites', get: p => fx(p, 'glycosites') },
            { label: 'Aromatic residues (F, W, Y, H; not the brace)', part: 'substrate-binding', tag: 'structure', unit: 'residues', get: p => fx(p, 'binding_aromatic') },
            { label: 'Asn + Gln', part: 'substrate-binding', tag: 'structure', unit: 'residues', get: p => fx(p, 'binding_amide') },
            { label: 'Charged residues (D, E, K, R)', part: 'substrate-binding', tag: 'structure', unit: 'residues', get: p => fx(p, 'binding_charged') },
            { label: 'Side chains H-bonded to the brace', part: 'brace', tag: 'structure', unit: 'residues', get: p => fx(p, 'brace_hbonds') },
            { label: 'Residues π-stacked on the brace', part: 'brace', tag: 'structure', unit: 'residues', get: p => fx(p, 'brace_stacking') },
            { label: 'Cation–π contacts to the brace', part: 'brace', tag: 'structure', unit: 'residues', get: p => fx(p, 'brace_cation_pi') },
            // the alignment columns the residues of the hole-hopping paths sit in, as a set: two
            // proteins differ by 1 - Jaccard (the share of their columns they do not share); with
            // the residue, a column counts as shared only when both have the same residue there
            { label: 'Hole-hopping positions', part: 'alignment columns of the path residues', tag: 'structure', kind: 'set',
              get: p => { const v = fx(p, 'hh_positions'); return v && v.map(x => x.split(':')[0]); } },
            { label: 'Hole-hopping positions and residues', part: 'alignment columns of the path residues', tag: 'structure', kind: 'set',
              get: p => fx(p, 'hh_positions') },
        );
        list.push(
            { label: 'Second-sphere motif', part: 'active site', tag: 'structure', kind: 'cat', get: p => cat((p.motif || [])[0]) },
            { label: 'Second-sphere motif with spacing', part: 'active site', tag: 'structure', kind: 'cat', get: p => cat((p.motif || [])[1] || (p.motif || [])[0]) },
            { label: 'Shielding residue', part: 'active site', tag: 'structure', kind: 'cat', get: p => cat(p.distal) },
            { label: 'N-His motif', part: 'brace', tag: 'sequence', kind: 'cat', get: p => cat(((p.brace || {}).n_his || [])[0]) },
            { label: 'Internal-His motif', part: 'brace', tag: 'sequence', kind: 'cat', get: p => cat(((p.brace || {}).internal_his || [])[0]) },
            // no active UniProt entry: unknown, not 'none'
            { label: 'CBM', part: 'UniProt', tag: 'annotation', kind: 'cat',
              get: p => ((p.cbms || {}).uniprot_known === false ? undefined : ((p.cbms || {}).uniprot || []).join(' + ') || 'none') },
            { label: 'Regioselectivity', part: 'papers or EC', tag: 'annotation', kind: 'cat', get: p => cat((p.regio || [])[0]) },
            { label: 'Taxon', part: 'NCBI', tag: 'annotation', kind: 'cat', get: p => (p.taxon || []).join(' · ') || undefined },
            // not biology, but it can make a grouping of its own
            { label: 'Experimental structure or model', part: 'data source', tag: 'annotation', kind: 'cat', caveat: true,
              get: p => (p.source === 'PDB' ? 'experimental' : 'model') },
        );
        // substrates reported in papers, one yes / no feature per substrate that at least
        // MIN_SUBSTRATE proteins of the catalogue are reported on; a protein without any reported
        // substrate has no value (not tested is not inactive), and "no" means not reported, which
        // need not mean inactive
        const reported = new Map();
        data.proteins.forEach(p => (p.substrate_papers || []).forEach(s => reported.set(s, (reported.get(s) || 0) + 1)));
        [...reported.entries()].filter(([, k]) => k >= MIN_SUBSTRATE).sort((a, b) => b[1] - a[1]).forEach(([s]) => {
            list.push({ label: `Substrate: ${s}`, part: 'reported in papers', tag: 'annotation', kind: 'cat',
                get: p => ((p.substrate_papers || []).length ? ((p.substrate_papers || []).includes(s) ? 'yes' : 'no') : undefined) });
        });
        return list.map((f, i) => ({ kind: 'num', ...f, id: i }));
    }

    // ---- scores ----

    function eta2(values, labels) {
        const n = values.length;
        if (n < 3) return null;
        const mean = values.reduce((a, b) => a + b, 0) / n;
        let sst = 0;
        values.forEach(v => { sst += (v - mean) ** 2; });
        if (sst < 1e-12) return null;
        const sums = new Map(), counts = new Map();
        values.forEach((v, i) => {
            sums.set(labels[i], (sums.get(labels[i]) || 0) + v);
            counts.set(labels[i], (counts.get(labels[i]) || 0) + 1);
        });
        let ssb = 0;
        sums.forEach((sum, g) => { ssb += counts.get(g) * (sum / counts.get(g) - mean) ** 2; });
        return ssb / sst;
    }

    function cramersV(cats, labels) {
        const n = cats.length;
        const rows = [...new Set(labels)], cols = [...new Set(cats)];
        if (n < 3 || rows.length < 2 || cols.length < 2) return null;
        const table = new Map();
        cats.forEach((c, i) => {
            const key = labels[i] + '\u0000' + c;
            table.set(key, (table.get(key) || 0) + 1);
        });
        const rowSum = new Map(), colSum = new Map();
        cats.forEach((c, i) => {
            rowSum.set(labels[i], (rowSum.get(labels[i]) || 0) + 1);
            colSum.set(c, (colSum.get(c) || 0) + 1);
        });
        let chi2 = 0;
        rows.forEach(r => cols.forEach(c => {
            const expected = rowSum.get(r) * colSum.get(c) / n;
            const observed = table.get(r + '\u0000' + c) || 0;
            chi2 += (observed - expected) ** 2 / expected;
        }));
        return Math.sqrt(chi2 / (n * (Math.min(rows.length, cols.length) - 1)));
    }

    // score of feature f for proteins `items` [{value}] under group labels `labels`
    function score(f, values, labels) {
        return f.kind === 'cat' ? cramersV(values, labels) : eta2(values, labels);
    }

    // Set features (kind 'set', values are arrays): two proteins differ by the Jaccard distance,
    // 1 - shared / all of their elements (0 for two empty sets); cached per feature by raw index
    function setDistance(f, i, j) {
        if (!f.sets) {
            f.sets = f.raw.map(v => (v === undefined ? undefined : new Set(v)));
            f.distances = new Map();
        }
        const key = i < j ? i * 65536 + j : j * 65536 + i;
        let d = f.distances.get(key);
        if (d === undefined) {
            const a = f.sets[i], b = f.sets[j];
            let shared = 0;
            a.forEach(x => { if (b.has(x)) shared++; });
            const all = a.size + b.size - shared;
            d = all ? 1 - shared / all : 0;
            f.distances.set(key, d);
        }
        return d;
    }

    // how different proteins i and j are in feature f (raw indices); undefined without both values
    function pairDiff(f, i, j) {
        const a = f.raw[i], b = f.raw[j];
        if (a === undefined || b === undefined) return undefined;
        if (f.kind === 'num') return Math.abs(a - b);
        if (f.kind === 'cat') return a === b ? 0 : 1;
        return setDistance(f, i, j);
    }

    // PERMANOVA R² of a set feature for groups `lab` over the proteins `idx`: 1 - the within-group
    // sum of squared distances / the total one (each sum over pairs, divided by the number of
    // proteins it covers), the same quantity as eta squared for distances
    function permanova(f, idx, lab) {
        const n = idx.length;
        if (n < 3) return null;
        const size = new Map(), within = new Map();
        idx.forEach(i => size.set(lab[i], (size.get(lab[i]) || 0) + 1));
        let total = 0;
        for (let a = 0; a < n; a++) {
            for (let b = a + 1; b < n; b++) {
                const d2 = setDistance(f, idx[a], idx[b]) ** 2;
                total += d2;
                if (lab[idx[a]] === lab[idx[b]]) within.set(lab[idx[a]], (within.get(lab[idx[a]]) || 0) + d2);
            }
        }
        total /= n;
        if (total < 1e-12) return null;
        let w = 0;
        within.forEach((sum, g) => { w += sum / size.get(g); });
        return 1 - w / total;
    }

    // a set feature's elements by how many of `items` (arrays) hold them, most common first
    function setShares(items) {
        const counts = new Map();
        items.forEach(v => v.forEach(x => counts.set(x, (counts.get(x) || 0) + 1)));
        return [...counts.entries()].sort((a, b) => b[1] - a[1] || String(a[0]).localeCompare(String(b[0]), undefined, { numeric: true }));
    }
    // '428:W' -> 'column 428 W', '428' -> 'column 428'
    const setElement = x => { const [c, aa] = String(x).split(':'); return `col ${c}${aa ? ' ' + aa : ''}`; };

    // how strongly two features move together across the same proteins (0..1)
    function association(fa, va, fb, vb) {
        if (fa.kind === 'set' || fb.kind === 'set') {
            // over the pairs of proteins: how different they are in one feature against the other
            const xs = [], ys = [];
            for (let i = 0; i < va.length; i++) {
                for (let j = i + 1; j < va.length; j++) {
                    const x = pairDiff(fa, i, j), y = pairDiff(fb, i, j);
                    if (x !== undefined && y !== undefined) { xs.push(x); ys.push(y); }
                }
            }
            return xs.length >= 6 ? Math.abs(pearson(xs, ys) || 0) : 0;
        }
        const both = va.map((a, i) => [a, vb[i]]).filter(([a, b]) => a !== undefined && b !== undefined);
        if (both.length < 4) return 0;
        const xs = both.map(p => p[0]), ys = both.map(p => p[1]);
        if (fa.kind === 'num' && fb.kind === 'num') {
            const mx = xs.reduce((a, b) => a + b, 0) / xs.length, my = ys.reduce((a, b) => a + b, 0) / ys.length;
            let sxy = 0, sxx = 0, syy = 0;
            xs.forEach((x, i) => { sxy += (x - mx) * (ys[i] - my); sxx += (x - mx) ** 2; syy += (ys[i] - my) ** 2; });
            return sxx && syy ? Math.abs(sxy / Math.sqrt(sxx * syy)) : 0;
        }
        if (fa.kind === 'cat' && fb.kind === 'cat') return cramersV(xs, ys) || 0;
        // a number against a category: the correlation ratio
        const [nums, cats] = fa.kind === 'num' ? [xs, ys] : [ys, xs];
        return Math.sqrt(eta2(nums, cats) || 0);
    }

    // a small seeded generator, so pressing the button twice gives the same chance level
    function mulberry32(seed) {
        return function () {
            seed |= 0; seed = seed + 0x6D2B79F5 | 0;
            let t = Math.imul(seed ^ seed >>> 15, 1 | seed);
            t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
            return ((t ^ t >>> 14) >>> 0) / 4294967296;
        };
    }

    // every feature with its values for `codes`, in that order (undefined = no value)
    function featureValues(data, codes, byCode) {
        const proteins = codes.map(c => byCode.get(c));
        return featureList(data).map(f => {
            const raw = proteins.map(p => {
                const v = p ? f.get(p) : undefined;
                if (f.kind === 'num') return typeof v === 'number' && isFinite(v) ? v : undefined;
                if (f.kind === 'set') return Array.isArray(v) ? v : undefined;
                return v;
            });
            return { ...f, raw };
        });
    }

    // hang each ranked feature under the first higher-ranked one it moves with
    function nest(scored) {
        const leaders = [];
        scored.forEach(f => {
            const leader = leaders.find(l => association(l, l.raw, f, f.raw) >= REDUNDANT);
            if (leader) {
                f.follows = leader;
                leader.followers.push(f);
            } else {
                f.followers = [];
                leaders.push(f);
            }
        });
        return leaders;
    }

    // ---- "Which features follow the tree?": no cut, every pair of proteins ----

    // average ranks (ties share their mean rank), for Spearman's correlation
    function ranks(values) {
        const order = values.map((v, i) => [v, i]).sort((a, b) => a[0] - b[0]);
        const out = new Array(values.length);
        for (let i = 0; i < order.length;) {
            let j = i;
            while (j + 1 < order.length && order[j + 1][0] === order[i][0]) j++;
            const mean = (i + j) / 2 + 1;
            for (let k = i; k <= j; k++) out[order[k][1]] = mean;
            i = j + 1;
        }
        return out;
    }

    function pearson(xs, ys) {
        const n = xs.length;
        const mx = xs.reduce((a, b) => a + b, 0) / n, my = ys.reduce((a, b) => a + b, 0) / n;
        let sxy = 0, sxx = 0, syy = 0;
        for (let i = 0; i < n; i++) {
            sxy += (xs[i] - mx) * (ys[i] - my);
            sxx += (xs[i] - mx) ** 2;
            syy += (ys[i] - my) ** 2;
        }
        return sxx && syy ? sxy / Math.sqrt(sxx * syy) : null;
    }

    // A Mantel test per feature: over every pair of proteins, Spearman's correlation between how
    // far apart the two are in the tree and how different they are in the feature (the absolute
    // difference of a number, 0 / 1 for the same or another category). Positive = proteins far
    // apart in the tree also differ more in the feature. The tree distances are permuted by
    // protein (rows and columns together, so each protein's pairs stay one unit) for the chance
    // level, the 95th percentile of the best correlation over all features, as for the groups.
    function analyseTree(data, cophenetic) {
        const byCode = new Map(data.proteins.map(p => [p.code, p]));
        const codes = cophenetic.codes;
        const n = codes.length;
        if (n < 4) return { mode: 'tree', error: 'Needs at least four compared proteins in the tree.' };
        const pairs = [];
        for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) pairs.push([i, j]);
        const tree = pairs.map(([i, j]) => cophenetic.matrix[i][j]);
        const features = featureValues(data, codes, byCode);

        features.forEach(f => {
            f.n = f.raw.filter(v => v !== undefined).length;
            // the pairs where both proteins have a value, and how different they are
            f.pairIdx = [];
            f.diff = [];
            pairs.forEach(([i, j], k) => {
                const a = f.raw[i], b = f.raw[j];
                if (a === undefined || b === undefined) return;
                f.pairIdx.push(k);
                f.diff.push(pairDiff(f, i, j));
            });
            f.diffRanks = ranks(f.diff);
            const r = f.pairIdx.length >= 6 && new Set(f.diff).size > 1
                ? pearson(ranks(f.pairIdx.map(k => tree[k])), f.diffRanks) : null;
            f.score = r;
        });
        const scored = features.filter(f => f.score != null);

        const pairIndex = (i, j) => { if (i > j) [i, j] = [j, i]; return i * n - i * (i + 1) / 2 + (j - i - 1); };
        const random = mulberry32(20260925);
        const order = codes.map((_, i) => i);
        const maxima = [];
        const exceed = new Map(scored.map(f => [f.id, 0]));
        for (let k = 0; k < PERMUTATIONS; k++) {
            for (let i = order.length - 1; i > 0; i--) {
                const j = Math.floor(random() * (i + 1));
                [order[i], order[j]] = [order[j], order[i]];
            }
            const permuted = pairs.map(([i, j]) => tree[pairIndex(order[i], order[j])]);
            let best = -1;
            scored.forEach(f => {
                const r = pearson(ranks(f.pairIdx.map(q => permuted[q])), f.diffRanks);
                if (r == null) return;
                if (r > best) best = r;
                if (r >= f.score - 1e-12) exceed.set(f.id, exceed.get(f.id) + 1);
            });
            maxima.push(best);
        }
        maxima.sort((a, b) => a - b);
        const chance = Math.max(0, maxima[Math.floor(0.95 * (maxima.length - 1))]);
        scored.forEach(f => {
            f.p = (exceed.get(f.id) + 1) / (PERMUTATIONS + 1);
            f.beats = f.score > chance + 1e-9;
        });
        scored.sort((a, b) => b.score - a.score || a.id - b.id);
        return { mode: 'tree', codes, pairs, tree, leaders: nest(scored), chance, n };
    }

    function analyse(data, grouping) {
        const byCode = new Map(data.proteins.map(p => [p.code, p]));
        const kept = grouping.groups;
        const singletons = kept.filter(g => g.codes.length < 2).flatMap(g => g.codes);
        if (kept.length < 2 || kept.every(g => g.codes.length < 2)) {
            return {
                error: grouping.user
                    ? (kept.length < 2
                        ? 'Only one group: add a second one on the browse page ("Add as group"), or compare it with proteins left ungrouped.'
                        : 'Every group is a single protein, so any feature would separate them: add more proteins to the groups.')
                    : kept.length < 2
                    ? 'Only one group: cut at a smaller distance (the line further right in the tree).'
                    : 'Every group is a single protein, so any feature would separate them: cut at a larger distance (the line further left in the tree).',
                singletons,
            };
        }
        const codes = kept.flatMap(g => g.codes);
        const labels = kept.flatMap((g, i) => g.codes.map(() => i));
        const features = featureValues(data, codes, byCode);

        // a feature scores on the proteins that have a value, when they fall in at least two
        // groups and at least one of those groups holds two of them
        function scoreOn(f, lab) {
            const idx = f.raw.map((v, i) => (v === undefined ? -1 : i)).filter(i => i >= 0);
            const perGroup = new Map();
            idx.forEach(i => perGroup.set(lab[i], (perGroup.get(lab[i]) || 0) + 1));
            const sizes = [...perGroup.values()];
            if (sizes.length < 2 || !sizes.some(k => k >= 2)) return null;
            if (f.kind === 'set') return permanova(f, idx, lab);
            return score(f, idx.map(i => f.raw[i]), idx.map(i => lab[i]));
        }
        features.forEach(f => { f.score = scoreOn(f, labels); f.n = f.raw.filter(v => v !== undefined).length; });
        const scored = features.filter(f => f.score != null);

        // chance level: the best score over all features on shuffled groups
        const random = mulberry32(20260925);
        const shuffled = labels.slice();
        const maxima = [];
        const exceed = new Map(scored.map(f => [f.id, 0]));
        for (let k = 0; k < PERMUTATIONS; k++) {
            for (let i = shuffled.length - 1; i > 0; i--) {
                const j = Math.floor(random() * (i + 1));
                [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]];
            }
            let best = 0;
            scored.forEach(f => {
                const s = scoreOn(f, shuffled);
                if (s == null) return;
                if (s > best) best = s;
                if (s >= f.score - 1e-12) exceed.set(f.id, exceed.get(f.id) + 1);
            });
            maxima.push(best);
        }
        maxima.sort((a, b) => a - b);
        const chance = maxima[Math.floor(0.95 * (maxima.length - 1))];
        scored.forEach(f => { f.p = (exceed.get(f.id) + 1) / (PERMUTATIONS + 1); });

        scored.sort((a, b) => b.score - a.score || a.id - b.id);
        const leaders = nest(scored);
        scored.forEach(f => { f.beats = f.score > chance + 1e-9; });
        return { mode: 'groups', groups: kept, singletons, codes, labels, leaders, chance, cut: grouping.cut,
                 user: !!grouping.user, n: codes.length };
    }

    // ---- "Which features discriminate the picked proteins from …": against a reference ----

    const REFERENCES = {
        comparison: 'the rest of the comparison',
        families: 'the rest of the catalogue in the same families',
        catalogue: 'the rest of the catalogue',
    };

    const median = xs => {
        const s = xs.slice().sort((a, b) => a - b);
        const k = (s.length - 1) / 2;
        return (s[Math.floor(k)] + s[Math.ceil(k)]) / 2;
    };

    function analysePicked(data, picked, compared, reference) {
        const byCode = new Map(data.proteins.map(p => [p.code, p]));
        const whole = !picked.length;
        const A = whole ? compared.filter(c => byCode.has(c)) : picked.filter(c => byCode.has(c));
        const inA = new Set(A);
        if (whole && reference === 'comparison') {
            return { error: 'Pick proteins in the tree (a row, or a branch for its subtree) to compare them with the rest of the comparison.' };
        }
        let B;
        if (reference === 'comparison') {
            B = compared.filter(c => byCode.has(c) && !inA.has(c));
        } else {
            B = data.proteins.map(p => p.code).filter(c => !inA.has(c));
            if (reference === 'families') {
                const families = new Set(A.map(c => byCode.get(c).family));
                B = B.filter(c => families.has(byCode.get(c).family));
            }
        }
        if (B.length < 2) {
            return { error: `Too few proteins to compare with: ${REFERENCES[reference]} holds ${B.length}.` };
        }
        const codes = A.concat(B);
        const labels = codes.map((_, i) => (i < A.length ? 1 : 0));
        const features = featureValues(data, codes, byCode);

        // per feature, the proteins with a value; numbers ranked once (ties share their rank),
        // categories as indices, so a shuffle only re-adds the picked proteins' part
        features.forEach(f => {
            f.idx = f.raw.map((v, i) => (v === undefined ? -1 : i)).filter(i => i >= 0);
            f.n = f.idx.length;
            if (f.kind === 'num') {
                f.rank = new Array(codes.length);
                const rk = ranks(f.idx.map(i => f.raw[i]));
                f.idx.forEach((i, k) => { f.rank[i] = rk[k]; });
            } else {
                f.cats = [...new Set(f.idx.map(i => f.raw[i]))];
                const at = new Map(f.cats.map((c, k) => [c, k]));
                f.catOf = new Array(codes.length);
                f.idx.forEach(i => { f.catOf[i] = at.get(f.raw[i]); });
                f.catTotal = f.cats.map(() => 0);
                f.idx.forEach(i => { f.catTotal[f.catOf[i]]++; });
            }
        });
        // score of feature f when `lab` marks the picked proteins; the signed AUC for numbers. A
        // feature counts only when at least half of the picked proteins have a value: otherwise
        // a lone value among them (a loop most proteins lack, a substrate most papers do not
        // report) easily looks extreme, and such features set the chance level
        const needed = Math.max(1, Math.ceil(A.length / 2));
        function scoreOn(f, lab) {
            let nA = 0;
            if (f.kind === 'set') {
                // ANOSIM-like: are the picked proteins closer to each other than to the reference?
                // |2 AUC - 1| over the within-pick distances against the pick-to-reference ones,
                // positive = more alike among themselves; needs two picked proteins with a value
                const a = f.idx.filter(i => lab[i]), b = f.idx.filter(i => !lab[i]);
                if (a.length < Math.max(2, needed) || b.length < 2) return null;
                const within = [], across = [];
                for (let x = 0; x < a.length; x++) {
                    for (let y = x + 1; y < a.length; y++) within.push(setDistance(f, a[x], a[y]));
                    for (let y = 0; y < b.length; y++) across.push(setDistance(f, a[x], b[y]));
                }
                const rk = ranks(within.concat(across));
                let sum = 0;
                for (let k = 0; k < within.length; k++) sum += rk[k];
                const higher = (sum - within.length * (within.length + 1) / 2) / (within.length * across.length);
                return 1 - 2 * higher;
            }
            if (f.kind === 'num') {
                let sum = 0;
                f.idx.forEach(i => { if (lab[i]) { nA++; sum += f.rank[i]; } });
                const nB = f.n - nA;
                if (nA < needed || nB < 2) return null;
                const auc = (sum - nA * (nA + 1) / 2) / (nA * nB);
                return 2 * auc - 1;
            }
            const inPicked = f.cats.map(() => 0);
            f.idx.forEach(i => { if (lab[i]) { nA++; inPicked[f.catOf[i]]++; } });
            const nB = f.n - nA;
            if (nA < needed || nB < 2 || f.cats.length < 2) return null;
            let tvd = 0;
            f.cats.forEach((_, k) => { tvd += Math.abs(inPicked[k] / nA - (f.catTotal[k] - inPicked[k]) / nB); });
            return tvd / 2;
        }
        features.forEach(f => {
            const s = scoreOn(f, labels);
            f.signed = s;
            f.score = s == null ? null : Math.abs(s);
        });
        const scored = features.filter(f => f.score != null && (f.kind !== 'num' || new Set(f.idx.map(i => f.raw[i])).size > 1));

        const random = mulberry32(20260925);
        const shuffled = labels.slice();
        const maxima = [];
        const exceed = new Map(scored.map(f => [f.id, 0]));
        for (let k = 0; k < PERMUTATIONS; k++) {
            for (let i = shuffled.length - 1; i > 0; i--) {
                const j = Math.floor(random() * (i + 1));
                [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]];
            }
            let best = 0;
            scored.forEach(f => {
                const s = scoreOn(f, shuffled);
                if (s == null) return;
                const a = Math.abs(s);
                if (a > best) best = a;
                if (a >= f.score - 1e-12) exceed.set(f.id, exceed.get(f.id) + 1);
            });
            maxima.push(best);
        }
        maxima.sort((a, b) => a - b);
        const chance = maxima[Math.floor(0.95 * (maxima.length - 1))];

        // what the score is about, for the list and the tooltip
        scored.forEach(f => {
            f.p = (exceed.get(f.id) + 1) / (PERMUTATIONS + 1);
            f.beats = f.score > chance + 1e-9;
            const a = f.idx.filter(i => labels[i]).map(i => f.raw[i]);
            const b = f.idx.filter(i => !labels[i]).map(i => f.raw[i]);
            f.valuesA = a;
            f.valuesB = b;
            if (f.kind === 'set') {
                // the mean overlap (Jaccard similarity) among the picked proteins and with the reference
                const ia = f.idx.filter(i => labels[i]), ib = f.idx.filter(i => !labels[i]);
                let w = 0, nw = 0, x = 0, nx = 0;
                ia.forEach((i, k) => {
                    ia.slice(k + 1).forEach(j => { w += 1 - setDistance(f, i, j); nw++; });
                    ib.forEach(j => { x += 1 - setDistance(f, i, j); nx++; });
                });
                f.overlap = [w / nw, x / nx];
                f.direction = `${f.signed > 0 ? 'more alike' : 'less alike'}: overlap ${fmt(w / nw)} vs ${fmt(x / nx)}`;
            } else if (f.kind === 'num') {
                // the medians, or the means when the medians are the same (whole numbers often tie)
                const mean = xs => xs.reduce((x, y) => x + y, 0) / xs.length;
                const same = fmt(median(a)) === fmt(median(b));
                const [va, vb] = same ? [mean(a), mean(b)] : [median(a), median(b)];
                f.direction = `${f.signed > 0 ? 'higher' : 'lower'}: ${same ? 'mean ' : ''}${fmt(va)} vs ${fmt(vb)}`;
            } else {
                // the category whose share differs most between the two sets
                const share = (xs, c) => xs.filter(x => x === c).length / xs.length;
                const top = f.cats.map(c => [c, share(a, c) - share(b, c)])
                    .sort((x, y) => Math.abs(y[1]) - Math.abs(x[1]))[0][0];
                f.direction = `${top}: ${a.filter(x => x === top).length}/${a.length} vs `
                    + `${Math.round(100 * share(b, top))} %`;
            }
        });
        scored.sort((a, b) => b.score - a.score || a.id - b.id);
        return {
            mode: 'picked', A, B, codes, labels, leaders: nest(scored), chance, n: codes.length,
            whole, reference, single: A.length === 1,
        };
    }

    // ---- drawing ----

    function fmt(v) {
        return Math.abs(v) >= 100 ? v.toFixed(0) : Math.abs(v) >= 10 ? v.toFixed(1) : v.toFixed(2);
    }

    // tree mode: one point per pair (tree distance, how different the two are) and a summary of
    // the closest against the farthest half of the pairs
    function treeDetail(result, f) {
        const box = el('div', 'gf-detail');
        box.appendChild(el('div', 'gf-detail-head', `${f.label} · ${f.part}` + (f.unit ? ` (${f.unit})` : '')));
        box.appendChild(el('div', 'gf-detail-meta',
            `Spearman r ${f.score >= 0 ? '+' : ''}${f.score.toFixed(2)} · `
            + `p ${f.p < 0.001 ? '< 0.001' : f.p.toFixed(3)} (tree shuffled by protein, this feature alone) · `
            + `${f.n} of ${result.n} proteins, ${f.pairIdx.length} pairs`));
        const dist = f.pairIdx.map(k => result.tree[k]);
        const sorted = dist.slice().sort((a, b) => a - b);
        const median = sorted[Math.floor((sorted.length - 1) / 2)];
        const near = f.diff.filter((_, i) => dist[i] <= median), far = f.diff.filter((_, i) => dist[i] > median);
        const mean = xs => (xs.length ? xs.reduce((a, b) => a + b, 0) / xs.length : NaN);
        box.appendChild(el('div', 'gf-detail-meta', f.kind === 'set'
            ? `Mean Jaccard distance ${fmt(mean(near))} in the closer half of the pairs, `
              + `${fmt(mean(far))} in the farther half (0 = the same positions, 1 = none shared).`
            : f.kind === 'cat'
            ? `Differ in ${Math.round(100 * mean(near))} % of the closer half of the pairs, `
              + `${Math.round(100 * mean(far))} % of the farther half.`
            : `Mean difference ${fmt(mean(near))} in the closer half of the pairs, `
              + `${fmt(mean(far))} in the farther half` + (f.unit ? ` (${f.unit}).` : '.')));
        if (f.kind !== 'cat') {
            const ns = 'http://www.w3.org/2000/svg', W = 250, H = 110, L = 30, B = 18;
            const svg = document.createElementNS(ns, 'svg');
            svg.setAttribute('width', W + L + 6);
            svg.setAttribute('height', H + B + 4);
            const put = (tag, attrs, text) => {
                const node = document.createElementNS(ns, tag);
                Object.entries(attrs).forEach(([k, v]) => node.setAttribute(k, v));
                if (text != null) node.textContent = text;
                svg.appendChild(node);
                return node;
            };
            const dx = Math.max(...dist) || 1, dy = Math.max(...f.diff) || 1;
            put('line', { x1: L, x2: L + W, y1: H, y2: H, stroke: '#cfd8d3' });
            put('line', { x1: L, x2: L, y1: 2, y2: H, stroke: '#cfd8d3' });
            put('text', { x: L + W, y: H + 13, 'text-anchor': 'end', 'font-size': 9, fill: '#7a8a82' }, 'tree distance →');
            put('text', { x: L - 4, y: 9, 'text-anchor': 'end', 'font-size': 9, fill: '#7a8a82' }, fmt(dy));
            put('text', { x: L - 4, y: H, 'text-anchor': 'end', 'font-size': 9, fill: '#7a8a82' }, '0');
            f.pairIdx.forEach((k, i) => {
                const [a, b] = result.pairs[k];
                const dot = put('circle', { cx: L + dist[i] / dx * W, cy: H - f.diff[i] / dy * (H - 6), r: 2.6,
                    fill: TAGS[f.tag].color, 'fill-opacity': 0.55 });
                const title = document.createElementNS(ns, 'title');
                title.textContent = `${result.codes[a]} – ${result.codes[b]}: difference ${fmt(f.diff[i])}`;
                dot.appendChild(title);
            });
            box.appendChild(svg);
        }
        return box;
    }

    // per group: a row of dots (numbers) or the category counts, in the group's shade
    // picked mode: the picked proteins (amber) above the reference (grey)
    const PICKED_COLOR = '#d9a400', REFERENCE_COLOR = '#9aa5a0';
    function pickedDetail(result, f) {
        const box = el('div', 'gf-detail');
        box.appendChild(el('div', 'gf-detail-head', `${f.label} · ${f.part}` + (f.unit ? ` (${f.unit})` : '')));
        box.appendChild(el('div', 'gf-detail-meta',
            (f.kind === 'set'
                ? `|2 AUC − 1| ${f.score.toFixed(2)}: two picked proteins are ${f.signed > 0 ? 'closer to each other' : 'farther apart'} than a picked and a reference protein in `
                  + `${Math.round(50 * (1 + Math.abs(f.signed)))} % of the comparisons (Jaccard distance; mean overlap `
                  + `${fmt(f.overlap[0])} among the picked, ${fmt(f.overlap[1])} with the reference)`
                : f.kind === 'num'
                ? `|2 AUC − 1| ${f.score.toFixed(2)}: the picked proteins are ${f.signed > 0 ? 'higher' : 'lower'} than the reference in `
                  + `${Math.round(50 * (1 + Math.abs(f.signed)))} % of the pairs`
                : `Total variation distance ${f.score.toFixed(2)} between the category shares`)
            + ` · p ${f.p < 0.001 ? '< 0.001' : f.p.toFixed(3)} (shuffled picks, this feature alone) · `
            + `${f.valuesA.length} picked and ${f.valuesB.length} reference proteins with a value`));
        const sets = [
            ['Picked', PICKED_COLOR, result.codes.map((c, i) => [c, f.raw[i], result.labels[i]]).filter(x => x[2] && x[1] !== undefined)],
            ['Reference', REFERENCE_COLOR, result.codes.map((c, i) => [c, f.raw[i], result.labels[i]]).filter(x => !x[2] && x[1] !== undefined)],
        ];
        if (f.kind === 'set') {
            // the positions most of the picked proteins hold, with the reference's share of each
            const picked = sets[0][2].map(x => x[1]), reference = sets[1][2].map(x => x[1]);
            const refCount = new Map(setShares(reference));
            const row = (color, text) => {
                const r = el('div', 'gf-strip');
                const sw = el('span', 'gf-swatch');
                sw.style.background = color;
                r.appendChild(sw);
                r.appendChild(el('span', 'gf-strip-codes', text));
                box.appendChild(r);
            };
            row(PICKED_COLOR, `Picked (${picked.length}): ` + (setShares(picked).slice(0, 8)
                .map(([x, k]) => `${setElement(x)} ${k}/${picked.length} (reference ${Math.round(100 * (refCount.get(x) || 0) / reference.length)} %)`)
                .join(', ') || '–'));
            row(REFERENCE_COLOR, `Reference (${reference.length}): ` + (setShares(reference).slice(0, 6)
                .map(([x, k]) => `${setElement(x)} ${Math.round(100 * k / reference.length)} %`).join(', ') || '–'));
        } else if (f.kind === 'num') {
            const all = sets.flatMap(s => s[2].map(x => x[1]));
            const lo = Math.min(...all), hi = Math.max(...all), span = hi - lo || 1;
            const W = 230, ns = 'http://www.w3.org/2000/svg';
            sets.forEach(([name, color, items], k) => {
                const row = el('div', 'gf-strip');
                const sw = el('span', 'gf-swatch');
                sw.style.background = color;
                row.appendChild(sw);
                const svg = document.createElementNS(ns, 'svg');
                svg.setAttribute('width', W + 12);
                svg.setAttribute('height', 14);
                items.forEach(([code, v]) => {
                    const c = document.createElementNS(ns, 'circle');
                    Object.entries({ cx: 6 + (v - lo) / span * W, cy: 7, r: k ? 3 : 4, fill: color,
                        'fill-opacity': k ? 0.45 : 0.9, stroke: 'white', 'stroke-width': 0.6 })
                        .forEach(([a, val]) => c.setAttribute(a, val));
                    svg.appendChild(c);
                });
                row.appendChild(svg);
                const values = items.map(x => x[1]);
                row.appendChild(el('span', 'gf-strip-codes', k
                    ? `${name}: ${items.length}, median ${fmt(median(values))}`
                    : items.length <= 8 ? items.map(([c, v]) => `${c} ${fmt(v)}`).join(', ')
                    : `${name}: ${items.length}, median ${fmt(median(values))}`));
                box.appendChild(row);
            });
            box.appendChild(el('div', 'gf-axis', `${fmt(lo)} … ${fmt(hi)}`));
        } else {
            sets.forEach(([name, color, items]) => {
                const counts = new Map();
                items.forEach(([, v]) => counts.set(v, (counts.get(v) || 0) + 1));
                const row = el('div', 'gf-strip');
                const sw = el('span', 'gf-swatch');
                sw.style.background = color;
                row.appendChild(sw);
                const top = [...counts.entries()].sort((a, b) => b[1] - a[1]);
                const text = top.slice(0, 6).map(([v, k]) => `${v} ${Math.round(100 * k / items.length)} %`).join(', ')
                    + (top.length > 6 ? ', …' : '');
                row.appendChild(el('span', 'gf-strip-codes', `${name} (${items.length}): ${text || '–'}`));
                box.appendChild(row);
            });
        }
        return box;
    }

    function detail(result, f) {
        if (result.mode === 'tree') return treeDetail(result, f);
        if (result.mode === 'picked') return pickedDetail(result, f);
        const box = el('div', 'gf-detail');
        const head = el('div', 'gf-detail-head');
        head.textContent = `${f.label} · ${f.part}` + (f.unit ? ` (${f.unit})` : '');
        box.appendChild(head);
        const meta = el('div', 'gf-detail-meta');
        meta.textContent = `${f.kind === 'cat' ? "Cramér's V" : f.kind === 'set' ? 'PERMANOVA R² (Jaccard distance)' : 'η²'} ${f.score.toFixed(2)} · `
            + `p ${f.p < 0.001 ? '< 0.001' : f.p.toFixed(3)} (shuffled groups, this feature alone) · `
            + `${f.n} of ${result.n} proteins with a value`;
        box.appendChild(meta);
        const values = result.groups.map((g, gi) => result.codes
            .map((code, i) => [code, f.raw[i], result.labels[i]])
            .filter(([, v, lab]) => lab === gi && v !== undefined));
        if (f.kind === 'set') {
            // per group, the positions most of its proteins hold
            result.groups.forEach((g, gi) => {
                const items = values[gi].map(x => x[1]);
                const row = el('div', 'gf-strip');
                const sw = el('span', 'gf-swatch');
                sw.style.background = g.color;
                row.appendChild(sw);
                row.appendChild(el('span', 'gf-strip-codes', setShares(items).slice(0, 6)
                    .map(([x, k]) => `${setElement(x)} ${k}/${items.length}`).join(', ') || '–'));
                box.appendChild(row);
            });
        } else if (f.kind === 'num') {
            const all = values.flat().map(x => x[1]);
            const lo = Math.min(...all), hi = Math.max(...all), span = hi - lo || 1;
            const W = 230, ns = 'http://www.w3.org/2000/svg';
            result.groups.forEach((g, gi) => {
                const row = el('div', 'gf-strip');
                const sw = el('span', 'gf-swatch');
                sw.style.background = g.color;
                row.appendChild(sw);
                const svg = document.createElementNS(ns, 'svg');
                svg.setAttribute('width', W + 12);
                svg.setAttribute('height', 14);
                const line = document.createElementNS(ns, 'line');
                Object.entries({ x1: 6, x2: W + 6, y1: 7, y2: 7, stroke: '#e3e8e5' }).forEach(([k, v]) => line.setAttribute(k, v));
                svg.appendChild(line);
                values[gi].forEach(([code, v]) => {
                    const c = document.createElementNS(ns, 'circle');
                    Object.entries({ cx: 6 + (v - lo) / span * W, cy: 7, r: 4, fill: g.color, 'fill-opacity': 0.85,
                        stroke: 'white', 'stroke-width': 0.8 }).forEach(([k, val]) => c.setAttribute(k, val));
                    const t = document.createElementNS(ns, 'title');
                    t.textContent = `${code}: ${fmt(v)}`;
                    c.appendChild(t);
                    svg.appendChild(c);
                });
                row.appendChild(svg);
                row.appendChild(el('span', 'gf-strip-codes', values[gi].map(([c, v]) => `${c} ${fmt(v)}`).join(', ')));
                box.appendChild(row);
            });
            box.appendChild(el('div', 'gf-axis', `${fmt(lo)} … ${fmt(hi)}`));
        } else {
            result.groups.forEach((g, gi) => {
                const counts = new Map();
                values[gi].forEach(([, v]) => counts.set(v, (counts.get(v) || 0) + 1));
                const row = el('div', 'gf-strip');
                const sw = el('span', 'gf-swatch');
                sw.style.background = g.color;
                row.appendChild(sw);
                row.appendChild(el('span', 'gf-strip-codes',
                    [...counts.entries()].sort((a, b) => b[1] - a[1]).map(([v, k]) => `${v} ×${k}`).join(', ') || '–'));
                box.appendChild(row);
            });
        }
        return box;
    }

    function render(root, result, showAll) {
        root.replaceChildren();
        if (result.error) {
            root.appendChild(el('p', 'ov-hint', result.error));
            return;
        }
        const tree = result.mode === 'tree';
        const picked = result.mode === 'picked';
        const head = el('div', 'gf-head');
        const pickedNames = result.A && (result.A.length <= 6 ? result.A.join(', ') : `${result.A.length} proteins`);
        head.textContent = tree
            ? `Following the tree: all ${result.n} proteins, ${result.pairs.length} pairs, no cut`
            : picked
            ? `${result.whole ? `All ${result.A.length} compared proteins` : `Picked: ${pickedNames}`}`
              + ` from ${REFERENCES[result.reference]} (${result.B.length})`
            : result.user
            ? `Your groups: ${result.groups.map(g => `${g.name} (${g.codes.length})`).join(' · ')}`
            : `At ${result.cut.toFixed(2)}: ${result.groups.length} groups, ${result.n} proteins`
              + (result.singletons.length ? ` · on their own: ${result.singletons.join(', ')}` : '');
        root.appendChild(head);
        const legend = el('div', 'gf-legend');
        Object.values(TAGS).forEach(t => {
            const item = el('span', 'gf-legend-item');
            const sw = el('span', 'gf-tag');
            sw.style.background = t.color;
            item.appendChild(sw);
            item.appendChild(document.createTextNode(t.label));
            legend.appendChild(item);
        });
        const chance = el('span', 'gf-legend-item');
        chance.appendChild(el('span', 'gf-chance-key'));
        chance.appendChild(document.createTextNode(`chance level ${result.chance.toFixed(2)}`));
        legend.appendChild(chance);
        root.appendChild(legend);
        // With many small groups random labels already separate the proteins perfectly on some
        // feature, so nothing can clear the bar: say so instead of listing perfect-looking scores.
        if (picked && result.single) {
            root.appendChild(el('p', 'gf-warning',
                'One protein: this describes how it differs from the reference, not a finding. With a '
                + 'single value many features look extreme, which is why the chance level is high.'));
        } else if (picked && result.chance >= 0.99) {
            root.appendChild(el('p', 'gf-warning',
                `Random picks of ${result.A.length} proteins already reach ${result.chance.toFixed(2)}, so no `
                + 'feature can be told apart from chance: pick more proteins or compare with a larger set.'));
        }
        if (!tree && !picked && result.chance >= 0.99) {
            root.appendChild(el('p', 'gf-warning',
                `${result.user ? 'With groups this small' : 'At this cut'} random groupings already reach ${result.chance.toFixed(2)}: with `
                + `${result.groups.length} groups for ${result.n} proteins almost any feature separates `
                + 'them, so no feature can be told apart from chance. '
                + (result.user ? 'Add more proteins to the groups (on the browse page).'
                               : 'Cut at a larger distance (the line further left in the tree).')));
        }

        const list = el('div', 'gf-list' + (picked ? ' gf-picked' : ''));
        const tooltip = root.__tooltip;
        const rows = [];
        result.leaders.slice(0, showAll ? undefined : SHOWN).forEach(f => {
            rows.push([f, false]);
            f.followers.forEach(x => rows.push([x, true]));
        });
        rows.forEach(([f, follower]) => {
            const row = el('div', 'gf-row' + (follower ? ' gf-follower' : '') + (f.beats ? '' : ' gf-below'));
            const tag = el('span', 'gf-tag');
            tag.style.background = TAGS[f.tag].color;
            tag.title = TAGS[f.tag].label;
            row.appendChild(tag);
            const name = el('span', 'gf-name');
            name.textContent = (follower ? '↳ ' : '') + f.label;
            name.appendChild(el('span', 'gf-part', ` · ${f.part}`));
            if (f.caveat) name.appendChild(el('span', 'gf-caveat', ' ⚠'));
            row.appendChild(name);
            if (picked) row.appendChild(el('span', 'gf-dir', f.direction));
            const bar = el('span', 'gf-bar');
            const fill = el('span', 'gf-fill');
            fill.style.width = `${(Math.max(0, f.score) * 100).toFixed(1)}%`;   // r < 0 draws no bar
            fill.style.background = TAGS[f.tag].color;
            bar.appendChild(fill);
            const mark = el('span', 'gf-chance');
            mark.style.left = `${(result.chance * 100).toFixed(1)}%`;
            bar.appendChild(mark);
            row.appendChild(bar);
            row.appendChild(el('span', 'gf-score', (tree && f.score > 0 ? '+' : '') + f.score.toFixed(2)));
            row.addEventListener('mouseenter', () => {
                tooltip.replaceChildren(detail(result, f));
                if (follower) tooltip.appendChild(el('div', 'gf-detail-meta',
                    `Moves with ${f.follows.label} · ${f.follows.part} across these proteins: one observation, not two.`));
                if (f.caveat) tooltip.appendChild(el('div', 'gf-detail-meta',
                    'Not biology: models and crystal structures can differ systematically (loops, termini, disorder).'));
                tooltip.hidden = false;
            });
            row.addEventListener('mousemove', event => {
                const r = tooltip.getBoundingClientRect();
                let x = event.clientX + 14, y = event.clientY + 12;
                if (x + r.width > window.innerWidth - 8) x = event.clientX - r.width - 14;
                if (y + r.height > window.innerHeight - 8) y = window.innerHeight - r.height - 8;
                tooltip.style.left = `${Math.max(8, x)}px`;
                tooltip.style.top = `${Math.max(8, y)}px`;
            });
            row.addEventListener('mouseleave', () => { tooltip.hidden = true; });
            list.appendChild(row);
        });
        root.appendChild(list);
        if (result.leaders.length > SHOWN) {
            const more = el('button', 'ov-view-btn gf-more', showAll ? 'Show fewer' : `Show all ${result.leaders.length}`);
            more.type = 'button';
            more.addEventListener('click', () => render(root, result, !showAll));
            root.appendChild(more);
        }
        root.appendChild(el('p', 'ov-hint gf-note', picked
            ? 'Numbers: |2 AUC − 1|, 1 = every picked protein lies above (or below) every reference protein, '
              + '0 = the two sets overlap evenly; beside each the medians (means when they tie), picked vs reference. Categories: '
              + 'the total variation distance between the two sets\' category shares, beside it the '
              + 'category that differs most (picked count vs reference share). Features below the chance '
              + 'line do no better than the best of all features for randomly picked proteins (same number). '
              + 'A feature is scored only when at least half of the picked proteins have a value. '
              + 'Hole-hopping positions: |2 AUC − 1| of the Jaccard distances among the picked proteins against those to the '
              + 'reference (1 = every two picked proteins share more of their path columns than any picked and reference protein; needs two picked proteins). '
              + 'Against the whole catalogue most answers are the family: the same-family reference avoids that.'
            : tree
            ? 'Spearman r over every pair of proteins between their distance in the tree (the height at '
              + 'which they join) and how different they are in the feature (the difference of a number; '
              + 'same or another category; for the hole-hopping positions the Jaccard distance, the share of '
              + 'their path columns the two do not share). +1 = the farther apart two proteins are in the tree, the more '
              + 'they differ; 0 = the feature ignores the tree; below 0 = distant proteins are more alike. '
              + 'No cut is needed. Features below the chance line do no better than the best of all '
              + 'features when the tree\'s proteins are shuffled. Pairs are not independent (every '
              + 'protein is in many), which is why the chance level comes from shuffling whole proteins. '
              + 'Features "from the sequence alone" are partly the distance itself.'
            : 'η² (numbers), Cramér\'s V (categories) and, for the hole-hopping positions, PERMANOVA R² '
            + 'of the Jaccard distances (η² for distances): 0 = the groups do not differ in this feature, '
            + '1 = it is constant inside each group and differs between them. Features below the chance '
            + 'line do no better than the best of all features on randomly shuffled groups. A protein '
            + 'on its own counts as a group: a feature it stands out in says what makes it different. The tree is '
            + 'built from sequence distance, so this shows what the groups share, not what causes them; '
            + 'features "from the sequence alone" are partly the distance itself.'));
    }

    function attach({ root, button, treeButton, pickButton, referenceSelect, sourceSelect, userGroups = [],
                      status, tree, data, compared = [] }) {
        if (!root || !button || !tree) return null;
        const tooltip = el('div', 'ov-tooltip gf-tooltip');
        tooltip.hidden = true;
        document.body.appendChild(tooltip);
        root.__tooltip = tooltip;
        let lastCut, lastMode, lastPick = '', currentCut = null;
        // "your groups": the named groups; one named group is compared with the ungrouped proteins
        const named = userGroups.filter(g => !g.ungrouped);
        const ownGroups = named.length >= 2 ? named : userGroups;
        const useOwn = () => !!sourceSelect && sourceSelect.value === 'yours' && ownGroups.length > 0;
        const canRun = () => useOwn() || currentCut != null;
        function grouping() {
            return useOwn() ? {groups: ownGroups.map(g => ({name: g.name, codes: g.codes})), cut: 'yours', user: true}
                            : tree.groups();
        }
        const pickKey = () => `${tree.picked().join(',')}|${referenceSelect ? referenceSelect.value : ''}`;
        // with nothing picked, the whole comparison is compared with the catalogue
        function syncReference() {
            if (!referenceSelect) return;
            const none = !tree.picked().length;
            const own = referenceSelect.querySelector('option[value="comparison"]');
            if (own) own.disabled = none;
            if (none && referenceSelect.value === 'comparison') referenceSelect.value = 'families';
            if (pickButton) pickButton.textContent = none
                ? 'Which features discriminate the compared proteins from' : 'Which features discriminate the picked proteins from';
        }
        function sync(stateCut) {
            currentCut = stateCut;
            button.disabled = !canRun();
            // a result goes stale when what it was computed for changes: the cut for groups, the
            // pick or the reference for picked proteins; the tree view has neither
            const ownNow = useOwn();
            const stale = (lastMode === 'groups' && (ownNow ? lastCut !== 'yours'
                              : (stateCut == null || lastCut === 'yours' || Math.abs(stateCut - lastCut) > 1e-9)))
                || (lastMode === 'picked' && pickKey() !== lastPick);
            if (stale) {
                status.textContent = 'The selection changed: press the button to update.';
                root.classList.add('gf-stale');
            } else {
                status.textContent = stateCut == null && lastMode !== 'picked' && !ownNow
                    ? 'Click proteins in the tree to pick them, or set "Group at" to compare groups.' : '';
                root.classList.remove('gf-stale');
            }
        }
        function run(compute, mode) {
            button.disabled = true;
            if (treeButton) treeButton.disabled = true;
            if (pickButton) pickButton.disabled = true;
            status.textContent = 'Scoring…';
            // let the status paint before the computation (the tree view shuffles every pair)
            setTimeout(() => {
                const result = compute();
                lastMode = mode;
                root.classList.remove('gf-stale');
                render(root, result, false);
                button.disabled = !canRun();
                if (treeButton) treeButton.disabled = false;
                if (pickButton) pickButton.disabled = false;
                status.textContent = '';
            }, 20);
        }
        button.addEventListener('click', () => {
            const groups = grouping();
            if (!groups) return;
            lastCut = groups.cut;
            run(() => analyse(data, groups), 'groups');
        });
        if (sourceSelect) sourceSelect.addEventListener('change', () => sync(currentCut));
        if (treeButton) {
            treeButton.addEventListener('click', () => run(() => analyseTree(data, tree.cophenetic()), 'tree'));
        }
        if (pickButton) {
            pickButton.addEventListener('click', () => {
                lastPick = pickKey();
                const reference = referenceSelect ? referenceSelect.value : 'comparison';
                run(() => analysePicked(data, tree.picked(), compared, reference), 'picked');
            });
            const changed = () => { syncReference(); sync(currentCut); };
            if (tree.onPickChange) tree.onPickChange(changed);
            if (referenceSelect) referenceSelect.addEventListener('change', changed);
            syncReference();
        }
        // with groups of its own, the comparison can be scored before any cut
        if (useOwn()) sync(currentCut);
        return { sync };
    }

    window.GroupingFeatures = { attach, analyse, analyseTree, analysePicked, featureList };
})();
