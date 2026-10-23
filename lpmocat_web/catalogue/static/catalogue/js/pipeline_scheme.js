// About page: the pipeline scheme (_pipeline_scheme.html). Draws an arrow from every box to the
// boxes in its data-to (side by side: left to right; in one column: top to bottom), and on hover
// or focus shows the box's details and highlights its arrows and the boxes they join. Below
// 1000 px the columns stack (info_pages.css) and the arrows between them are CSS.
(function () {
    const flow = document.getElementById('pipelineFlow');
    const tip = document.getElementById('pipelineFlowTip');
    if (!flow || !tip) return;
    const svg = flow.querySelector('.flow-links');
    const SVG = 'http://www.w3.org/2000/svg';
    const wide = window.matchMedia('(min-width: 1000px)');
    const nodes = [...flow.querySelectorAll('[data-node]')];
    const byId = new Map(nodes.map(n => [n.dataset.node, n]));
    let links = [];     // {froms: [box ids], to, path, arrow}

    function box(el) {
        const r = el.getBoundingClientRect();
        const f = flow.getBoundingClientRect();
        const left = r.left - f.left, top = r.top - f.top;
        return { left, top, right: left + r.width, bottom: top + r.height,
                 cx: left + r.width / 2, cy: top + r.height / 2, tall: r.height > 120 };
    }

    const clamp = (v, lo, hi) => Math.max(lo, Math.min(hi, v));

    function pathBetween(a, b) {
        if (b.left >= a.right - 1) {
            // a group meets the other box at that box's height
            const y1 = a.tall ? clamp(b.cy, a.top + 16, a.bottom - 16) : a.cy;
            const y2 = b.tall ? clamp(y1, b.top + 16, b.bottom - 16) : b.cy;
            // the curve ends in a short straight run, so the arrowhead lines up with the stroke
            // (a steep curve that turns flat only at its last pixels leaves the head askew)
            const x1 = a.right, x2 = b.left - 3, tail = Math.min(16, (x2 - x1) / 2);
            const dx = (x2 - tail - x1) / 2;
            return `M${x1},${y1} C${x1 + dx},${y1} ${x2 - tail - dx},${y2} ${x2 - tail},${y2} H${x2}`;
        }
        const down = b.top >= a.bottom - 1;
        const y1 = down ? a.bottom : a.top, y2 = down ? b.top - 3 : b.bottom + 3;
        const tail = (down ? 1 : -1) * Math.min(16, Math.abs(y2 - y1) / 2);
        const dy = (y2 - tail - y1) / 2;
        return `M${a.cx},${y1} C${a.cx},${y1 + dy} ${b.cx},${y2 - tail - dy} ${b.cx},${y2 - tail} V${y2}`;
    }

    function draw() {
        svg.replaceChildren();
        links = [];
        if (!wide.matches) return;
        svg.setAttribute('width', flow.scrollWidth);
        svg.setAttribute('height', flow.scrollHeight);
        const defs = document.createElementNS(SVG, 'defs');
        for (const [id, cls] of [['flowArrow', ''], ['flowArrowOn', 'on']]) {
            const marker = document.createElementNS(SVG, 'marker');
            // a fixed size: a highlighted (thicker) line would otherwise scale its arrowhead up
            Object.entries({ id, viewBox: '0 0 10 10', refX: 7, refY: 5, markerWidth: 9,
                             markerHeight: 9, markerUnits: 'userSpaceOnUse', orient: 'auto' })
                .forEach(([k, v]) => marker.setAttribute(k, v));
            const head = document.createElementNS(SVG, 'path');
            head.setAttribute('d', 'M0,1 L9,5 L0,9 z');
            head.setAttribute('class', `flow-arrowhead ${cls}`);
            marker.appendChild(head);
            defs.appendChild(marker);
        }
        svg.appendChild(defs);
        // Arrows by target; three or more boxes on the left of one target join on a bus: a stub
        // from each box to a vertical line, and one arrow from that line to the target
        const incoming = new Map();
        for (const from of nodes) {
            for (const to of (from.dataset.to || '').split(' ').filter(id => byId.has(id))) {
                if (!incoming.has(to)) incoming.set(to, []);
                incoming.get(to).push(from);
            }
        }
        for (const [to, sources] of incoming) {
            const b = box(byId.get(to));
            const boxes = sources.map(box);
            if (sources.length >= 3 && boxes.every(a => b.left >= a.right - 1)) {
                const busX = (Math.max(...boxes.map(a => a.right)) + b.left) / 2;
                const top = Math.min(...boxes.map(a => a.cy)), bottom = Math.max(...boxes.map(a => a.cy));
                const y = clamp(b.cy, top, bottom);
                const froms = sources.map(n => n.dataset.node);
                sources.forEach((source, i) =>
                    addLink(`M${boxes[i].right},${boxes[i].cy} H${busX}`, [source.dataset.node], to, false));
                addLink(`M${busX},${top} V${bottom}`, froms, to, false);
                const dx = (b.left - 3 - busX) / 2;
                addLink(`M${busX},${y} C${busX + dx},${y} ${b.left - 3 - dx},${b.cy} ${b.left - 3},${b.cy}`,
                        froms, to, true);
            } else {
                sources.forEach((source, i) =>
                    addLink(pathBetween(boxes[i], b), [source.dataset.node], to, true));
            }
        }
    }

    function addLink(d, froms, to, arrow) {
        const path = document.createElementNS(SVG, 'path');
        path.setAttribute('d', d);
        path.setAttribute('class', 'flow-link');
        if (arrow) path.setAttribute('marker-end', 'url(#flowArrow)');
        svg.appendChild(path);
        links.push({ froms, to, path, arrow });
    }

    // A box inside a group also counts as the group, so its way in and out lights up
    function highlight(node, on) {
        const ids = new Set([node.dataset.node, node.dataset.parent].filter(Boolean));
        node.classList.toggle('active', on);
        for (const link of links) {
            const out = link.froms.some(id => ids.has(id));
            if (!out && !ids.has(link.to)) continue;
            link.path.classList.toggle('on', on);
            if (link.arrow) link.path.setAttribute('marker-end', `url(#${on ? 'flowArrowOn' : 'flowArrow'})`);
            // the boxes at the other end (for a shared bus arrow: the target only)
            for (const id of out ? [link.to] : link.froms) {
                const other = byId.get(id);
                if (other && other !== node) other.classList.toggle('linked', on);
            }
            if (on) svg.appendChild(link.path);     // on top of the others
        }
    }

    // A box's own details (a group holds the details of its boxes too)
    const detailOf = node => node.querySelector(':scope > .flow-detail');

    // Beside the box, on its right if there is room, else on its left, else under or over it
    function showTip(node) {
        tip.innerHTML = detailOf(node).innerHTML;
        if (node.dataset.docs) {
            const link = document.createElement('a');
            link.className = 'flow-doc';
            link.href = node.dataset.docs;
            link.textContent = `${node.dataset.docsLabel} \u2192`;
            tip.appendChild(link);
        }
        tip.hidden = false;
        const r = node.getBoundingClientRect();
        const vw = document.documentElement.clientWidth, vh = window.innerHeight;
        const width = tip.offsetWidth, height = tip.offsetHeight;
        let x, y = r.top;
        if (r.right + 10 + width <= vw - 8) {
            x = r.right + 10;
        } else if (r.left - 10 - width >= 8) {
            x = r.left - 10 - width;
        } else {
            x = clamp(r.left, 8, vw - width - 8);
            y = r.bottom + 6 + height <= vh - 8 ? r.bottom + 6 : r.top - 6 - height;
        }
        tip.style.left = `${x}px`;
        tip.style.top = `${clamp(y, 8, vh - height - 8)}px`;
    }

    // The details can hold a link (the step's page of the pipeline documentation), so they
    // stay while the pointer moves from the box onto them, and while the link has the focus
    const LEAVE_DELAY = 200;    // ms to cross the gap between the box and its details
    let current = null, leaving = null;
    function enter(node) {
        clearTimeout(leaving);
        if (current === node) return;
        if (current) leave(current);
        current = node;
        highlight(node, true);
        showTip(node);
    }
    function leave(node) {
        clearTimeout(leaving);
        if (current !== node) return;
        highlight(node, false);
        tip.hidden = true;
        current = null;
    }
    function leaveSoon(node) {
        clearTimeout(leaving);
        leaving = setTimeout(() => {
            if (!tip.matches(':hover') && !tip.contains(document.activeElement)) leave(node);
        }, LEAVE_DELAY);
    }

    for (const node of nodes) {
        if (!detailOf(node)) continue;
        node.addEventListener('mouseenter', () => enter(node));
        node.addEventListener('mouseleave', () => leaveSoon(node));
        node.addEventListener('focus', () => enter(node));
        node.addEventListener('blur', () => leaveSoon(node));
        // a step box opens its documentation page on Enter (a resource box is a link itself)
        node.addEventListener('keydown', event => {
            if (event.key === 'Enter' && node.tagName !== 'A' && node.dataset.docs) {
                window.location.href = node.dataset.docs;
            }
        });
    }
    tip.addEventListener('mouseenter', () => clearTimeout(leaving));
    tip.addEventListener('mouseleave', () => { if (current) leaveSoon(current); });
    tip.addEventListener('focusout', () => { if (current) leaveSoon(current); });
    window.addEventListener('scroll', () => { if (current) leave(current); }, { passive: true });

    new ResizeObserver(draw).observe(flow);
    wide.addEventListener('change', draw);
    flow.querySelectorAll('img').forEach(img => img.addEventListener('load', draw));
    draw();
})();
