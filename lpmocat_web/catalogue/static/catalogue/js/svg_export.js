/* SVG downloads of the charts: the overview cards (browse page and selection overview) and the
   selection overview's sequence logo and similarity heatmap. The card's title and subtitle are
   HTML beside the chart, so they are written into the image on top of it.
   window.SvgExport.plotly(div, {title, sub, filename}) for a Plotly chart (redrawn off-screen),
   window.SvgExport.svg(node, {title, sub, legend, filename}) for a drawing of our own. */
(function () {
    'use strict';

    const SVG_NS = 'http://www.w3.org/2000/svg';
    const FONT = 'Arial, sans-serif';
    const TITLE_PX = 13;
    const SUB_PX = 11;
    const LINE = 15;            // px per subtitle line
    const CHAR_W = 0.55;        // average character width of Arial, in font sizes
    const MAX_TEXT_W = 900;     // px: a wide drawing (the logo) keeps readable lines

    const escapeHtml = text => text.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');

    function filenameOf(title, prefix) {
        const slug = String(title || 'chart').toLowerCase().replace(/[^a-z0-9]+/g, '_').replace(/^_+|_+$/g, '');
        return (prefix || 'lpmocat_') + slug;
    }

    function save(text, filename) {
        const blob = new Blob([text], { type: 'image/svg+xml' });
        const link = document.createElement('a');
        link.href = URL.createObjectURL(blob);
        link.download = filename.endsWith('.svg') ? filename : filename + '.svg';
        document.body.appendChild(link);
        link.click();
        link.remove();
        setTimeout(() => URL.revokeObjectURL(link.href), 1000);
    }

    // A subtitle as words: plain text, or a caption node whose colour swatches (.ov-swatch)
    // become coloured squares. -> [{text, color}]
    function wordsOf(sub) {
        const words = [];
        const walk = node => {
            if (node.nodeType === Node.TEXT_NODE) {
                node.textContent.split(/\s+/).filter(Boolean).forEach(text => words.push({ text }));
            } else if (node.nodeType === Node.ELEMENT_NODE && getComputedStyle(node).display !== 'none') {
                if (node.classList.contains('ov-swatch')) {
                    words.push({ text: '■', color: node.style.background || node.style.backgroundColor });
                } else node.childNodes.forEach(walk);
            }
        };
        if (!sub) return words;
        if (typeof sub === 'string') sub.split(/\s+/).filter(Boolean).forEach(text => words.push({ text }));
        else walk(sub);
        return words;
    }

    // words -> lines of at most `width` px
    function wrap(words, width) {
        const maxChars = Math.max(20, Math.floor(width / (SUB_PX * CHAR_W)));
        const lines = [];
        let line = [], length = 0;
        words.forEach(word => {
            if (line.length && length + 1 + word.text.length > maxChars) {
                lines.push(line);
                line = [];
                length = 0;
            }
            length += (line.length ? 1 : 0) + word.text.length;
            line.push(word);
        });
        if (line.length) lines.push(line);
        return lines;
    }

    // A Plotly chart: the same figure with a taller top margin holding the title and subtitle
    function plotly(div, { title, sub, filename }) {
        if (!div || !div._fullLayout || !window.Plotly) return Promise.resolve();
        const full = div._fullLayout;
        const width = Math.round(full.width);
        const lines = wrap(wordsOf(sub), width - 20);
        const head = 14 + LINE + lines.length * LINE;
        const layout = Object.assign({}, div.layout, {
            width,
            height: Math.round(full.height) + head,
            // the margins as drawn (_size: widened by automargin for long labels)
            margin: { l: full._size.l, r: full._size.r, b: full._size.b, t: full._size.t + head,
                      pad: full.margin.pad, autoexpand: false },
            paper_bgcolor: 'white',
            title: { text: '' },
        });
        // title and subtitle as annotations at the top of the image (paper y = 1 is the plot top)
        const plotH = layout.height - layout.margin.t - layout.margin.b;
        const yAt = px => 1 + (layout.margin.t - px) / plotH;
        const texts = [{ text: `<b>${escapeHtml(title || '')}</b>`, size: TITLE_PX, top: 10 }]
            .concat(lines.map((line, i) => ({
                text: line.map(w => w.color ? `<span style="color:${w.color}">${w.text}</span>` : escapeHtml(w.text)).join(' '),
                size: SUB_PX, top: 10 + LINE + 3 + i * LINE,
            })));
        layout.annotations = (div.layout.annotations || []).concat(texts.map(t => ({
            text: t.text, xref: 'paper', yref: 'paper', x: -(layout.margin.l - 10) / (width - layout.margin.l - layout.margin.r),
            y: yAt(t.top), xanchor: 'left', yanchor: 'top', showarrow: false, align: 'left',
            font: { family: FONT, size: t.size, color: t.size === TITLE_PX ? '#222' : '#555' },
        })));
        return window.Plotly.downloadImage({ data: div.data, layout }, {
            format: 'svg', width, height: layout.height, filename: filename || filenameOf(title),
        });
    }

    // One of our own drawings: a clean copy (no invisible hover areas, no cursor) under a title,
    // its subtitle and an optional legend [[label, colour], ...]
    function svg(node, { title, sub, legend, filename }) {
        if (!node) return;
        const copy = node.cloneNode(true);
        copy.querySelectorAll('[opacity="0"], [display="none"]').forEach(n => n.remove());
        const width = Number(copy.getAttribute('width')) || node.getBoundingClientRect().width;
        const height = Number(copy.getAttribute('height')) || node.getBoundingClientRect().height;
        const textW = Math.min(Math.max(width, 400), MAX_TEXT_W) - 20;
        const lines = wrap(wordsOf(sub), textW);
        const legendLines = legend && legend.length ? wrap(legend.map(([text, color]) => ({ text, color })), textW) : [];
        const head = 14 + LINE + (lines.length + legendLines.length) * LINE + 6;
        const outW = Math.max(width, textW + 20);

        const out = document.createElementNS(SVG_NS, 'svg');
        out.setAttribute('xmlns', SVG_NS);
        out.setAttribute('width', outW);
        out.setAttribute('height', height + head);
        out.setAttribute('viewBox', `0 0 ${outW} ${height + head}`);
        out.setAttribute('font-family', FONT);
        const background = document.createElementNS(SVG_NS, 'rect');
        Object.entries({ x: 0, y: 0, width: outW, height: height + head, fill: 'white' })
            .forEach(([k, v]) => background.setAttribute(k, v));
        out.appendChild(background);
        const addLine = (words, y, size, color, weight) => {
            const text = document.createElementNS(SVG_NS, 'text');
            text.setAttribute('x', 10);
            text.setAttribute('y', y);
            text.setAttribute('font-size', size);
            text.setAttribute('fill', color);
            if (weight) text.setAttribute('font-weight', weight);
            words.forEach((word, i) => {
                const span = document.createElementNS(SVG_NS, 'tspan');
                span.textContent = (i ? ' ' : '') + word.text;
                if (word.color) span.setAttribute('fill', word.color);
                text.appendChild(span);
            });
            out.appendChild(text);
        };
        addLine([{ text: title || '' }], 10 + TITLE_PX, TITLE_PX, '#222', 'bold');
        lines.concat(legendLines).forEach((line, i) => addLine(line, 10 + LINE + 3 + SUB_PX + i * LINE, SUB_PX, '#555'));
        copy.setAttribute('x', 0);
        copy.setAttribute('y', head);
        out.appendChild(copy);
        save(new XMLSerializer().serializeToString(out), filename || filenameOf(title));
    }

    window.SvgExport = { plotly, svg, filenameOf };
})();
