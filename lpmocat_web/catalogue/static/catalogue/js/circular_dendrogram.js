class CircularDendrogram {
    constructor(containerId, options = {}) {
        this.containerId = containerId;
        this.container = d3.select(`#${containerId}`);
        
        this.config = {
            width: options.width || 800,
            height: options.height || 800,
            margin: options.margin || { top: 20, right: 20, bottom: 20, left: 20 },
            innerRadius: options.innerRadius || 100,
            labelOffset: options.labelOffset || 10,
            duration: options.duration || 750,
            showBranchLengths: options.showBranchLengths || false,
            showLabels: options.showLabels !== false,
            initialScale: options.initialScale || 1.35
        };

        this.radius = Math.min(this.config.width, this.config.height) / 2 - Math.max(...Object.values(this.config.margin));
        this.outerRadius = this.radius - 60;

        // family -> colour, from dendrogram.CLADE_COLORS (static_dendrogram_circular)
        this.cladeColors = options.cladeColors || {};

        this.defaultColor = '#569976';
        
        this.initializeSVG();
        this.applyScale();
        this.setupTooltip();
        this.setupResize();
    }

    initializeSVG() {
        this.container.selectAll("*").remove();

        const containerRect = this.container.node().getBoundingClientRect();
        const width = Math.max(containerRect.width || this.config.width, 600);
        const height = Math.max(containerRect.height || this.config.height, 600);
        
        this.config.width = width;
        this.config.height = height;
        this.radius = Math.min(width, height) / 2 - 100;
        this.outerRadius = this.radius - 60;

        this.svg = this.container
            .append("svg")
            .attr("width", "100%")
            .attr("height", "100%")
            .attr("viewBox", `0 0 ${width} ${height}`)
            .attr("preserveAspectRatio", "xMidYMid meet")
            .style("font-family", "Arial, sans-serif");

        // zoomLayer carries the zoom; mainGroup keeps the tree's own centring and scale
        this.zoomLayer = this.svg.append("g").attr("class", "zoom-layer");
        this.mainGroup = this.zoomLayer
            .append("g")
            .attr("class", "main-group")
            .attr("transform", `translate(${width / 2}, ${height / 2})`);

        this.linkGroup = this.mainGroup.append("g").attr("class", "links");
        // rings of the browse page's groups ("Add as group") under the search and overview halos
        this.groupMarkGroup = this.mainGroup.append("g").attr("class", "group-marks");
        this.hitGroup = this.mainGroup.append("g").attr("class", "search-hits");
        this.nodeGroup = this.mainGroup.append("g").attr("class", "nodes");
        this.labelGroup = this.mainGroup.append("g").attr("class", "labels");
        this.setLabelsVisible(this.config.showLabels);
        this.setupZoom(width, height);
    }

    // Zoom without taking the page's scrolling: the wheel zooms only with Ctrl/Cmd (a trackpad
    // pinch sends that too), dragging pans only when zoomed in, and double-click does nothing
    // (clicks select). The tree cannot be pushed out of its window.
    setupZoom(width, height) {
        this.zoom = d3.zoom()
            .scaleExtent([1, 20])
            .extent([[0, 0], [width, height]])
            .translateExtent([[0, 0], [width, height]])
            // d3 speeds the wheel up ten times with Ctrl (for pinches); here Ctrl is always held,
            // so a mouse-wheel notch (deltaY 100) zooms by 1.3 instead of 4
            .wheelDelta(event => -event.deltaY * (event.deltaMode === 1 ? 0.1 : event.deltaMode ? 2 : 0.004))
            .filter(event => {
                if (event.type === 'wheel') return event.ctrlKey || event.metaKey;
                if (event.type === 'dblclick') return false;
                if (event.type === 'mousedown') return !event.button && !event.shiftKey && this.zoomScale() > 1;
                return true;
            })
            .on("zoom", event => {
                this.zoomK = event.transform.k;
                this.zoomLayer.attr("transform", event.transform);
                this.svg.classed("zoomed", event.transform.k > 1);
                this.applyZoomToMarks();
            });
        this.svg.call(this.zoom).on("dblclick.zoom", null);
    }

    zoomScale() {
        return d3.zoomTransform(this.svg.node()).k;
    }

    // Zooming spreads the tree out, but dots, halos, names and the gap before a name keep their
    // size on the screen (the name's font size follows --zoom in main_page.css).
    markSize(size) {
        return size / (this.zoomK || 1);
    }

    nodeRadius(d) {
        const selected = d.data && d.data.isLeaf && this.selectedProteins && this.selectedProteins.has(d.data.name);
        const hovered = d.data && d.data.isLeaf && this.hoveredLeaf === d.data.name;
        return this.markSize((selected ? 5 : 3) + (hovered ? 2 : 0));
    }

    // A leaf's dot and name are marked together when the cursor is on either (null = none).
    setLeafHover(name) {
        if (this.hoveredLeaf === name) return;
        this.hoveredLeaf = name;
        const mine = d => d.data && d.data.isLeaf && d.data.name === name;
        this.labelGroup.selectAll(".label").classed("leaf-hover", mine);
        this.nodeGroup.selectAll(".node")
            .classed("leaf-hover", mine)
            .filter(d => d.data && d.data.isLeaf)
            .attr("r", d => this.nodeRadius(d));
    }

    labelTransform(d) {
        const radius = this.config.showBranchLengths ? (d.radius || d.y) : d.y;
        const angle = d.x - 90;
        const labelRadius = radius + this.markSize(this.config.labelOffset);
        return `rotate(${angle}) translate(${labelRadius}, 0) ${angle > 90 && angle < 270 ? "rotate(180)" : ""}`;
    }

    applyZoomToMarks() {
        this.svg.style("--zoom", this.zoomK || 1);
        this.nodeGroup.selectAll(".node").attr("r", d => this.nodeRadius(d));
        this.hitGroup.selectAll(".search-hit-halo")
            .attr("r", (d, i, nodes) => this.markSize(nodes[i].classList.contains("best") ? 11 : 8));
        this.hitGroup.selectAll(".overview-halo").attr("r", this.markSize(8));
        if (this.groupMarkGroup) this.groupMarkGroup.selectAll(".group-ring").attr("r", this.markSize(7.5));
        this.labelGroup.selectAll(".label").attr("transform", d => this.labelTransform(d));
    }

    zoomBy(factor) {
        this.svg.transition().duration(300).call(this.zoom.scaleBy, factor);
    }

    resetZoom() {
        this.svg.transition().duration(400).call(this.zoom.transform, d3.zoomIdentity);
    }

    // Fit the branch under `node` (its lines, dots and, when shown, names) into the window.
    zoomToBranch(node) {
        const radiusOf = d => this.config.showBranchLengths ? (d.radius || d.y) : d.y;
        const point = (r, x) => {
            const a = (x - 90) / 180 * Math.PI;
            return [r * Math.cos(a), r * Math.sin(a)];
        };
        const points = [];
        const labelEnds = [];
        node.descendants().forEach(d => {
            points.push(point(radiusOf(d), d.x));
            (d.children || []).forEach(child => points.push(point(radiusOf(d), child.x)));
            if (!d.children && this.config.showLabels) labelEnds.push(d);
        });
        if (node.parent) points.push(point(radiusOf(node.parent), node.x));

        // mainGroup coordinates -> the svg's (zoom-free) coordinates. A name keeps its screen size
        // (about 70 px with its gap), so its end depends on the zoom: fitted twice.
        const { width, height, initialScale } = this.config;
        const fit = k => {
            const all = points.concat(labelEnds.map(d => point(radiusOf(d) + 70 / (k * initialScale), d.x)));
            const xs = all.map(p => width / 2 + initialScale * p[0]);
            const ys = all.map(p => height / 2 + initialScale * p[1]);
            const box = [Math.min(...xs), Math.max(...xs), Math.min(...ys), Math.max(...ys)];
            const fitted = 0.9 * Math.min(width / (box[1] - box[0] + 20), height / (box[3] - box[2] + 20));
            return { box, k: Math.max(1, Math.min(20, fitted)) };
        };
        let { k } = fit(1);
        ({ k } = fit(k));
        const [x0, x1, y0, y1] = fit(k).box;
        const transform = d3.zoomIdentity
            .translate(width / 2, height / 2)
            .scale(k)
            .translate(-(x0 + x1) / 2, -(y0 + y1) / 2);
        this.svg.transition().duration(600).call(this.zoom.transform, transform);
    }

    // Protein names on or off; the dots and branches stay clickable, and Export SVG follows it.
    setLabelsVisible(visible) {
        this.config.showLabels = visible;
        this.labelGroup.attr("display", visible ? null : "none");
    }

    // Fixed scale, no d3.zoom: wheel zooming and dragging moved the tree while the page scrolled.
    applyScale() {
        this.mainGroup.attr("transform",
            `translate(${this.config.width / 2}, ${this.config.height / 2}) scale(${this.config.initialScale})`
        );
    }

    setupTooltip() {
        this.tooltip = d3.select("body")
            .append("div")
            // .attr("class", "dendrogram-tooltip")
            .style("position", "absolute")
            .style("visibility", "hidden")
            .style("background", "rgba(0, 0, 0, 0.9)")
            .style("color", "white")
            .style("padding", "12px")
            .style("border-radius", "8px")
            .style("font-size", "12px")
            .style("max-width", "320px")
            .style("z-index", "1000")
            .style("pointer-events", "auto")
            .style("box-shadow", "0 4px 12px rgba(0,0,0,0.3)")
            .style("border", "1px solid rgba(255,255,255,0.1)");
            
        this.tooltip.on("mouseenter", () => this.cancelHideTooltip());
        this.tooltip.on("mouseleave", () => this.scheduleHideTooltip());
    }

    setupResize() {
        const resizeHandler = () => {
            if (this.lastMatrix && this.lastLabels) {
                this.initializeSVG();
                this.applyScale();
                this.render(this.lastMatrix, this.lastLabels, this.proteinInfo, Array.from(this.selectedProteins || []));
            }
        };

        let resizeTimeout;
        window.addEventListener('resize', () => {
            clearTimeout(resizeTimeout);
            resizeTimeout = setTimeout(resizeHandler, 250);
        });
    }

    matrixToTree(distanceMatrix, labels) {
        if (typeof HierarchicalClustering !== 'undefined') {
            const clustering = new HierarchicalClustering('average');
            const tree = clustering.cluster(distanceMatrix, labels);
            return clustering.optimizeTree(tree);
        }
        
        return this.simpleMatrixToTree(distanceMatrix, labels);
    }

    simpleMatrixToTree(distanceMatrix, labels) {
        let nodes = labels.map((label, i) => ({
            name: label,
            id: i,
            isLeaf: true,
            distance: 0,
            children: null
        }));

        let clusters = [...nodes];
        let nodeId = labels.length;

        while (clusters.length > 1) {
            let minDist = Infinity;
            let minI = 0, minJ = 1;

            for (let i = 0; i < clusters.length; i++) {
                for (let j = i + 1; j < clusters.length; j++) {
                    const dist = this.getDistance(clusters[i], clusters[j], distanceMatrix, labels);
                    if (dist < minDist) {
                        minDist = dist;
                        minI = i;
                        minJ = j;
                    }
                }
            }

            const newNode = {
                name: `internal_${nodeId}`,
                id: nodeId++,
                isLeaf: false,
                distance: minDist / 2,
                children: [clusters[minI], clusters[minJ]]
            };

            const newClusters = clusters.filter((_, idx) => idx !== minI && idx !== minJ);
            newClusters.push(newNode);
            clusters = newClusters;
        }

        return clusters[0];
    }

    getDistance(cluster1, cluster2, matrix, labels) {
        const getLeafIndices = (cluster) => {
            if (cluster.isLeaf) {
                return [labels.indexOf(cluster.name)];
            }
            return cluster.children.flatMap(child => getLeafIndices(child));
        };

        const indices1 = getLeafIndices(cluster1);
        const indices2 = getLeafIndices(cluster2);

        let totalDist = 0;
        let count = 0;

        for (const i of indices1) {
            for (const j of indices2) {
                totalDist += matrix[i][j];
                count++;
            }
        }

        return count > 0 ? totalDist / count : 0;
    }

    render(distanceMatrix, labels, proteinInfo = {}, selectedProteins = []) {
        this.proteinInfo = proteinInfo;
        this.selectedProteins = new Set(selectedProteins);

        const treeData = this.matrixToTree(distanceMatrix, labels);
        
        const root = d3.hierarchy(treeData, d => d.children);
        
        const cluster = d3.cluster()
            .size([360, this.outerRadius])
            .separation((a, b) => 1);

        cluster(root);

        if (this.config.showBranchLengths) {
            this.scaleBranchLengths(root);
        }
        // Smaller names with branch lengths: leaves nearer the centre have less room
        this.svg.classed("branch-lengths", this.config.showBranchLengths);

        this.lastRoot = root;

        this.renderLinks(root);
        this.renderNodes(root);
        this.renderLabels(root);
        this.renderSearchHits();
        this.renderGroupMarks();

        setTimeout(() => {
            this.ensureVisibility();
        }, this.config.duration + 100);

        return this;
    }

    scaleBranchLengths(root) {
        const maxDistance = this.getMaxDistance(root);
        const scale = this.outerRadius / maxDistance;

        const setRadius = (node, parentRadius = 0) => {
            const distance = node.data.distance || 0;
            node.radius = parentRadius + (distance * scale);
            if (node.children) {
                node.children.forEach(child => setRadius(child, node.radius));
            }
        };

        setRadius(root);
    }

    getMaxDistance(node) {
        if (!node.children) {
            return node.data.distance || 0;
        }
        return (node.data.distance || 0) + Math.max(...node.children.map(child => this.getMaxDistance(child)));
    }

    renderLinks(root) {
        const links = root.links();

        const linkPath = this.linkGroup
            .selectAll(".link")
            .data(links);

        linkPath.exit().remove();

        const linkEnter = linkPath
            .enter()
            .append("path")
            .attr("class", "link")
            .attr("fill", "none")
            .attr("stroke", "#999")
            .attr("stroke-width", 1.5);

        linkPath.merge(linkEnter)
            .transition()
            .duration(this.config.duration)
            .attr("d", d => this.linkArc(d))
            .attr("stroke", d => this.getLinkColor(d.target));

        // A branch line acts for everything below it: one protein, or a whole branch.
        linkPath.merge(linkEnter)
            .on("click", (event, d) => {
                event.stopPropagation();
                if (event.shiftKey && !d.target.data.isLeaf) {
                    this.zoomToBranch(d.target);
                } else if (d.target.data.isLeaf) {
                    this.toggleProteinSelection(d.target.data.name);
                } else {
                    this.toggleBranchSelection(d.target);
                    this.showBranchTooltip(event, d.target, true);
                }
            })
            .on("mouseover", (event, d) => {
                this.highlightBranch(d.target);
                if (!d.target.data.isLeaf) this.showBranchTooltip(event, d.target);
            })
            .on("mousemove", (event, d) => {
                if (!d.target.data.isLeaf && !this.tooltipFor) this.moveTooltip(event);
            })
            .on("mouseout", () => {
                this.highlightBranch(null);
                this.hideBranchTooltip();
            });
    }

    linkArc(d) {
        const sourceRadius = this.config.showBranchLengths ? 
            (d.source.radius || d.source.y) : d.source.y;
        const targetRadius = this.config.showBranchLengths ? 
            (d.target.radius || d.target.y) : d.target.y;

        const sourceAngle = (d.source.x - 90) / 180 * Math.PI;
        const targetAngle = (d.target.x - 90) / 180 * Math.PI;

        const sx = sourceRadius * Math.cos(sourceAngle);
        const sy = sourceRadius * Math.sin(sourceAngle);
        const tx = targetRadius * Math.cos(targetAngle);
        const ty = targetRadius * Math.sin(targetAngle);

        const mx = sourceRadius * Math.cos(targetAngle);
        const my = sourceRadius * Math.sin(targetAngle);

        return `M${sx},${sy}A${sourceRadius},${sourceRadius} 0 0,${targetAngle > sourceAngle ? 1 : 0} ${mx},${my}L${tx},${ty}`;
    }

    getLinkColor(node) {
        if (node.data.isLeaf && this.proteinInfo[node.data.name]) {
            const clade = this.proteinInfo[node.data.name].clade;
            return this.cladeColors[clade] || this.defaultColor;
        }
        return "#999";
    }

    renderNodes(root) {
        const nodes = root.descendants();

        const nodeCircle = this.nodeGroup
            .selectAll(".node")
            .data(nodes);

        nodeCircle.exit().remove();

        const nodeEnter = nodeCircle
            .enter()
            .append("circle")
            .attr("class", "node")
            .attr("r", d => this.nodeRadius(d))
            .attr("fill", d => d.data.isLeaf ? this.getNodeColor(d) : "#999")
            .attr("stroke", d => {
                if (d.data.isLeaf && this.selectedProteins && this.selectedProteins.has(d.data.name)) {
                    return "#FFD700";
                }
                return "#fff";
            })
            .attr("stroke-width", d => {
                if (d.data.isLeaf && this.selectedProteins && this.selectedProteins.has(d.data.name)) {
                    return 2;
                }
                return 1;
            })
            .style("cursor", d => d.data.isLeaf ? "pointer" : "pointer");

        nodeCircle.merge(nodeEnter)
            .transition()
            .duration(this.config.duration)
            .attr("transform", d => {
                const radius = this.config.showBranchLengths ? 
                    (d.radius || d.y) : d.y;
                const angle = (d.x - 90) / 180 * Math.PI;
                const x = radius * Math.cos(angle);
                const y = radius * Math.sin(angle);
                return `translate(${x}, ${y})`;
            })
            .attr("r", d => this.nodeRadius(d))
            .attr("fill", d => d.data.isLeaf ? this.getNodeColor(d) : "#999")
            .attr("stroke", d => {
                if (d.data.isLeaf && this.selectedProteins && this.selectedProteins.has(d.data.name)) {
                    return "#FFD700";
                }
                return "#fff";
            })
            .attr("stroke-width", d => {
                if (d.data.isLeaf && this.selectedProteins && this.selectedProteins.has(d.data.name)) {
                    return 2;
                }
                return 1;
            });

        nodeCircle.merge(nodeEnter)
            .on("click", (event, d) => {
                event.stopPropagation();
                if (event.shiftKey && !d.data.isLeaf) {
                    this.zoomToBranch(d);
                } else if (!d.data.isLeaf) {
                    this.toggleBranchSelection(d);
                    this.showBranchTooltip(event, d, true);
                } else {
                    this.toggleProteinSelection(d.data.name);
                }
            })
            // A leaf dot opens the protein's tooltip as its name does (the names can be hidden).
            .on("mouseover", (event, d) => {
                if (d.data.isLeaf) {
                    this.setLeafHover(d.data.name);
                    this.hoverLabel(event, d);
                    return;
                }
                this.highlightBranch(d);
                this.showBranchTooltip(event, d);
            })
            .on("mousemove", (event, d) => {
                if (!d.data.isLeaf && !this.tooltipFor) this.moveTooltip(event);
            })
            .on("mouseout", (event, d) => {
                if (d.data.isLeaf) {
                    this.setLeafHover(null);
                    this.cancelShowTooltip();
                    this.scheduleHideTooltip();
                    return;
                }
                this.highlightBranch(null);
                this.hideBranchTooltip();
            });
    }

    getNodeColor(node) {
        if (node.data.isLeaf && this.proteinInfo[node.data.name]) {
            const clade = this.proteinInfo[node.data.name].clade;
            const baseColor = this.cladeColors[clade] || this.defaultColor;
            
            if (this.selectedProteins && this.selectedProteins.has(node.data.name)) {
                return baseColor;
            } else {
                return baseColor + '80';
            }
        }
        return this.defaultColor;
    }

    // Labels are ink, not the family colour (several family colours are too light for text on
    // white); the node dot and branch line beside a label carry its family.
    getLabelColor(node) {
        const selected = node && node.data && this.selectedProteins && this.selectedProteins.has(node.data.name);
        return selected ? "#1a1a1a" : "#595959";
    }

    renderLabels(root) {
        const leaves = root.leaves();

        const label = this.labelGroup
            .selectAll(".label")
            .data(leaves);

        label.exit().remove();

        const labelEnter = label
            .enter()
            .append("text")
            .attr("class", "label")
            .attr("dy", "0.31em")
            .attr("font-size", "11px")
            .attr("font-family", "Arial, sans-serif")
            .style("cursor", "pointer")
            // .text(d => d.data.name);
            .text(d => {
            if (this.proteinInfo && this.proteinInfo[d.data.name]) {
                return this.proteinInfo[d.data.name].code;
            }
            return d.data.name;
            });

        const labelUpdate = label.merge(labelEnter);

        labelUpdate
            .attr("font-size", this.config.showBranchLengths ? "8px" : "11px")  // for Export SVG; the screen size is in main_page.css
            .transition()
            .duration(this.config.duration)
            .attr("transform", d => this.labelTransform(d))
            .attr("text-anchor", d => {
                const angle = d.x - 90;
                return angle > 90 && angle < 270 ? "end" : "start";
            })
            .attr("fill", d => this.getLabelColor(d))
            .attr("font-weight", d => {
                if (this.selectedProteins && this.selectedProteins.has(d.data.name)) {
                    return "bold";
                }
                return "normal";
            })
            .style("cursor", "pointer");

        labelUpdate
            .on("mouseover", (event, d) => {
                this.setLeafHover(d.data.name);
                this.hoverLabel(event, d);
            })
            // Delayed, so the cursor can move onto the tooltip (its links) without losing it.
            .on("mouseout", () => {
                this.setLeafHover(null);
                this.cancelShowTooltip();
                this.scheduleHideTooltip();
            })
            .on("click", (event, d) => {
                event.stopPropagation();
                this.toggleProteinSelection(d.data.name);
            });
    }

    showTooltip(event, d) {
        if (!this.proteinInfo[d.data.name]) return;

        const info = this.proteinInfo[d.data.name];
        const proteinUrl = `/catalogue/compare/proteins/${info.code}/`;
        
        // Handle nan values for substrates, regioselectivity, and structure
        const substrates = (info.substrates && info.substrates.toLowerCase() !== 'nan') ? info.substrates : '';
        const regioselectivity = (info.regioselectivity && info.regioselectivity.toLowerCase() !== 'nan') ? info.regioselectivity : '';
        const associated_CBMs = (info.associated_CBMs && info.associated_CBMs.trim() !== '') ? info.associated_CBMs : '';
        
        // Process structure field to display "AlphaFold3 model" or "PDB code" with link
        let structure = '';
        let uniprotLinks = '';
        if (info.selected_structure && info.selected_structure.toLowerCase() !== 'nan') {
            if (info.selected_structure.toUpperCase().startsWith('AF-')) {
                structure = 'AlphaFold model';
                // Add UniProt links after AlphaFold3 model if available
                if (info.uniprot && info.uniprot.trim() && info.uniprot.toLowerCase() !== 'nan' && info.uniprot.toLowerCase() !== 'none' && info.uniprot.toLowerCase() !== 'null') {
                    const uniprotList = info.uniprot.split(/[,;]/).map(entry => entry.trim()).filter(entry => {
                        return entry && entry.toLowerCase() !== 'nan' && entry.toLowerCase() !== 'none' && entry.toLowerCase() !== 'null' && entry.length > 0;
                    });
                    if (uniprotList.length > 0) {
                        const uniprotLinkHtml = uniprotList.map(entry => {
                            return `<a href="https://www.uniprot.org/uniprotkb/${entry}/entry" target="_blank" style="color: #569976; text-decoration: underline;">${entry}</a>`;
                        }).join(', ');
                        uniprotLinks = ` ${uniprotLinkHtml}`;
                    }
                }
                // No active UniProt entry: link the UniParc record of the sequence
                if (!uniprotLinks && info.uniparc) {
                    uniprotLinks = info.uniparc.split(/[,;]/).map(u => u.trim()).filter(Boolean).map(u =>
                        ` <a href="https://www.uniprot.org/uniparc/${u}/entry" target="_blank" style="color: #569976; text-decoration: underline;">UniParc ${u}</a>`).join(',');
                }
            } else if (info.selected_structure.toUpperCase().startsWith('C3D-')) {
                // CAZyme3D model: C3D-<GenBank accession, version dot written as '_'>
                const genbank = info.selected_structure.slice(4).replace(/_(\d+)$/, '.$1');
                structure = `CAZyme3D model <a href="https://pro.unl.edu/CAZyme3D/get_entry.php?cazyid=${genbank}" target="_blank" style="color: #569976; text-decoration: underline;">${genbank}</a>`;
            } else {
                const pdbUrl = `https://www.rcsb.org/structure/${info.selected_structure}`;
                structure = `PDB code <a href="${pdbUrl}" target="_blank" style="color: #569976; text-decoration: underline;">${info.selected_structure}</a>`;
            }
        }
        
        const sequenceLength = info.sequence_length || 0;
        
        // Format charge display
        let chargeDisplay = 'N/A';
        let chargeNote = '';
        if (info.total_charge !== null && info.total_charge !== undefined) {
            chargeDisplay = `${info.total_charge.toFixed(0)}`
        }
        
        const otherNames = (info.other_names && info.other_names.trim()) ? info.other_names : '';

        const html = `
            <strong>${info.code}</strong>${otherNames ? `<br/><span style="color:#888; font-size:0.9em;">${otherNames}</span>` : ''}<br/>
            <strong>Clade:</strong> ${info.clade}<br/>
            <strong>Taxonomy:</strong> ${info.taxonomy || info.kingdom || ''}<br/>
            <strong>Sequence Length:</strong> ${sequenceLength}${info.catalytic_length ? ` (catalytic domain ${info.catalytic_length})` : ''} residues<br/>
            ${info.mass_kda != null ? `<strong>Mass:</strong> ${Number(info.mass_kda).toFixed(1)} kDa<br/>` : ''}
            <strong>Net charge:</strong> ${chargeDisplay}<br/>
            <strong>Structure:</strong> ${structure}${uniprotLinks}<br/>
            <strong>Organism:</strong> ${info.organism}<br/>
            <strong>EC:</strong> ${info.ec_numbers || ''}<br/>
            <strong>Substrates:</strong> ${substrates}<br/>
            <strong>Regioselectivity:</strong> ${regioselectivity}<br/>
            <strong>Associated CBMs:</strong> ${associated_CBMs}<br/>
            <a href="${proteinUrl}" target="_blank" style="color: #569976; text-decoration: underline; font-weight: bold;">
                Show this protein →
            </a>
        `;

        this.cancelHideTooltip();
        this.tooltipFor = d.data.name;
        this.tooltip
            .style("visibility", "visible")
            .html(html);

        this.placeProteinTooltip(event, d);
    }

    // A protein tooltip holds links, so it is easy to reach: it opens beside the label on the
    // outer side of the circle (no branches to cross) and stays there, it stays a moment after
    // the cursor leaves the label, and it switches to another protein only when the cursor rests
    // on that label, so labels crossed on the way to it don't replace it.
    hoverLabel(event, d) {
        this.cancelShowTooltip();
        const open = this.tooltip.style("visibility") === "visible" && this.tooltipFor;
        if (!open || this.tooltipFor === d.data.name) {
            this.showTooltip(event, d);
            return;
        }
        const at = { pageX: event.pageX, pageY: event.pageY };
        this.showTooltipTimer = setTimeout(() => this.showTooltip(at, d), 250);
    }

    cancelShowTooltip() {
        clearTimeout(this.showTooltipTimer);
        this.showTooltipTimer = null;
    }

    placeProteinTooltip(event, d) {
        const node = this.tooltip.node();
        // Measured at the page's left edge: near the right edge it shrinks to the space left
        this.tooltip.style("left", window.scrollX + "px");
        const width = node.offsetWidth;
        const angle = d.x - 90;
        const leftSpace = event.pageX - 16 - width >= window.scrollX + 4;
        const rightSpace = event.pageX + 16 + width <= window.scrollX + document.documentElement.clientWidth - 4;
        const leftHalf = angle > 90 && angle < 270;
        const toLeft = leftHalf ? leftSpace || !rightSpace : !rightSpace && leftSpace;
        let x = toLeft ? event.pageX - 16 - width : event.pageX + 16;
        let y = event.pageY - 16;
        const bottom = window.scrollY + window.innerHeight - 8;
        if (y + node.offsetHeight > bottom) y = Math.max(window.scrollY + 8, bottom - node.offsetHeight);
        x = Math.max(window.scrollX + 4, x);
        this.tooltip.style("left", x + "px").style("top", y + "px");
    }

    moveTooltip(event) {
        this.tooltip
            .style("top", (event.pageY - 10) + "px")
            .style("left", (event.pageX + 10) + "px");
    }

    hideTooltip() {
        this.cancelHideTooltip();
        this.cancelShowTooltip();
        this.tooltipFor = null;
        this.tooltip.style("visibility", "hidden");
    }

    // Hide a branch tooltip; an open protein tooltip stays (it hides on its own timer).
    hideBranchTooltip() {
        if (!this.tooltipFor) this.hideTooltip();
    }

    // Hide shortly unless the cursor is on the tooltip; showing a tooltip in the meantime
    // cancels it (otherwise leaving one protein hid the tooltip of the next one).
    scheduleHideTooltip() {
        this.cancelHideTooltip();
        this.hideTooltipTimer = setTimeout(() => {
            this.hideTooltipTimer = null;
            if (!this.tooltip.node().matches(':hover')) this.hideTooltip();
        }, 350);
    }

    cancelHideTooltip() {
        clearTimeout(this.hideTooltipTimer);
        this.hideTooltipTimer = null;
    }

    // The protein checkboxes for every leaf under an internal node (a branch / subtree).
    branchCheckboxes(node) {
        return this.proteinCheckboxes(node.leaves().map(leaf => leaf.data.name));
    }

    // the boxes of `names`, without the proteins of a group (browse page, "Add as group"), which
    // cannot be selected again
    proteinCheckboxes(names) {
        const grouped = name => !!(window.isGroupedProtein && window.isGroupedProtein(name));
        return names
            .filter(name => !grouped(name))
            .map(name => document.querySelector(`input[name="filter_protein_ids"][value="${CSS.escape(name)}"]`))
            .filter(Boolean);
    }

    // Select every protein in the branch; if they are all selected already, deselect them.
    toggleBranchSelection(node) {
        this.toggleNamesSelection(node.leaves().map(leaf => leaf.data.name));
    }

    // Select the proteins `names`; if they are all selected already, deselect them.
    toggleNamesSelection(names) {
        const checkboxes = this.proteinCheckboxes(names);
        const select = !checkboxes.every(cb => cb.checked);
        const changed = checkboxes.filter(cb => cb.checked !== select);
        changed.forEach(cb => { cb.checked = select; });
        // One change event is enough: the page handler re-reads every checkbox.
        if (changed.length) {
            changed[0].dispatchEvent(new Event('change', { bubbles: true }));
        }
    }

    // force: shown even over an open protein tooltip (a click on the branch)
    showBranchTooltip(event, node, force = false) {
        if (this.tooltipFor && !force) return;
        const leaves = node.leaves();
        const families = {};
        leaves.forEach(leaf => {
            const info = this.proteinInfo[leaf.data.name];
            const family = info ? info.clade : 'Unknown';
            families[family] = (families[family] || 0) + 1;
        });
        const familyText = Object.entries(families)
            .sort((a, b) => b[1] - a[1])
            .map(([family, count]) => `${family} &times;${count}`)
            .join(', ');
        const allSelected = this.branchCheckboxes(node).every(cb => cb.checked);

        this.cancelHideTooltip();
        this.cancelShowTooltip();
        this.tooltipFor = null;
        this.tooltip
            .style("visibility", "visible")
            .html(`<strong>Branch: ${leaves.length} proteins</strong><br/>${familyText}<br/>
                   <span style="color: #9fd4b4;">Click to ${allSelected ? 'deselect' : 'select'} all &middot; Shift+click to zoom in</span>`);
        this.moveTooltip(event);
    }

    // Emphasise the lines and labels under `node` (null clears it).
    highlightBranch(node) {
        const inBranch = new Set(node ? node.descendants() : []);
        this.linkGroup.selectAll(".link").classed("branch-hover", d => inBranch.has(d.target));
        this.labelGroup.selectAll(".label").classed("branch-hover", d => inBranch.has(d));
    }

    // The clade's proteins in the drawn tree
    cladeLeaves(clade) {
        return (this.lastRoot ? this.lastRoot.leaves() : [])
            .map(leaf => leaf.data.name)
            .filter(name => this.proteinInfo[name] && this.proteinInfo[name].clade === clade);
    }

    // Clade legend (dendrogram.py): hovering a clade highlights its proteins, a click selects
    // them (or deselects them when all are selected), like a branch.
    bindLegend(legend) {
        this.legend = legend;
        if (!legend) return;
        legend.querySelectorAll('[data-family]').forEach(item => {
            const clade = item.dataset.family;
            const on = () => this.setOverviewHighlight(this.cladeLeaves(clade));
            const off = () => this.setOverviewHighlight(null);
            item.addEventListener('mouseenter', on);
            item.addEventListener('mouseleave', off);
            item.addEventListener('focus', on);
            item.addEventListener('blur', off);
            item.addEventListener('click', event => {
                event.stopPropagation();
                this.toggleNamesSelection(this.cladeLeaves(clade));
            });
        });
        this.syncLegend();
    }

    // A clade whose proteins are all selected is marked in the legend.
    syncLegend() {
        if (!this.legend) return;
        this.legend.querySelectorAll('[data-family]').forEach(item => {
            const names = this.cladeLeaves(item.dataset.family);
            const all = names.length > 0 && names.every(name => this.selectedProteins && this.selectedProteins.has(name));
            item.classList.toggle('all-selected', all);
            item.title = `${item.dataset.family}: ${names.length} protein${names.length === 1 ? '' : 's'} `
                + `- hover to highlight, click to ${all ? 'deselect' : 'select'} all`;
        });
    }

    toggleBranchLengths() {
        this.config.showBranchLengths = !this.config.showBranchLengths;
        this.render(this.lastMatrix, this.lastLabels, this.proteinInfo);
        return this.config.showBranchLengths;
    }

    setData(matrix, labels) {
        this.lastMatrix = matrix;
        this.lastLabels = labels;
    }

    updateSelection(selectedProteins) {
        this.selectedProteins = new Set(selectedProteins);
        this.syncLegend();
        
        this.nodeGroup.selectAll(".node")
            .transition()
            .duration(200)
            .attr("r", d => this.nodeRadius(d))
            .attr("fill", d => d.data && d.data.isLeaf ? this.getNodeColor(d) : "#999")
            .attr("stroke", d => {
                if (d.data && d.data.isLeaf && this.selectedProteins.has(d.data.name)) {
                    return "#FFD700";
                }
                return "#fff";
            })
            .attr("stroke-width", d => {
                if (d.data && d.data.isLeaf && this.selectedProteins.has(d.data.name)) {
                    return 2;
                }
                return 1;
            });

        this.labelGroup.selectAll(".label")
            .transition()
            .duration(200)
            .attr("transform", d => d && d.data ? this.labelTransform(d) : "")
            .attr("text-anchor", d => {
                if (!d || !d.data) return "start";
                const angle = d.x - 90;
                return angle > 90 && angle < 270 ? "end" : "start";
            })
            .attr("fill", d => this.getLabelColor(d))
            .attr("font-weight", d => {
                if (d && d.data && this.selectedProteins.has(d.data.name)) {
                    return "bold";
                }
                return "normal";
            })
            .style("cursor", "pointer")
            .style("visibility", "visible");
    }

    forceRerender() {
        if (this.lastMatrix && this.lastLabels && this.proteinInfo) {
            this.nodeGroup.selectAll("*").remove();
            this.labelGroup.selectAll("*").remove();
            this.linkGroup.selectAll("*").remove();
            
            this.render(this.lastMatrix, this.lastLabels, this.proteinInfo, Array.from(this.selectedProteins || []));
        }
    }

    ensureVisibility() {
        if (this.container) {
            this.container.style("visibility", "visible");
            this.container.style("opacity", "1");
        }
        
        this.labelGroup.selectAll(".label")
            .style("visibility", "visible")
            .style("opacity", "1")
            .style("display", "block");
            
        this.nodeGroup.selectAll(".node")
            .style("visibility", "visible")
            .style("opacity", "1")
            .style("display", "block");
            
        this.linkGroup.selectAll(".link")
            .style("visibility", "visible")
            .style("opacity", "1")
            .style("display", "block");

        if (this.lastRoot) {
            this.nodeGroup.selectAll(".node")
                .attr("transform", d => {
                    const radius = this.config.showBranchLengths ? 
                        (d.radius || d.y) : d.y;
                    const angle = (d.x - 90) / 180 * Math.PI;
                    const x = radius * Math.cos(angle);
                    const y = radius * Math.sin(angle);
                    return `translate(${x}, ${y})`;
                });

            this.labelGroup.selectAll(".label")
                .attr("transform", d => this.labelTransform(d));
        }
    }

    isRendered() {
        return this.lastRoot && 
               this.nodeGroup.selectAll(".node").size() > 0 &&
               this.labelGroup.selectAll(".label").size() > 0 &&
               this.linkGroup.selectAll(".link").size() > 0;
    }

    safeUpdate(selectedProteins) {
        if (this.isRendered()) {
            this.updateSelection(selectedProteins);
        } else {
            console.warn('Dendrogram not fully rendered, forcing re-render...');
            this.forceRerender();
        }
    }

    // Search-bar matches: a halo behind each matching leaf and a bold label.
    // codes = matching protein codes (best first when ranked), or null = no active search.
    setSearchHits(codes, ranked = false) {
        this.searchHits = codes ? new Map(codes.map((code, i) => [code, i])) : null;
        this.searchRanked = ranked;
        this.renderSearchHits();
    }

    renderSearchHits() {
        if (!this.hitGroup) return;
        const hits = this.searchHits;
        const isHit = d => !!(hits && d && d.data && hits.has(d.data.name));
        const isBest = d => this.searchRanked && isHit(d) && hits.get(d.data.name) === 0;
        const leaves = hits && this.lastRoot ? this.lastRoot.leaves().filter(isHit) : [];

        const halos = this.hitGroup.selectAll(".search-hit-halo").data(leaves, d => d.data.name);
        halos.exit().remove();
        halos.enter()
            .append("circle")
            .attr("class", "search-hit-halo")
            .merge(halos)
            .classed("best", isBest)
            .attr("r", d => this.markSize(isBest(d) ? 11 : 8))
            .attr("transform", d => {
                const radius = this.config.showBranchLengths ? (d.radius || d.y) : d.y;
                const angle = (d.x - 90) / 180 * Math.PI;
                return `translate(${radius * Math.cos(angle)}, ${radius * Math.sin(angle)})`;
            });

        this.labelGroup.selectAll(".label")
            .classed("search-match", isHit);
    }

    // The browse page's groups: a ring in the group's colour around each grouped protein and a bold
    // name, kept while the group exists. window.dendrogramGroupMarks = {code: colour} (set by the
    // page), read here so the rings survive a redraw of the tree.
    renderGroupMarks() {
        if (!this.groupMarkGroup) return;
        const marks = window.dendrogramGroupMarks || {};
        const marked = d => !!(d && d.data && marks[d.data.name]);
        const leaves = this.lastRoot ? this.lastRoot.leaves().filter(marked) : [];
        const rings = this.groupMarkGroup.selectAll(".group-ring").data(leaves, d => d.data.name);
        rings.exit().remove();
        rings.enter()
            .append("circle")
            .attr("class", "group-ring")
            .merge(rings)
            .attr("r", this.markSize(7.5))
            .attr("stroke", d => marks[d.data.name])
            .attr("fill", d => marks[d.data.name])
            .attr("transform", d => {
                const radius = this.config.showBranchLengths ? (d.radius || d.y) : d.y;
                const angle = (d.x - 90) / 180 * Math.PI;
                return `translate(${radius * Math.cos(angle)}, ${radius * Math.sin(angle)})`;
            });
        this.labelGroup.selectAll(".label").classed("group-match", marked);
    }

    // Proteins under the pointer in the catalogue overview charts (catalogue_overview.js) or the
    // clade legend: a halo behind each and a bold label (the other leaves are not faded: that made
    // large groups hard to read). Separate from the search classes, so the search state is back as
    // it was when codes = null clears it.
    setOverviewHighlight(codes) {
        if (!this.hitGroup) return;
        const hits = codes && codes.length ? new Set(codes) : null;
        const isHit = d => !!(hits && d && d.data && hits.has(d.data.name));
        const leaves = hits && this.lastRoot ? this.lastRoot.leaves().filter(isHit) : [];

        const halos = this.hitGroup.selectAll(".overview-halo").data(leaves, d => d.data.name);
        halos.exit().remove();
        halos.enter()
            .append("circle")
            .attr("class", "overview-halo")
            .merge(halos)
            .attr("r", this.markSize(8))
            .attr("transform", d => {
                const radius = this.config.showBranchLengths ? (d.radius || d.y) : d.y;
                const angle = (d.x - 90) / 180 * Math.PI;
                return `translate(${radius * Math.cos(angle)}, ${radius * Math.sin(angle)})`;
            });

        this.labelGroup.selectAll(".label")
            .classed("overview-match", isHit);
    }

    toggleProteinSelection(proteinName) {
        if (window.isGroupedProtein && window.isGroupedProtein(proteinName)) {
            if (window.notifyGroupedProtein) window.notifyGroupedProtein(proteinName);
            return;
        }
        const checkbox = document.querySelector(`input[name="filter_protein_ids"][value="${proteinName}"]`);
        if (checkbox) {
            checkbox.checked = !checkbox.checked;

            // this.showSelectionNotification(proteinName, checkbox.checked);
            
            checkbox.dispatchEvent(new Event('change', { bubbles: true }));
        }
    }


    // Always the whole tree: the zoom is undone for the copy and put back.
    exportSVG() {
        const svgNode = this.svg.node();
        const zoomed = d3.zoomTransform(svgNode);
        this.svg.interrupt().call(this.zoom.transform, d3.zoomIdentity);
        const serializer = new XMLSerializer();
        const svgString = serializer.serializeToString(svgNode);
        this.svg.call(this.zoom.transform, zoomed);
        
        const blob = new Blob([svgString], {type: "image/svg+xml"});
        const url = URL.createObjectURL(blob);
        
        const link = document.createElement("a");
        link.href = url;
        link.download = "dendrogram.svg";
        link.click();
        
        URL.revokeObjectURL(url);
    }

    destroy() {
        if (this.tooltip) {
            this.cancelHideTooltip();
            this.tooltip.remove();
        }
        this.container.selectAll("*").remove();
    }
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = CircularDendrogram;
} else {
    window.CircularDendrogram = CircularDendrogram;
}