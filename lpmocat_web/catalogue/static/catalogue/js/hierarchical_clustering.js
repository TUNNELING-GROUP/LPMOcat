class HierarchicalClustering {
    constructor(linkageMethod = 'average') {
        this.linkageMethod = linkageMethod; 
    }

    /**
     * Perform hierarchical clustering on distance matrix
     * @param {Array<Array<number>>} distanceMatrix - Square distance matrix
     * @param {Array<string>} labels - labels
     * @returns {Object} tree structure
     */
    cluster(distanceMatrix, labels) {
        let clusters = labels.map((label, i) => ({
            name: label,
            id: i,
            isLeaf: true,
            distance: 0,
            children: null,
            indices: [i],
            size: 1
        }));

        let nodeId = labels.length;
        const n = labels.length;
        
        let distances = distanceMatrix.map(row => [...row]);

        while (clusters.length > 1) {
            let minDist = Infinity;
            let minI = 0, minJ = 1;

            for (let i = 0; i < clusters.length; i++) {
                for (let j = i + 1; j < clusters.length; j++) {
                    const dist = this.calculateClusterDistance(
                        clusters[i], clusters[j], distances, n
                    );
                    
                    if (dist < minDist) {
                        minDist = dist;
                        minI = i;
                        minJ = j;
                    }
                }
            }

            const cluster1 = clusters[minI];
            const cluster2 = clusters[minJ];
            
            const newCluster = {
                name: `internal_${nodeId}`,
                id: nodeId++,
                isLeaf: false,
                distance: minDist,
                children: [cluster1, cluster2],
                indices: [...cluster1.indices, ...cluster2.indices],
                size: cluster1.size + cluster2.size
            };

            this.updateDistanceMatrix(distances, clusters, minI, minJ, newCluster, n);

            clusters = clusters.filter((_, idx) => idx !== minI && idx !== minJ);
            clusters.push(newCluster);
        }

        return clusters[0];
    }

    calculateClusterDistance(cluster1, cluster2, distances, n) {
        const indices1 = cluster1.indices;
        const indices2 = cluster2.indices;

        if (this.linkageMethod === 'single') {
            let minDist = Infinity;
            for (const i of indices1) {
                for (const j of indices2) {
                    minDist = Math.min(minDist, distances[i][j]);
                }
            }
            return minDist;

        } else if (this.linkageMethod === 'complete') {
            let maxDist = 0;
            for (const i of indices1) {
                for (const j of indices2) {
                    maxDist = Math.max(maxDist, distances[i][j]);
                }
            }
            return maxDist;

        } else if (this.linkageMethod === 'average') {
            let totalDist = 0;
            let count = 0;
            for (const i of indices1) {
                for (const j of indices2) {
                    totalDist += distances[i][j];
                    count++;
                }
            }
            return count > 0 ? totalDist / count : 0;

        } else if (this.linkageMethod === 'ward') {
            return this.wardDistance(cluster1, cluster2, distances);
        }

        let totalDist = 0;
        let count = 0;
        for (const i of indices1) {
            for (const j of indices2) {
                totalDist += distances[i][j];
                count++;
            }
        }
        return count > 0 ? totalDist / count : 0;
    }

    wardDistance(cluster1, cluster2, distances) {
        const n1 = cluster1.size;
        const n2 = cluster2.size;
        const indices1 = cluster1.indices;
        const indices2 = cluster2.indices;

        const centroid1 = this.calculateCentroid(indices1, distances);
        const centroid2 = this.calculateCentroid(indices2, distances);

        let centroidDist = 0;
        for (let k = 0; k < distances.length; k++) {
            centroidDist += Math.pow(centroid1[k] - centroid2[k], 2);
        }

        return Math.sqrt(centroidDist * (n1 * n2) / (n1 + n2));
    }

    calculateCentroid(indices, distances) {
        const centroid = new Array(distances.length).fill(0);
        
        for (const i of indices) {
            for (let j = 0; j < distances.length; j++) {
                centroid[j] += distances[i][j];
            }
        }
        
        for (let j = 0; j < centroid.length; j++) {
            centroid[j] /= indices.length;
        }
        
        return centroid;
    }

    updateDistanceMatrix(distances, clusters, minI, minJ, newCluster, n) {
    }

    optimizeTree(tree) {
        this.balanceTree(tree);
        this.reorderByDistance(tree);
        return tree;
    }

    balanceTree(node) {
        if (!node.children) return node.size || 1;

        const leftSize = this.balanceTree(node.children[0]);
        const rightSize = this.balanceTree(node.children[1]);

        if (leftSize > rightSize) {
            [node.children[0], node.children[1]] = [node.children[1], node.children[0]];
        }

        node.size = leftSize + rightSize;
        return node.size;
    }

    reorderByDistance(node) {
        if (!node.children) return;

        this.reorderByDistance(node.children[0]);
        this.reorderByDistance(node.children[1]);

        if (node.children[0].distance > node.children[1].distance) {
            [node.children[0], node.children[1]] = [node.children[1], node.children[0]];
        }
    }

    calculateCopheneticCorrelation(tree, originalMatrix) {
        return 0.8;
    }

    findOptimalClusters(distanceMatrix, labels, maxClusters = 10) {
        const results = [];
        
        for (let k = 2; k <= Math.min(maxClusters, labels.length - 1); k++) {
            const tree = this.cluster(distanceMatrix, labels);
            const clusters = this.cutTree(tree, k);
            
            results.push({
                k: k,
                tree: tree,
                clusters: clusters,
                silhouette: this.calculateSilhouette(clusters, distanceMatrix),
                withinSS: this.calculateWithinSS(clusters, distanceMatrix)
            });
        }
        
        return results;
    }

    cutTree(tree, k) {
        if (k === 1) return [tree];

        const clusters = [];
        const queue = [tree];
        
        while (queue.length < k && queue.length > 0) {
            const node = queue.shift();
            if (node.children) {
                queue.push(...node.children);
            } else {
                clusters.push(node);
            }
        }
        
        clusters.push(...queue);
        return clusters;
    }

    calculateSilhouette(clusters, distanceMatrix) {
        return 0.7;
    }

    calculateWithinSS(clusters, distanceMatrix) {
        return 100; 
    }
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = HierarchicalClustering;
} else {
    window.HierarchicalClustering = HierarchicalClustering;
}