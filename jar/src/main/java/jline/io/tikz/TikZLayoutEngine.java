/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io.tikz;

import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.nodes.*;
import jline.util.matrix.Matrix;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Computes automatic layout for queueing network visualization.
 * Uses a layered graph drawing approach (Sugiyama-style).
 */
public class TikZLayoutEngine {

    private final Network model;
    private final TikZOptions options;
    private final Map<Node, double[]> positions;
    private List<List<Node>> layers;

    public TikZLayoutEngine(Network model, TikZOptions options) {
        this.model = model;
        this.options = options;
        this.positions = new HashMap<Node, double[]>();
        this.layers = new ArrayList<List<Node>>();
    }

    /**
     * Computes the layout for all nodes.
     */
    public void computeLayout() {
        List<Node> nodes = model.getNodes();
        if (nodes.isEmpty()) {
            return;
        }

        // Build adjacency information
        Map<Node, List<Node>> successors = buildSuccessorMap(nodes);
        Map<Node, List<Node>> predecessors = buildPredecessorMap(nodes, successors);

        // Step 1: Assign nodes to layers using topological sort
        assignLayers(nodes, predecessors, successors);

        // Step 2: Order nodes within layers to minimize crossings
        minimizeCrossings(successors, predecessors);

        // Step 3: Assign x,y coordinates
        assignCoordinates();
    }

    /**
     * Gets the computed position for a node.
     *
     * @param node The node
     * @return Array [x, y] in cm, or null if not computed
     */
    public double[] getPosition(Node node) {
        return positions.get(node);
    }

    /**
     * Gets all computed positions.
     */
    public Map<Node, double[]> getAllPositions() {
        return Collections.unmodifiableMap(positions);
    }

    /**
     * Gets the layers (for debugging/testing).
     */
    public List<List<Node>> getLayers() {
        return Collections.unmodifiableList(layers);
    }

    /**
     * Builds a map from each node to its successors based on the connection matrix.
     */
    private Map<Node, List<Node>> buildSuccessorMap(List<Node> nodes) {
        Map<Node, List<Node>> successors = new HashMap<Node, List<Node>>();
        for (Node node : nodes) {
            successors.put(node, new ArrayList<Node>());
        }

        NetworkStruct sn = model.getStruct();
        if (sn == null || sn.connmatrix == null) {
            return successors;
        }

        Matrix connMatrix = sn.connmatrix;
        int n = nodes.size();

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (connMatrix.get(i, j) > 0) {
                    Node from = nodes.get(i);
                    Node to = nodes.get(j);
                    successors.get(from).add(to);
                }
            }
        }

        return successors;
    }

    /**
     * Builds a map from each node to its predecessors.
     */
    private Map<Node, List<Node>> buildPredecessorMap(List<Node> nodes, Map<Node, List<Node>> successors) {
        Map<Node, List<Node>> predecessors = new HashMap<Node, List<Node>>();
        for (Node node : nodes) {
            predecessors.put(node, new ArrayList<Node>());
        }

        for (Map.Entry<Node, List<Node>> entry : successors.entrySet()) {
            Node from = entry.getKey();
            for (Node to : entry.getValue()) {
                predecessors.get(to).add(from);
            }
        }

        return predecessors;
    }

    /**
     * Assigns nodes to layers using a modified topological sort.
     * Sources and nodes with no predecessors go to layer 0.
     */
    private void assignLayers(List<Node> nodes, Map<Node, List<Node>> predecessors, Map<Node, List<Node>> successors) {
        layers.clear();
        Map<Node, Integer> nodeLayer = new HashMap<Node, Integer>();
        Set<Node> assigned = new HashSet<Node>();

        // Find all sources and nodes with no predecessors
        List<Node> layer0 = new ArrayList<Node>();
        for (Node node : nodes) {
            if (node instanceof Source || predecessors.get(node).isEmpty()) {
                layer0.add(node);
                nodeLayer.put(node, 0);
                assigned.add(node);
            }
        }

        if (layer0.isEmpty() && !nodes.isEmpty()) {
            // No clear starting point - pick first node
            Node first = nodes.get(0);
            layer0.add(first);
            nodeLayer.put(first, 0);
            assigned.add(first);
        }

        layers.add(layer0);

        // BFS-style layer assignment
        int currentLayerIdx = 0;
        while (assigned.size() < nodes.size()) {
            List<Node> currentLayer = layers.get(currentLayerIdx);
            List<Node> nextLayer = new ArrayList<Node>();

            // Find all unassigned nodes whose predecessors are all assigned
            for (Node node : nodes) {
                if (!assigned.contains(node)) {
                    List<Node> preds = predecessors.get(node);
                    boolean allPredsAssigned = true;
                    for (Node pred : preds) {
                        if (!assigned.contains(pred)) {
                            allPredsAssigned = false;
                            break;
                        }
                    }
                    if (allPredsAssigned) {
                        nextLayer.add(node);
                    }
                }
            }

            if (nextLayer.isEmpty()) {
                // Handle cycles - pick an unassigned node
                for (Node node : nodes) {
                    if (!assigned.contains(node)) {
                        nextLayer.add(node);
                        break;
                    }
                }
            }

            if (!nextLayer.isEmpty()) {
                for (Node node : nextLayer) {
                    nodeLayer.put(node, currentLayerIdx + 1);
                    assigned.add(node);
                }
                layers.add(nextLayer);
            }

            currentLayerIdx++;

            // Safety check to avoid infinite loop
            if (currentLayerIdx > nodes.size()) {
                break;
            }
        }

        // Move sinks to the last layer
        if (layers.size() > 1) {
            List<Node> lastLayer = layers.get(layers.size() - 1);
            for (int i = 0; i < layers.size() - 1; i++) {
                List<Node> layer = layers.get(i);
                Iterator<Node> it = layer.iterator();
                while (it.hasNext()) {
                    Node node = it.next();
                    if (node instanceof Sink) {
                        it.remove();
                        if (!lastLayer.contains(node)) {
                            lastLayer.add(node);
                        }
                    }
                }
            }
            // Remove empty layers
            layers = layers.stream().filter(l -> !l.isEmpty()).collect(Collectors.toList());
        }
    }

    /**
     * Minimizes edge crossings by reordering nodes within layers.
     * Uses the barycentric heuristic.
     */
    private void minimizeCrossings(Map<Node, List<Node>> successors, Map<Node, List<Node>> predecessors) {
        // Perform multiple passes to improve ordering
        for (int pass = 0; pass < 4; pass++) {
            // Forward pass: order based on predecessors
            for (int i = 1; i < layers.size(); i++) {
                reorderLayerByBarycenter(layers.get(i), predecessors, layers.get(i - 1));
            }

            // Backward pass: order based on successors
            for (int i = layers.size() - 2; i >= 0; i--) {
                reorderLayerByBarycenter(layers.get(i), successors, layers.get(i + 1));
            }
        }
    }

    /**
     * Reorders a layer based on barycentric positions of connected nodes.
     */
    private void reorderLayerByBarycenter(List<Node> layer, Map<Node, List<Node>> connections, List<Node> referenceLayer) {
        // Compute position of each node in reference layer
        Map<Node, Integer> refPositions = new HashMap<Node, Integer>();
        for (int i = 0; i < referenceLayer.size(); i++) {
            refPositions.put(referenceLayer.get(i), i);
        }

        // Compute barycenter for each node in current layer
        Map<Node, Double> barycenters = new HashMap<Node, Double>();
        for (Node node : layer) {
            List<Node> connected = connections.get(node);
            if (connected == null || connected.isEmpty()) {
                barycenters.put(node, Double.MAX_VALUE);
            } else {
                double sum = 0;
                int count = 0;
                for (Node conn : connected) {
                    Integer pos = refPositions.get(conn);
                    if (pos != null) {
                        sum += pos;
                        count++;
                    }
                }
                barycenters.put(node, count > 0 ? sum / count : Double.MAX_VALUE);
            }
        }

        // Sort layer by barycenter
        Collections.sort(layer, new Comparator<Node>() {
            @Override
            public int compare(Node a, Node b) {
                return Double.compare(barycenters.get(a), barycenters.get(b));
            }
        });
    }

    /**
     * Assigns x,y coordinates to all nodes based on their layer positions.
     */
    private void assignCoordinates() {
        positions.clear();

        for (int layerIdx = 0; layerIdx < layers.size(); layerIdx++) {
            List<Node> layer = layers.get(layerIdx);
            double x = layerIdx * options.getLayerSpacing();

            // Center the layer vertically
            double totalHeight = (layer.size() - 1) * options.getNodeSpacing();
            double startY = totalHeight / 2.0;

            for (int nodeIdx = 0; nodeIdx < layer.size(); nodeIdx++) {
                Node node = layer.get(nodeIdx);
                double y = startY - nodeIdx * options.getNodeSpacing();
                positions.put(node, new double[]{x, y});
            }
        }
    }
}
