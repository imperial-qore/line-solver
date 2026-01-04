/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io.tikz;

import jline.lang.nodes.Node;

import java.util.*;

/**
 * Routes edges around nodes to achieve a planar graph layout.
 *
 * This is a general-purpose edge routing algorithm for queueing network diagrams.
 * It handles:
 * - Forward edges: routes around any intermediate nodes
 * - Backward edges (feedback loops): routes around the entire network
 * - Self-loops: routes above/below the node
 * - Vertical edges: routes with horizontal offset
 * - Multiple edges to same target: staggers entry points
 *
 * The algorithm ensures a planar layout where edges don't cross nodes
 * and parallel edges are visually separated.
 */
public class TikZEdgeRouter {

    private final Map<Node, double[]> positions;
    private final TikZOptions options;

    // Node bounding box sizes (half-width, half-height)
    private static final double NODE_HALF_WIDTH = 1.2;  // cm
    private static final double NODE_HALF_HEIGHT = 0.9; // cm
    private static final double MARGIN = 0.25;          // extra margin around nodes (reduced by 50%)
    private static final double ROUTE_SPACING = 0.5;    // spacing between parallel routes (reduced by 50%)
    private static final double ENTRY_OFFSET = 0.2;     // vertical offset for multiple entries to same node (reduced by 50%)

    // Track used routing channels to avoid overlapping edges
    private final Set<String> usedChannels = new HashSet<String>();
    private int backwardEdgeCount = 0;

    // Track edges going to the same target for offset calculation
    private final Map<Node, Integer> targetEdgeCounts = new HashMap<Node, Integer>();

    // Global bounds of all nodes
    private double globalMinX = Double.MAX_VALUE;
    private double globalMaxX = Double.MIN_VALUE;
    private double globalMinY = Double.MAX_VALUE;
    private double globalMaxY = Double.MIN_VALUE;

    public TikZEdgeRouter(Map<Node, double[]> positions, TikZOptions options) {
        this.positions = positions;
        this.options = options;
        computeGlobalBounds();
    }

    private void computeGlobalBounds() {
        for (double[] pos : positions.values()) {
            globalMinX = Math.min(globalMinX, pos[0] - NODE_HALF_WIDTH);
            globalMaxX = Math.max(globalMaxX, pos[0] + NODE_HALF_WIDTH);
            globalMinY = Math.min(globalMinY, pos[1] - NODE_HALF_HEIGHT);
            globalMaxY = Math.max(globalMaxY, pos[1] + NODE_HALF_HEIGHT);
        }
    }

    /**
     * Computes waypoints for an edge to avoid overlapping nodes.
     * Handles forward edges, backward edges, self-loops, and vertical edges.
     *
     * @param from Source node
     * @param to Target node
     * @param allNodes All nodes in the network
     * @return List of waypoints (x, y pairs). Empty if straight line is OK.
     */
    public List<double[]> computeWaypoints(Node from, Node to, List<Node> allNodes) {
        double[] fromPos = positions.get(from);
        double[] toPos = positions.get(to);

        if (fromPos == null || toPos == null) {
            return Collections.emptyList();
        }

        // Handle self-loops (edge from node to itself)
        if (from == to) {
            return routeSelfLoop(fromPos);
        }

        // Check if this is a backward edge (feedback loop)
        boolean isBackward = toPos[0] < fromPos[0] - NODE_HALF_WIDTH;

        if (isBackward) {
            // Route backward edges around the entire graph
            return routeBackwardEdge(fromPos, toPos, allNodes);
        }

        // Check if this is a vertical edge (same X position)
        boolean isVertical = Math.abs(toPos[0] - fromPos[0]) < NODE_HALF_WIDTH;

        if (isVertical) {
            // Route vertical edges with horizontal offset to avoid overlap
            return routeVerticalEdge(fromPos, toPos, allNodes);
        }

        // For forward edges, check for obstacles
        List<Node> obstacles = findObstacles(from, to, fromPos, toPos, allNodes);

        if (obstacles.isEmpty()) {
            return Collections.emptyList(); // Straight line is fine
        }

        // Route forward edge around obstacles
        return routeForwardEdge(fromPos, toPos, obstacles);
    }

    /**
     * Routes a self-loop (edge from a node to itself).
     * Returns empty waypoints - self-loops use special TikZ rendering.
     */
    private List<double[]> routeSelfLoop(double[] nodePos) {
        // Return empty - self-loops are rendered with special TikZ syntax
        return Collections.emptyList();
    }

    /**
     * Checks if this edge is a self-loop.
     */
    public boolean isSelfLoop(Node from, Node to) {
        return from == to;
    }

    /**
     * Routes a vertical edge (nodes at same X position).
     */
    private List<double[]> routeVerticalEdge(double[] fromPos, double[] toPos, List<Node> allNodes) {
        List<double[]> waypoints = new ArrayList<double[]>();

        // Route with horizontal offset
        double offsetX = NODE_HALF_WIDTH + MARGIN + 0.5;
        double midY = (fromPos[1] + toPos[1]) / 2.0;

        waypoints.add(new double[]{fromPos[0] + offsetX, fromPos[1]});
        waypoints.add(new double[]{toPos[0] + offsetX, toPos[1]});

        return waypoints;
    }

    /**
     * Routes a backward edge (feedback loop) around the entire network.
     * Uses orthogonal routing: right, down/up, left, down/up to target.
     * Each backward edge uses a unique channel and entry point to avoid overlaps.
     *
     * To ensure a planar layout (no edge crossings), edges are routed based on
     * the source node's Y position relative to the target:
     * - If source is above target, route above the network
     * - If source is below target, route below the network
     * This prevents crossings when multiple nodes at different Y positions
     * have backward edges to the same target.
     */
    private List<double[]> routeBackwardEdge(double[] fromPos, double[] toPos, List<Node> allNodes) {
        List<double[]> waypoints = new ArrayList<double[]>();

        backwardEdgeCount++;

        // Route based on source Y position relative to target to ensure planarity:
        // - Sources above target route above the network
        // - Sources below target route below the network
        // This prevents crossing when multiple backward edges go to the same target
        boolean routeAbove = fromPos[1] >= toPos[1];

        // Each backward edge gets its own unique channel depth
        int channelIndex = backwardEdgeCount;

        double routeY;
        if (routeAbove) {
            routeY = globalMaxY + MARGIN + (ROUTE_SPACING * channelIndex);
        } else {
            routeY = globalMinY - MARGIN - (ROUTE_SPACING * channelIndex);
        }

        // Exit point: each edge exits at a slightly different X to avoid overlap at source
        double exitXOffset = (backwardEdgeCount - 1) * 0.3;
        double exitX = fromPos[0] + NODE_HALF_WIDTH + MARGIN + exitXOffset;

        // Entry point: each edge enters at a slightly different X to spread them out
        double entryXOffset = (backwardEdgeCount - 1) * 0.3;
        double entryX = toPos[0] - NODE_HALF_WIDTH - MARGIN - entryXOffset;

        // Target Y: offset each edge vertically so they arrive at different points
        double targetYOffset = (backwardEdgeCount - 1) * ENTRY_OFFSET;
        if (!routeAbove) {
            targetYOffset = -targetYOffset; // Edges from below arrive at lower points
        }
        double targetY = toPos[1] + targetYOffset;

        // Waypoint 1: Exit horizontally from source
        waypoints.add(new double[]{exitX, fromPos[1]});

        // Waypoint 2: Go to routing channel (vertical move)
        waypoints.add(new double[]{exitX, routeY});

        // Waypoint 3: Travel horizontally along the routing channel
        waypoints.add(new double[]{entryX, routeY});

        // Waypoint 4: Go up/down toward target at offset Y
        waypoints.add(new double[]{entryX, targetY});

        return waypoints;
    }

    /**
     * Routes a forward edge around obstacles.
     */
    private List<double[]> routeForwardEdge(double[] fromPos, double[] toPos, List<Node> obstacles) {
        List<double[]> waypoints = new ArrayList<double[]>();

        // Determine if we should go above or below the obstacles
        double avgY = 0;
        for (Node obs : obstacles) {
            avgY += positions.get(obs)[1];
        }
        avgY /= obstacles.size();

        double midY = (fromPos[1] + toPos[1]) / 2.0;
        boolean routeAbove = midY >= avgY;

        // Compute vertical offset (reduced by 50% for tighter layout)
        double offset = (NODE_HALF_HEIGHT + MARGIN + 0.15) * (routeAbove ? 1 : -1);
        double routeY = avgY + offset;

        // Find the x-range of obstacles
        double minX = Double.MAX_VALUE;
        double maxX = Double.MIN_VALUE;
        for (Node obs : obstacles) {
            double[] pos = positions.get(obs);
            minX = Math.min(minX, pos[0] - NODE_HALF_WIDTH - MARGIN);
            maxX = Math.max(maxX, pos[0] + NODE_HALF_WIDTH + MARGIN);
        }

        // Create orthogonal waypoints
        // Waypoint 1: Move to entry x, then up/down
        double entryX = Math.max(fromPos[0] + 0.3, minX - 0.3);
        waypoints.add(new double[]{entryX, routeY});

        // Waypoint 2: Move horizontally past obstacles
        double exitX = Math.min(toPos[0] - 0.3, maxX + 0.3);
        if (exitX > entryX + 0.1) {
            waypoints.add(new double[]{exitX, routeY});
        }

        return waypoints;
    }

    /**
     * Finds nodes that lie on the path between from and to.
     */
    private List<Node> findObstacles(Node from, Node to, double[] fromPos, double[] toPos, List<Node> allNodes) {
        List<Node> obstacles = new ArrayList<Node>();

        for (Node node : allNodes) {
            if (node == from || node == to) {
                continue;
            }

            double[] nodePos = positions.get(node);
            if (nodePos == null) {
                continue;
            }

            if (lineIntersectsNode(fromPos, toPos, nodePos)) {
                obstacles.add(node);
            }
        }

        // Sort obstacles by x position
        obstacles.sort(new Comparator<Node>() {
            @Override
            public int compare(Node a, Node b) {
                double[] posA = positions.get(a);
                double[] posB = positions.get(b);
                return Double.compare(posA[0], posB[0]);
            }
        });

        return obstacles;
    }

    /**
     * Checks if a line segment intersects a node's bounding box.
     */
    private boolean lineIntersectsNode(double[] lineStart, double[] lineEnd, double[] nodeCenter) {
        double halfW = NODE_HALF_WIDTH + MARGIN;
        double halfH = NODE_HALF_HEIGHT + MARGIN;

        double left = nodeCenter[0] - halfW;
        double right = nodeCenter[0] + halfW;
        double bottom = nodeCenter[1] - halfH;
        double top = nodeCenter[1] + halfH;

        double x1 = lineStart[0], y1 = lineStart[1];
        double x2 = lineEnd[0], y2 = lineEnd[1];

        // Quick reject
        if ((x1 < left && x2 < left) || (x1 > right && x2 > right)) {
            return false;
        }
        if ((y1 < bottom && y2 < bottom) || (y1 > top && y2 > top)) {
            return false;
        }

        // Check if either endpoint is inside the box
        if (pointInBox(x1, y1, left, right, bottom, top) ||
            pointInBox(x2, y2, left, right, bottom, top)) {
            return true;
        }

        // Check line intersection with box edges
        return lineIntersectsSegment(x1, y1, x2, y2, left, bottom, left, top) ||
               lineIntersectsSegment(x1, y1, x2, y2, right, bottom, right, top) ||
               lineIntersectsSegment(x1, y1, x2, y2, left, bottom, right, bottom) ||
               lineIntersectsSegment(x1, y1, x2, y2, left, top, right, top);
    }

    private boolean pointInBox(double x, double y, double left, double right, double bottom, double top) {
        return x >= left && x <= right && y >= bottom && y <= top;
    }

    private boolean lineIntersectsSegment(double x1, double y1, double x2, double y2,
                                          double x3, double y3, double x4, double y4) {
        double denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
        if (Math.abs(denom) < 1e-10) {
            return false;
        }

        double t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom;
        double u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / denom;

        return t >= 0 && t <= 1 && u >= 0 && u <= 1;
    }

    /**
     * Generates TikZ code for a routed edge with waypoints.
     */
    public static String renderRoutedEdge(String fromId, String toId,
                                          String fromAnchor, String toAnchor,
                                          List<double[]> waypoints, double prob,
                                          TikZOptions options) {
        StringBuilder sb = new StringBuilder();

        if (waypoints.isEmpty()) {
            // Straight edge
            if (options.isShowRoutingProb() && !Double.isNaN(prob) &&
                prob >= options.getMinProbToShow() && prob < 1.0 - options.getMinProbToShow()) {
                sb.append(String.format("\\draw[conn] (%s%s) -- node[problabel] {%.2f} (%s%s);\n",
                        fromId, fromAnchor, prob, toId, toAnchor));
            } else {
                sb.append(String.format("\\draw[conn] (%s%s) -- (%s%s);\n",
                        fromId, fromAnchor, toId, toAnchor));
            }
        } else {
            // Edge with waypoints - use |- and -| for orthogonal routing
            sb.append(String.format("\\draw[conn] (%s%s)", fromId, fromAnchor));
            for (double[] wp : waypoints) {
                sb.append(String.format(" -- (%.2f,%.2f)", wp[0], wp[1]));
            }
            if (options.isShowRoutingProb() && !Double.isNaN(prob) &&
                prob >= options.getMinProbToShow() && prob < 1.0 - options.getMinProbToShow()) {
                sb.append(String.format(" -- node[problabel] {%.2f} (%s%s);\n", prob, toId, toAnchor));
            } else {
                sb.append(String.format(" -- (%s%s);\n", toId, toAnchor));
            }
        }

        return sb.toString();
    }
}
