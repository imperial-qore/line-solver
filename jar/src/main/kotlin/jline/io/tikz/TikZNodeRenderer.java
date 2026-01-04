/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io.tikz;

import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.*;

/**
 * Renders queueing network nodes as TikZ code.
 */
public class TikZNodeRenderer {

    private final TikZOptions options;

    public TikZNodeRenderer(TikZOptions options) {
        this.options = options;
    }

    /**
     * Returns the TikZ preamble with style definitions.
     */
    public String getPreamble() {
        StringBuilder sb = new StringBuilder();
        sb.append("\\documentclass[tikz,border=").append(options.getBorderPadding()).append("pt]{standalone}\n");
        sb.append("\\usepackage{tikz}\n");
        sb.append("\\usetikzlibrary{arrows.meta,positioning,shapes.geometric,shapes.misc,calc,decorations.pathreplacing}\n");
        sb.append("\n");
        sb.append("\\tikzset{\n");
        sb.append("    % Queue: Rectangle buffer\n");
        sb.append("    queue/.style={\n");
        sb.append("        rectangle,\n");
        sb.append("        draw=black,\n");
        sb.append("        minimum width=1.1cm,\n");
        sb.append("        minimum height=0.72cm,\n");
        sb.append("        fill=white\n");
        sb.append("    },\n");
        sb.append("    % Server circle\n");
        sb.append("    server/.style={\n");
        sb.append("        circle,\n");
        sb.append("        draw=black,\n");
        sb.append("        minimum size=0.6cm,\n");
        sb.append("        fill=white\n");
        sb.append("    },\n");
        sb.append("    % Delay: Vertical rectangle (infinite server)\n");
        sb.append("    delay/.style={\n");
        sb.append("        rectangle,\n");
        sb.append("        draw=black,\n");
        sb.append("        minimum width=0.8cm,\n");
        sb.append("        minimum height=1.5cm,\n");
        sb.append("        fill=green!10\n");
        sb.append("    },\n");
        sb.append("    % Source: Small circle (white)\n");
        sb.append("    source/.style={\n");
        sb.append("        circle,\n");
        sb.append("        draw=black,\n");
        sb.append("        minimum size=0.3cm,\n");
        sb.append("        fill=white\n");
        sb.append("    },\n");
        sb.append("    % Sink: Small circle (black)\n");
        sb.append("    sink/.style={\n");
        sb.append("        circle,\n");
        sb.append("        draw=black,\n");
        sb.append("        minimum size=0.3cm,\n");
        sb.append("        fill=black\n");
        sb.append("    },\n");
        sb.append("    % Fork: Diamond\n");
        sb.append("    fork/.style={\n");
        sb.append("        diamond,\n");
        sb.append("        draw=black,\n");
        sb.append("        minimum size=1cm,\n");
        sb.append("        fill=orange!20,\n");
        sb.append("        aspect=1.5\n");
        sb.append("    },\n");
        sb.append("    % Join: Diamond\n");
        sb.append("    joinnode/.style={\n");
        sb.append("        diamond,\n");
        sb.append("        draw=black,\n");
        sb.append("        minimum size=1cm,\n");
        sb.append("        fill=purple!20,\n");
        sb.append("        aspect=1.5\n");
        sb.append("    },\n");
        sb.append("    % Router: Hexagon\n");
        sb.append("    router/.style={\n");
        sb.append("        regular polygon,\n");
        sb.append("        regular polygon sides=6,\n");
        sb.append("        draw=black,\n");
        sb.append("        minimum size=1cm,\n");
        sb.append("        fill=cyan!10\n");
        sb.append("    },\n");
        sb.append("    % ClassSwitch: Trapezium\n");
        sb.append("    classswitch/.style={\n");
        sb.append("        trapezium,\n");
        sb.append("        draw=black,\n");
        sb.append("        trapezium left angle=70,\n");
        sb.append("        trapezium right angle=110,\n");
        sb.append("        minimum width=1.5cm,\n");
        sb.append("        minimum height=0.8cm,\n");
        sb.append("        fill=pink!20\n");
        sb.append("    },\n");
        sb.append("    % Cache: Stacked rectangle\n");
        sb.append("    cache/.style={\n");
        sb.append("        rectangle,\n");
        sb.append("        draw=black,\n");
        sb.append("        minimum width=1.5cm,\n");
        sb.append("        minimum height=1.2cm,\n");
        sb.append("        fill=gray!20\n");
        sb.append("    },\n");
        sb.append("    % Logger: Rectangle with lines\n");
        sb.append("    logger/.style={\n");
        sb.append("        rectangle,\n");
        sb.append("        draw=black,\n");
        sb.append("        minimum width=1.2cm,\n");
        sb.append("        minimum height=0.8cm,\n");
        sb.append("        fill=brown!10\n");
        sb.append("    },\n");
        sb.append("    % Place (Petri net): Circle\n");
        sb.append("    place/.style={\n");
        sb.append("        circle,\n");
        sb.append("        draw=black,\n");
        sb.append("        minimum size=0.8cm,\n");
        sb.append("        fill=white\n");
        sb.append("    },\n");
        sb.append("    % Transition (Petri net): Rectangle\n");
        sb.append("    transition/.style={\n");
        sb.append("        rectangle,\n");
        sb.append("        draw=black,\n");
        sb.append("        minimum width=0.2cm,\n");
        sb.append("        minimum height=1cm,\n");
        sb.append("        fill=black\n");
        sb.append("    },\n");
        sb.append("    % Connection arrow\n");
        sb.append("    conn/.style={\n");
        sb.append("        ->,\n");
        sb.append("        >=Stealth,\n");
        sb.append("        thick\n");
        sb.append("    },\n");
        sb.append("    % Probability label\n");
        sb.append("    problabel/.style={\n");
        sb.append("        font=\\footnotesize,\n");
        sb.append("        midway,\n");
        sb.append("        above,\n");
        sb.append("        sloped\n");
        sb.append("    },\n");
        sb.append("    % Node name label\n");
        sb.append("    nodelabel/.style={\n");
        sb.append("        font=\\small\n");
        sb.append("    }\n");
        sb.append("}\n");
        return sb.toString();
    }

    /**
     * Renders a node at the given position.
     *
     * @param node The node to render
     * @param x    X coordinate (cm)
     * @param y    Y coordinate (cm)
     * @return TikZ code for the node
     */
    public String renderNode(Node node, double x, double y) {
        String nodeId = sanitizeId(node.getName());

        // Check Delay before Queue since Delay extends Queue
        if (node instanceof Delay) {
            return renderDelay((Delay) node, nodeId, x, y);
        } else if (node instanceof Queue) {
            return renderQueue((Queue) node, nodeId, x, y);
        } else if (node instanceof Source) {
            return renderSource((Source) node, nodeId, x, y);
        } else if (node instanceof Sink) {
            return renderSink((Sink) node, nodeId, x, y);
        } else if (node instanceof Fork) {
            return renderFork((Fork) node, nodeId, x, y);
        } else if (node instanceof Join) {
            return renderJoin((Join) node, nodeId, x, y);
        } else if (node instanceof Router) {
            return renderRouter((Router) node, nodeId, x, y);
        } else if (node instanceof ClassSwitch) {
            return renderClassSwitch((ClassSwitch) node, nodeId, x, y);
        } else if (node instanceof Cache) {
            return renderCache((Cache) node, nodeId, x, y);
        } else if (node instanceof Logger) {
            return renderLogger((Logger) node, nodeId, x, y);
        } else if (node instanceof Place) {
            return renderPlace((Place) node, nodeId, x, y);
        } else if (node instanceof Transition) {
            return renderTransition((Transition) node, nodeId, x, y);
        } else {
            // Generic node fallback
            return renderGenericNode(node, nodeId, x, y);
        }
    }

    private String renderQueue(Queue queue, String nodeId, double x, double y) {
        StringBuilder sb = new StringBuilder();
        // Buffer rectangle
        sb.append(String.format("\\node[queue] (%s) at (%.2f,%.2f) {};\n", nodeId, x, y));

        // Node name label
        if (options.isShowNodeNames()) {
            sb.append(String.format("\\node[nodelabel,above=2pt of %s] {%s};\n", nodeId, escapeLatex(queue.getName())));
        }

        // Scheduling strategy label
        if (options.isShowScheduling()) {
            SchedStrategy sched = queue.getSchedStrategy();
            if (sched != null) {
                sb.append(String.format("\\node[font=\\tiny,below=2pt of %s] {%s};\n", nodeId, sched.name()));
            }
        }

        // Server circle (attached directly to buffer, no gap)
        if (options.isShowServerCount()) {
            int servers = queue.getNumberOfServers();
            String serverLabel = servers == Integer.MAX_VALUE ? "$\\infty$" : String.valueOf(servers);
            // Position server so its west edge touches the buffer's east edge
            sb.append(String.format("\\node[server,anchor=west] (%s_server) at (%s.east) {%s};\n", nodeId, nodeId, serverLabel));
        }

        return sb.toString();
    }

    private String renderDelay(Delay delay, String nodeId, double x, double y) {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("\\node[delay] (%s) at (%.2f,%.2f) {$\\infty$};\n", nodeId, x, y));

        if (options.isShowNodeNames()) {
            sb.append(String.format("\\node[nodelabel,above=2pt of %s] {%s};\n", nodeId, escapeLatex(delay.getName())));
        }

        return sb.toString();
    }

    private String renderSource(Source source, String nodeId, double x, double y) {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("\\node[source] (%s) at (%.2f,%.2f) {};\n", nodeId, x, y));

        if (options.isShowNodeNames()) {
            sb.append(String.format("\\node[nodelabel,above=2pt of %s] {%s};\n", nodeId, escapeLatex(source.getName())));
        }

        // Incoming arrow to indicate external arrivals
        sb.append(String.format("\\draw[conn] ([xshift=-0.6cm]%s.west) -- (%s.west);\n", nodeId, nodeId));

        return sb.toString();
    }

    private String renderSink(Sink sink, String nodeId, double x, double y) {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("\\node[sink] (%s) at (%.2f,%.2f) {};\n", nodeId, x, y));

        if (options.isShowNodeNames()) {
            sb.append(String.format("\\node[nodelabel,above=2pt of %s] {%s};\n", nodeId, escapeLatex(sink.getName())));
        }

        // Outgoing arrow to indicate departures
        sb.append(String.format("\\draw[conn] (%s.east) -- ([xshift=0.6cm]%s.east);\n", nodeId, nodeId));

        return sb.toString();
    }

    private String renderFork(Fork fork, String nodeId, double x, double y) {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("\\node[fork] (%s) at (%.2f,%.2f) {};\n", nodeId, x, y));

        if (options.isShowNodeNames()) {
            sb.append(String.format("\\node[nodelabel,above=2pt of %s] {%s};\n", nodeId, escapeLatex(fork.getName())));
        }

        return sb.toString();
    }

    private String renderJoin(Join join, String nodeId, double x, double y) {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("\\node[joinnode] (%s) at (%.2f,%.2f) {};\n", nodeId, x, y));

        if (options.isShowNodeNames()) {
            sb.append(String.format("\\node[nodelabel,above=2pt of %s] {%s};\n", nodeId, escapeLatex(join.getName())));
        }

        return sb.toString();
    }

    private String renderRouter(Router router, String nodeId, double x, double y) {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("\\node[router] (%s) at (%.2f,%.2f) {};\n", nodeId, x, y));

        if (options.isShowNodeNames()) {
            sb.append(String.format("\\node[nodelabel,above=2pt of %s] {%s};\n", nodeId, escapeLatex(router.getName())));
        }

        return sb.toString();
    }

    private String renderClassSwitch(ClassSwitch cs, String nodeId, double x, double y) {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("\\node[classswitch] (%s) at (%.2f,%.2f) {};\n", nodeId, x, y));

        if (options.isShowNodeNames()) {
            sb.append(String.format("\\node[nodelabel,above=2pt of %s] {%s};\n", nodeId, escapeLatex(cs.getName())));
        }

        return sb.toString();
    }

    private String renderCache(Cache cache, String nodeId, double x, double y) {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("\\node[cache] (%s) at (%.2f,%.2f) {};\n", nodeId, x, y));

        if (options.isShowNodeNames()) {
            sb.append(String.format("\\node[nodelabel,above=2pt of %s] {%s};\n", nodeId, escapeLatex(cache.getName())));
        }

        // Draw cache level lines (three horizontal lines representing cache blocks)
        sb.append(String.format("\\draw ([yshift=-0.3cm]%s.north west) -- ([yshift=-0.3cm]%s.north east);\n", nodeId, nodeId));
        sb.append(String.format("\\draw (%s.west) -- (%s.east);\n", nodeId, nodeId));  // Middle line
        sb.append(String.format("\\draw ([yshift=0.3cm]%s.south west) -- ([yshift=0.3cm]%s.south east);\n", nodeId, nodeId));

        return sb.toString();
    }

    private String renderLogger(Logger logger, String nodeId, double x, double y) {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("\\node[logger] (%s) at (%.2f,%.2f) {};\n", nodeId, x, y));

        if (options.isShowNodeNames()) {
            sb.append(String.format("\\node[nodelabel,above=2pt of %s] {%s};\n", nodeId, escapeLatex(logger.getName())));
        }

        return sb.toString();
    }

    private String renderPlace(Place place, String nodeId, double x, double y) {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("\\node[place] (%s) at (%.2f,%.2f) {};\n", nodeId, x, y));

        if (options.isShowNodeNames()) {
            sb.append(String.format("\\node[nodelabel,above=2pt of %s] {%s};\n", nodeId, escapeLatex(place.getName())));
        }

        return sb.toString();
    }

    private String renderTransition(Transition transition, String nodeId, double x, double y) {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("\\node[transition] (%s) at (%.2f,%.2f) {};\n", nodeId, x, y));

        if (options.isShowNodeNames()) {
            sb.append(String.format("\\node[nodelabel,above=2pt of %s] {%s};\n", nodeId, escapeLatex(transition.getName())));
        }

        return sb.toString();
    }

    private String renderGenericNode(Node node, String nodeId, double x, double y) {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("\\node[draw,circle,minimum size=0.8cm] (%s) at (%.2f,%.2f) {};\n", nodeId, x, y));

        if (options.isShowNodeNames()) {
            sb.append(String.format("\\node[nodelabel,above=2pt of %s] {%s};\n", nodeId, escapeLatex(node.getName())));
        }

        return sb.toString();
    }

    /**
     * Renders a connection between two nodes.
     *
     * @param fromNode Source node
     * @param toNode   Destination node
     * @param prob     Routing probability (or NaN to hide)
     * @return TikZ code for the connection
     */
    public String renderConnection(Node fromNode, Node toNode, double prob) {
        String fromId = sanitizeId(fromNode.getName());
        String toId = sanitizeId(toNode.getName());

        StringBuilder sb = new StringBuilder();

        // Determine anchor points based on node types
        String fromAnchor = getOutputAnchor(fromNode);
        String toAnchor = getInputAnchor(toNode);

        if (options.isShowRoutingProb() && !Double.isNaN(prob) && prob >= options.getMinProbToShow() && prob < 1.0 - options.getMinProbToShow()) {
            // Show probability label (hide if essentially 0 or 1)
            sb.append(String.format("\\draw[conn] (%s%s) -- node[problabel] {%.2f} (%s%s);\n",
                    fromId, fromAnchor, prob, toId, toAnchor));
        } else {
            sb.append(String.format("\\draw[conn] (%s%s) -- (%s%s);\n",
                    fromId, fromAnchor, toId, toAnchor));
        }

        return sb.toString();
    }

    /**
     * Renders a curved connection for parallel paths (e.g., fork-join).
     */
    public String renderCurvedConnection(Node fromNode, Node toNode, double prob, double bendAngle) {
        String fromId = sanitizeId(fromNode.getName());
        String toId = sanitizeId(toNode.getName());

        // Get proper anchor points
        String fromAnchor = getOutputAnchor(fromNode);
        String toAnchor = getInputAnchor(toNode);

        StringBuilder sb = new StringBuilder();
        String bendDir = bendAngle > 0 ? "bend left" : "bend right";
        double absBend = Math.abs(bendAngle);

        if (options.isShowRoutingProb() && !Double.isNaN(prob) && prob >= options.getMinProbToShow() && prob < 1.0 - options.getMinProbToShow()) {
            sb.append(String.format("\\draw[conn] (%s%s) to[%s=%.0f] node[problabel] {%.2f} (%s%s);\n",
                    fromId, fromAnchor, bendDir, absBend, prob, toId, toAnchor));
        } else {
            sb.append(String.format("\\draw[conn] (%s%s) to[%s=%.0f] (%s%s);\n",
                    fromId, fromAnchor, bendDir, absBend, toId, toAnchor));
        }

        return sb.toString();
    }

    /**
     * Gets the output anchor point for a node type.
     * For Queue nodes, returns the server node suffix.
     */
    private String getOutputAnchor(Node node) {
        if (node instanceof Queue && options.isShowServerCount()) {
            return "_server.east";  // Connect from server circle
        }
        return ".east";
    }

    /**
     * Gets the input anchor point for a node type.
     */
    private String getInputAnchor(Node node) {
        return ".west";
    }

    /**
     * Sanitizes a node name for use as a TikZ node identifier.
     */
    private String sanitizeId(String name) {
        // Replace spaces and special characters with underscores
        return name.replaceAll("[^a-zA-Z0-9]", "_");
    }

    /**
     * Escapes special LaTeX characters in text.
     */
    private String escapeLatex(String text) {
        return text
                .replace("\\", "\\textbackslash{}")
                .replace("_", "\\_")
                .replace("&", "\\&")
                .replace("%", "\\%")
                .replace("$", "\\$")
                .replace("#", "\\#")
                .replace("{", "\\{")
                .replace("}", "\\}")
                .replace("~", "\\textasciitilde{}")
                .replace("^", "\\textasciicircum{}");
    }
}
