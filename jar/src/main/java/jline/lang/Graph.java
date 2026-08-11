/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import java.io.Serializable;
import java.util.List;

/**
 * A lightweight directed-graph view of a {@link Network} topology, returned by
 * {@link Network#getGraph()}. Nodes are the network node names and edges carry
 * the per-class routing probability (weight) between nodes.
 */
public class Graph implements Serializable {

    private static final long serialVersionUID = 1L;

    /**
     * A directed edge between two nodes for a given job class.
     */
    public static class Edge implements Serializable {

        private static final long serialVersionUID = 1L;

        public final String source;
        public final String destination;
        public final String jobClass;
        public final double weight;

        public Edge(String source, String destination, String jobClass, double weight) {
            this.source = source;
            this.destination = destination;
            this.jobClass = jobClass;
            this.weight = weight;
        }

        public String getSource() {
            return source;
        }

        public String getDestination() {
            return destination;
        }

        public String getJobClass() {
            return jobClass;
        }

        public double getWeight() {
            return weight;
        }

        @Override
        public String toString() {
            return source + " -> " + destination
                + " [class=" + jobClass + ", weight=" + weight + "]";
        }
    }

    private final List<String> nodes;
    private final List<Edge> edges;

    public Graph(List<String> nodes, List<Edge> edges) {
        this.nodes = nodes;
        this.edges = edges;
    }

    /**
     * Returns the list of node names in the graph.
     */
    public List<String> getNodes() {
        return nodes;
    }

    /**
     * Returns the list of directed edges in the graph.
     */
    public List<Edge> getEdges() {
        return edges;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Graph with ").append(nodes.size()).append(" nodes and ")
          .append(edges.size()).append(" edges\n");
        sb.append("Nodes: ").append(nodes).append('\n');
        sb.append("Edges:\n");
        for (Edge e : edges) {
            sb.append("  ").append(e).append('\n');
        }
        return sb.toString();
    }
}
