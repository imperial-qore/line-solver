/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.wf;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import jline.GlobalConstants;
import jline.util.Pair;
import jline.util.matrix.Matrix;

/**
 * Workflow Loop Pattern Detection.
 *
 * Implements algorithms for detecting loop and iterative patterns in workflow
 * traces.
 */
public final class Wf_loop_detector {
    private Wf_loop_detector() {}

    public static List<Integer> detectLoops(Matrix linkMatrix,
                                            List<Integer> serviceNodes,
                                            List<Integer> routerNodes) {
        return detectLoops(linkMatrix, serviceNodes, routerNodes, Collections.<Integer>emptyList());
    }

    public static List<Integer> detectLoops(Matrix linkMatrix,
                                            List<Integer> serviceNodes,
                                            List<Integer> routerNodes,
                                            List<Integer> joinNodes) {
        List<Integer> loopNodes = new ArrayList<Integer>();
        Set<Integer> routerSet = new HashSet<Integer>(routerNodes);

        Map<Integer, List<Pair<Integer, Double>>> adjacency = buildAdjacencyMap(linkMatrix);

        for (Integer serviceNode : serviceNodes) {
            if (isInSimpleLoop(serviceNode, adjacency, routerSet)) {
                loopNodes.add(serviceNode);
            }
        }

        if (!joinNodes.isEmpty()) {
            List<Integer> complexLoops = findComplexLoops(linkMatrix, serviceNodes, routerNodes, joinNodes);
            loopNodes.addAll(complexLoops);
        }

        // distinct
        LinkedHashSet<Integer> seen = new LinkedHashSet<Integer>(loopNodes);
        return new ArrayList<Integer>(seen);
    }

    private static Map<Integer, List<Pair<Integer, Double>>> buildAdjacencyMap(Matrix linkMatrix) {
        Map<Integer, List<Pair<Integer, Double>>> adjacency = new HashMap<Integer, List<Pair<Integer, Double>>>();

        for (int i = 0; i < linkMatrix.getNumRows(); i++) {
            int start = (int) linkMatrix.get(i, 0);
            int end = (int) linkMatrix.get(i, 1);
            double prob = linkMatrix.get(i, 2);
            List<Pair<Integer, Double>> list = adjacency.get(start);
            if (list == null) {
                list = new ArrayList<Pair<Integer, Double>>();
                adjacency.put(start, list);
            }
            list.add(new Pair<Integer, Double>(end, prob));
        }
        return adjacency;
    }

    private static boolean isInSimpleLoop(int serviceNode,
                                          Map<Integer, List<Pair<Integer, Double>>> adjacency,
                                          Set<Integer> routerSet) {
        List<Pair<Integer, Double>> neighbors = adjacency.get(serviceNode);
        if (neighbors == null) return false;

        for (Pair<Integer, Double> p : neighbors) {
            int routerNode = p.getLeft();
            if (routerSet.contains(routerNode)) {
                List<Pair<Integer, Double>> routerNeighbors = adjacency.get(routerNode);
                if (routerNeighbors == null) continue;
                for (Pair<Integer, Double> q : routerNeighbors) {
                    if (q.getLeft() == serviceNode) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    private static List<Integer> findComplexLoops(Matrix linkMatrix,
                                                  List<Integer> serviceNodes,
                                                  List<Integer> routerNodes,
                                                  List<Integer> joinNodes) {
        List<Integer> complexLoops = new ArrayList<Integer>();
        Set<Integer> serviceSet = new HashSet<Integer>(serviceNodes);
        Set<Integer> routerSet = new HashSet<Integer>(routerNodes);
        Set<Integer> joinSet = new HashSet<Integer>(joinNodes);

        Map<Integer, Set<Integer>> graph = buildDirectedGraph(linkMatrix);
        List<List<Integer>> sccs = findStronglyConnectedComponents(graph);

        for (List<Integer> scc : sccs) {
            if (scc.size() > 1) {
                List<Integer> serviceNodesInSCC = new ArrayList<Integer>();
                boolean hasRouterOrJoin = false;
                for (Integer n : scc) {
                    if (serviceSet.contains(n)) serviceNodesInSCC.add(n);
                    if (routerSet.contains(n) || joinSet.contains(n)) hasRouterOrJoin = true;
                }
                if (!serviceNodesInSCC.isEmpty() && hasRouterOrJoin) {
                    complexLoops.addAll(serviceNodesInSCC);
                }
            }
        }
        return complexLoops;
    }

    private static Map<Integer, Set<Integer>> buildDirectedGraph(Matrix linkMatrix) {
        Map<Integer, Set<Integer>> graph = new HashMap<Integer, Set<Integer>>();
        for (int i = 0; i < linkMatrix.getNumRows(); i++) {
            int start = (int) linkMatrix.get(i, 0);
            int end = (int) linkMatrix.get(i, 1);
            Set<Integer> set = graph.get(start);
            if (set == null) {
                set = new HashSet<Integer>();
                graph.put(start, set);
            }
            set.add(end);
        }
        return graph;
    }

    private static List<List<Integer>> findStronglyConnectedComponents(Map<Integer, Set<Integer>> graph) {
        final List<List<Integer>> sccs = new ArrayList<List<Integer>>();
        final Map<Integer, Integer> indices = new HashMap<Integer, Integer>();
        final Map<Integer, Integer> lowLinks = new HashMap<Integer, Integer>();
        final Set<Integer> onStack = new HashSet<Integer>();
        final Deque<Integer> stack = new ArrayDeque<Integer>();
        final int[] indexCounter = new int[]{0};

        for (Integer node : graph.keySet()) {
            if (!indices.containsKey(node)) {
                strongConnect(node, graph, indices, lowLinks, onStack, stack, sccs, indexCounter);
            }
        }
        return sccs;
    }

    private static void strongConnect(int node,
                                      Map<Integer, Set<Integer>> graph,
                                      Map<Integer, Integer> indices,
                                      Map<Integer, Integer> lowLinks,
                                      Set<Integer> onStack,
                                      Deque<Integer> stack,
                                      List<List<Integer>> sccs,
                                      int[] indexCounter) {
        indices.put(node, indexCounter[0]);
        lowLinks.put(node, indexCounter[0]);
        indexCounter[0]++;
        stack.push(node);
        onStack.add(node);

        Set<Integer> neighbors = graph.get(node);
        if (neighbors != null) {
            for (Integer neighbor : neighbors) {
                if (!indices.containsKey(neighbor)) {
                    strongConnect(neighbor, graph, indices, lowLinks, onStack, stack, sccs, indexCounter);
                    lowLinks.put(node, Math.min(lowLinks.get(node), lowLinks.get(neighbor)));
                } else if (onStack.contains(neighbor)) {
                    lowLinks.put(node, Math.min(lowLinks.get(node), indices.get(neighbor)));
                }
            }
        }

        if (lowLinks.get(node).equals(indices.get(node))) {
            List<Integer> scc = new ArrayList<Integer>();
            int w;
            do {
                w = stack.pop();
                onStack.remove(w);
                scc.add(w);
            } while (w != node);
            sccs.add(scc);
        }
    }

    public static double getLoopProbability(int serviceNode,
                                            Matrix linkMatrix,
                                            List<Integer> routerNodes) {
        Set<Integer> routerSet = new HashSet<Integer>(routerNodes);

        for (int i = 0; i < linkMatrix.getNumRows(); i++) {
            int start = (int) linkMatrix.get(i, 0);
            int end = (int) linkMatrix.get(i, 1);

            if (start == serviceNode && routerSet.contains(end)) {
                for (int j = 0; j < linkMatrix.getNumRows(); j++) {
                    int loopStart = (int) linkMatrix.get(j, 0);
                    int loopEnd = (int) linkMatrix.get(j, 1);
                    double loopProb = linkMatrix.get(j, 2);

                    if (loopStart == end && loopEnd == serviceNode) {
                        return loopProb;
                    }
                }
            }
        }
        return 0.0;
    }

    public static boolean validateLoopPattern(int loopNode,
                                              Matrix linkMatrix,
                                              List<Integer> routerNodes) {
        Set<Integer> routerSet = new HashSet<Integer>(routerNodes);
        Map<Integer, List<Pair<Integer, Double>>> adjacency = buildAdjacencyMap(linkMatrix);
        List<Pair<Integer, Double>> neighbors = adjacency.get(loopNode);
        if (neighbors == null) return false;

        for (Pair<Integer, Double> p : neighbors) {
            int routerNode = p.getLeft();
            if (routerSet.contains(routerNode)) {
                List<Pair<Integer, Double>> routerNeighbors = adjacency.get(routerNode);
                if (routerNeighbors == null) continue;
                for (Pair<Integer, Double> q : routerNeighbors) {
                    if (q.getLeft() == loopNode) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    public static double getExpectedLoopIterations(double loopProbability) {
        if (loopProbability >= 1.0) return GlobalConstants.Inf;
        return 1.0 / (1.0 - loopProbability);
    }

    public static Map<String, Object> getLoopStats(List<Integer> loopNodes,
                                                   Matrix linkMatrix,
                                                   List<Integer> routerNodes) {
        Map<String, Object> stats = new HashMap<String, Object>();
        stats.put("numLoops", loopNodes.size());

        List<Double> probabilities = new ArrayList<Double>();
        for (Integer n : loopNodes) {
            probabilities.add(getLoopProbability(n, linkMatrix, routerNodes));
        }

        double avgProb = 0.0;
        double maxProb = 0.0;
        double minProb = 0.0;
        if (!probabilities.isEmpty()) {
            double sum = 0.0;
            maxProb = Double.NEGATIVE_INFINITY;
            minProb = Double.POSITIVE_INFINITY;
            for (Double p : probabilities) {
                sum += p;
                if (p > maxProb) maxProb = p;
                if (p < minProb) minProb = p;
            }
            avgProb = sum / probabilities.size();
        }
        stats.put("avgLoopProbability", avgProb);
        stats.put("maxLoopProbability", maxProb);
        stats.put("minLoopProbability", minProb);

        List<Double> iterations = new ArrayList<Double>();
        for (Double p : probabilities) {
            double it = getExpectedLoopIterations(p);
            if (Double.isFinite(it)) iterations.add(it);
        }
        double avgIt = 0.0;
        double maxIt = 0.0;
        if (!iterations.isEmpty()) {
            double sum = 0.0;
            maxIt = Double.NEGATIVE_INFINITY;
            for (Double v : iterations) {
                sum += v;
                if (v > maxIt) maxIt = v;
            }
            avgIt = sum / iterations.size();
        }
        stats.put("avgExpectedIterations", avgIt);
        stats.put("maxExpectedIterations", maxIt);

        return stats;
    }
}
