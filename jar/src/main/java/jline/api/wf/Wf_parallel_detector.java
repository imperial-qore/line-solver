/**
 * @file Workflow Parallel Pattern Detection
 *
 * @since LINE 3.0
 */
package jline.api.wf;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;

import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class Wf_parallel_detector {
    private Wf_parallel_detector() {}

    /**
     * Detect parallel patterns in a workflow network.
     */
    public static List<List<Integer>> detectParallel(
            Matrix linkMatrix,
            List<Integer> serviceNodes,
            List<Integer> forkNodes,
            List<Integer> joinNodes) {
        List<List<Integer>> parallelPatterns = new ArrayList<List<Integer>>();

        List<Pair<Integer, Integer>> forkJoinPairs = findForkJoinPairs(linkMatrix, forkNodes, joinNodes);

        for (Pair<Integer, Integer> p : forkJoinPairs) {
            int forkNode = p.getLeft();
            int joinNode = p.getRight();
            List<Integer> parallelServices = findParallelServices(linkMatrix, serviceNodes, forkNode, joinNode);
            if (parallelServices.size() > 1) {
                parallelPatterns.add(parallelServices);
            }
        }

        return parallelPatterns;
    }

    private static List<Pair<Integer, Integer>> findForkJoinPairs(
            Matrix linkMatrix, List<Integer> forkNodes, List<Integer> joinNodes) {
        List<Pair<Integer, Integer>> pairs = new ArrayList<Pair<Integer, Integer>>();
        Set<Integer> forkSet = new HashSet<Integer>(forkNodes);
        Set<Integer> joinSet = new HashSet<Integer>(joinNodes);

        Map<Integer, List<Integer>> adjacency = buildAdjacencyMap(linkMatrix);

        for (Integer fork : forkNodes) {
            for (Integer join : joinNodes) {
                if (isValidForkJoinPair(fork, join, adjacency, forkSet, joinSet)) {
                    pairs.add(new Pair<Integer, Integer>(fork, join));
                }
            }
        }

        return pairs;
    }

    private static Map<Integer, List<Integer>> buildAdjacencyMap(Matrix linkMatrix) {
        Map<Integer, List<Integer>> adjacency = new HashMap<Integer, List<Integer>>();

        for (int i = 0; i < linkMatrix.getNumRows(); i++) {
            int start = (int) linkMatrix.get(i, 0);
            int end = (int) linkMatrix.get(i, 1);

            List<Integer> list = adjacency.get(start);
            if (list == null) {
                list = new ArrayList<Integer>();
                adjacency.put(start, list);
            }
            list.add(end);
        }

        return adjacency;
    }

    private static boolean isValidForkJoinPair(
            int fork, int join, Map<Integer, List<Integer>> adjacency,
            Set<Integer> forkSet, Set<Integer> joinSet) {
        Queue<Integer> queue = new LinkedList<Integer>();
        Set<Integer> visited = new HashSet<Integer>();
        Map<Integer, Integer> pathCount = new HashMap<Integer, Integer>();

        queue.offer(fork);
        pathCount.put(fork, 1);

        while (!queue.isEmpty()) {
            int current = queue.poll();
            if (visited.contains(current)) continue;
            visited.add(current);

            List<Integer> neighbors = adjacency.get(current);
            if (neighbors == null) continue;
            for (Integer neighbor : neighbors) {
                if (neighbor == join) {
                    int currentPaths = pathCount.containsKey(current) ? pathCount.get(current) : 0;
                    int joinCount = pathCount.containsKey(join) ? pathCount.get(join) : 0;
                    pathCount.put(join, joinCount + currentPaths);
                } else if (!visited.contains(neighbor) && !forkSet.contains(neighbor) && !joinSet.contains(neighbor)) {
                    queue.offer(neighbor);
                    int currentPaths = pathCount.containsKey(current) ? pathCount.get(current) : 0;
                    int neighborCount = pathCount.containsKey(neighbor) ? pathCount.get(neighbor) : 0;
                    pathCount.put(neighbor, neighborCount + currentPaths);
                }
            }
        }

        Integer joinPaths = pathCount.get(join);
        return joinPaths != null && joinPaths > 1;
    }

    private static List<Integer> findParallelServices(
            Matrix linkMatrix, List<Integer> serviceNodes, int forkNode, int joinNode) {
        List<Integer> parallelServices = new ArrayList<Integer>();
        Set<Integer> serviceSet = new HashSet<Integer>(serviceNodes);

        Set<Integer> reachableFromFork = findReachableNodes(linkMatrix, forkNode, joinNode);
        Set<Integer> canReachJoin = findNodesThatCanReach(linkMatrix, joinNode, forkNode);

        Set<Integer> parallelNodes = new HashSet<Integer>(reachableFromFork);
        parallelNodes.retainAll(canReachJoin);

        for (Integer node : parallelNodes) {
            if (serviceSet.contains(node)) {
                parallelServices.add(node);
            }
        }

        return parallelServices;
    }

    private static Set<Integer> findReachableNodes(Matrix linkMatrix, int startNode, int endNode) {
        Set<Integer> reachable = new HashSet<Integer>();
        Queue<Integer> queue = new LinkedList<Integer>();
        Set<Integer> visited = new HashSet<Integer>();

        queue.offer(startNode);

        while (!queue.isEmpty()) {
            int current = queue.poll();
            if (visited.contains(current) || current == endNode) continue;
            visited.add(current);

            for (int i = 0; i < linkMatrix.getNumRows(); i++) {
                int start = (int) linkMatrix.get(i, 0);
                int end = (int) linkMatrix.get(i, 1);

                if (start == current && end != endNode) {
                    reachable.add(end);
                    queue.offer(end);
                }
            }
        }

        return reachable;
    }

    private static Set<Integer> findNodesThatCanReach(Matrix linkMatrix, int targetNode, int startNode) {
        Set<Integer> canReach = new HashSet<Integer>();

        Map<Integer, List<Integer>> reverseAdj = new HashMap<Integer, List<Integer>>();
        for (int i = 0; i < linkMatrix.getNumRows(); i++) {
            int start = (int) linkMatrix.get(i, 0);
            int end = (int) linkMatrix.get(i, 1);
            List<Integer> list = reverseAdj.get(end);
            if (list == null) {
                list = new ArrayList<Integer>();
                reverseAdj.put(end, list);
            }
            list.add(start);
        }

        Queue<Integer> queue = new LinkedList<Integer>();
        Set<Integer> visited = new HashSet<Integer>();

        queue.offer(targetNode);

        while (!queue.isEmpty()) {
            int current = queue.poll();
            if (visited.contains(current) || current == startNode) continue;
            visited.add(current);

            List<Integer> predecessors = reverseAdj.get(current);
            if (predecessors == null) continue;
            for (Integer pred : predecessors) {
                if (pred != startNode) {
                    canReach.add(pred);
                    queue.offer(pred);
                }
            }
        }

        return canReach;
    }

    public static boolean validateParallelPattern(
            List<Integer> pattern, Matrix linkMatrix,
            List<Integer> forkNodes, List<Integer> joinNodes) {
        if (pattern.size() < 2) return false;

        Integer forkNode = findCommonSource(pattern, linkMatrix, forkNodes);
        Integer joinNode = findCommonTarget(pattern, linkMatrix, joinNodes);

        return forkNode != null && joinNode != null;
    }

    private static Integer findCommonSource(List<Integer> pattern, Matrix linkMatrix, List<Integer> forkNodes) {
        Set<Integer> forkSet = new HashSet<Integer>(forkNodes);
        Set<Integer> sources = new HashSet<Integer>();

        for (Integer node : pattern) {
            for (int i = 0; i < linkMatrix.getNumRows(); i++) {
                int start = (int) linkMatrix.get(i, 0);
                int end = (int) linkMatrix.get(i, 1);

                if (end == node && forkSet.contains(start)) {
                    sources.add(start);
                }
            }
        }

        return sources.size() == 1 ? sources.iterator().next() : null;
    }

    private static Integer findCommonTarget(List<Integer> pattern, Matrix linkMatrix, List<Integer> joinNodes) {
        Set<Integer> joinSet = new HashSet<Integer>(joinNodes);
        Set<Integer> targets = new HashSet<Integer>();

        for (Integer node : pattern) {
            for (int i = 0; i < linkMatrix.getNumRows(); i++) {
                int start = (int) linkMatrix.get(i, 0);
                int end = (int) linkMatrix.get(i, 1);

                if (start == node && joinSet.contains(end)) {
                    targets.add(end);
                }
            }
        }

        return targets.size() == 1 ? targets.iterator().next() : null;
    }

    public static Map<String, Object> getParallelStats(List<List<Integer>> patterns) {
        Map<String, Object> stats = new HashMap<String, Object>();

        stats.put("numPatterns", patterns.size());
        int total = 0;
        int max = 0;
        for (List<Integer> p : patterns) {
            total += p.size();
            if (p.size() > max) max = p.size();
        }
        stats.put("totalParallelNodes", total);
        stats.put("avgParallelism", patterns.isEmpty() ? 0.0 : ((double) total) / patterns.size());
        stats.put("maxParallelism", max);

        return stats;
    }
}
