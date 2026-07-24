/**
 * @file Workflow Branch Pattern Detection
 *
 * Implements algorithms for detecting branching patterns in workflow traces.
 *
 * @since LINE 3.0
 */
package jline.api.wf;

import jline.util.Pair;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Automatic branch structure detector for workflow networks.
 */
public final class Wf_branch_detector {

    private Wf_branch_detector() {}

    /** Branch pattern with probabilities. */
    public static final class BranchPattern {
        private final List<Integer> branchNodes;
        private final List<Double> probabilities;
        private final Integer forkNode;
        private final Integer joinNode;

        public BranchPattern(List<Integer> branchNodes, List<Double> probabilities, Integer forkNode, Integer joinNode) {
            this.branchNodes = branchNodes;
            this.probabilities = probabilities;
            this.forkNode = forkNode;
            this.joinNode = joinNode;
        }

        public List<Integer> getBranchNodes() { return branchNodes; }
        public List<Double> getProbabilities() { return probabilities; }
        public Integer getForkNode() { return forkNode; }
        public Integer getJoinNode() { return joinNode; }
    }

    public static List<BranchPattern> detectBranches(Matrix linkMatrix, List<Integer> serviceNodes, List<Integer> joinNodes) {
        List<BranchPattern> branchPatterns = new ArrayList<BranchPattern>();
        Set<Integer> serviceSet = new HashSet<Integer>(serviceNodes);
        Set<Integer> joinSet = new HashSet<Integer>(joinNodes);

        Map<Integer, List<Pair<Integer, Double>>> adjacency = buildAdjacencyMap(linkMatrix);
        Map<Integer, List<Integer>> reverseAdjacency = buildReverseAdjacencyMap(linkMatrix);

        List<Integer> branchPoints = findBranchPoints(adjacency, serviceSet, joinSet);
        for (Integer branchPoint : branchPoints) {
            BranchPattern pattern = analyzeBranchPattern(branchPoint, adjacency, reverseAdjacency, serviceSet, joinSet);
            if (pattern != null && pattern.getBranchNodes().size() > 1) {
                branchPatterns.add(pattern);
            }
        }
        return branchPatterns;
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

    private static Map<Integer, List<Integer>> buildReverseAdjacencyMap(Matrix linkMatrix) {
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
        return reverseAdj;
    }

    private static List<Integer> findBranchPoints(Map<Integer, List<Pair<Integer, Double>>> adjacency,
                                                  Set<Integer> serviceSet, Set<Integer> joinSet) {
        List<Integer> branchPoints = new ArrayList<Integer>();
        for (Map.Entry<Integer, List<Pair<Integer, Double>>> e : adjacency.entrySet()) {
            List<Pair<Integer, Double>> neighbors = e.getValue();
            if (neighbors.size() > 1) {
                int serviceTargets = 0;
                for (Pair<Integer, Double> p : neighbors) {
                    if (serviceSet.contains(p.getLeft())) serviceTargets++;
                }
                if (serviceTargets > 1) branchPoints.add(e.getKey());
            }
        }
        return branchPoints;
    }

    private static BranchPattern analyzeBranchPattern(int branchPoint,
                                                     Map<Integer, List<Pair<Integer, Double>>> adjacency,
                                                     Map<Integer, List<Integer>> reverseAdjacency,
                                                     Set<Integer> serviceSet,
                                                     Set<Integer> joinSet) {
        List<Pair<Integer, Double>> neighbors = adjacency.get(branchPoint);
        if (neighbors == null) return null;

        List<Pair<Integer, Double>> branchTargets = new ArrayList<Pair<Integer, Double>>();
        for (Pair<Integer, Double> p : neighbors) {
            if (serviceSet.contains(p.getLeft())) branchTargets.add(p);
        }
        if (branchTargets.size() < 2) return null;

        List<Integer> branchNodeIds = new ArrayList<Integer>();
        for (Pair<Integer, Double> p : branchTargets) branchNodeIds.add(p.getLeft());
        Integer commonJoin = findCommonJoinPoint(branchNodeIds, adjacency, reverseAdjacency, joinSet);

        double totalProb = 0.0;
        for (Pair<Integer, Double> p : branchTargets) totalProb += p.getRight();
        if (Math.abs(totalProb - 1.0) > 0.01) return null;

        List<Integer> nodes = new ArrayList<Integer>();
        List<Double> probs = new ArrayList<Double>();
        for (Pair<Integer, Double> p : branchTargets) {
            nodes.add(p.getLeft());
            probs.add(p.getRight());
        }
        return new BranchPattern(nodes, probs, branchPoint, commonJoin);
    }

    private static Integer findCommonJoinPoint(List<Integer> branchNodes,
                                              Map<Integer, List<Pair<Integer, Double>>> adjacency,
                                              Map<Integer, List<Integer>> reverseAdjacency,
                                              Set<Integer> joinSet) {
        if (branchNodes.isEmpty()) return null;
        Set<Integer> commonReachable = findReachableNodes(branchNodes.get(0), adjacency, joinSet);
        for (int i = 1; i < branchNodes.size(); i++) {
            Set<Integer> r = findReachableNodes(branchNodes.get(i), adjacency, joinSet);
            commonReachable.retainAll(r);
        }
        Set<Integer> joinPoints = new HashSet<Integer>(commonReachable);
        joinPoints.retainAll(joinSet);
        if (!joinPoints.isEmpty()) return joinPoints.iterator().next();
        return commonReachable.isEmpty() ? null : commonReachable.iterator().next();
    }

    private static Set<Integer> findReachableNodes(int startNode,
                                                  Map<Integer, List<Pair<Integer, Double>>> adjacency,
                                                  Set<Integer> stopSet) {
        Set<Integer> reachable = new HashSet<Integer>();
        LinkedList<Integer> queue = new LinkedList<Integer>();
        Set<Integer> visited = new HashSet<Integer>();
        queue.offer(startNode);
        while (!queue.isEmpty()) {
            int current = queue.poll();
            if (visited.contains(current)) continue;
            visited.add(current);
            List<Pair<Integer, Double>> neighbors = adjacency.get(current);
            if (neighbors == null) continue;
            for (Pair<Integer, Double> p : neighbors) {
                int neighbor = p.getLeft();
                reachable.add(neighbor);
                if (!stopSet.contains(neighbor)) queue.offer(neighbor);
            }
        }
        return reachable;
    }

    public static boolean validateBranchPattern(BranchPattern pattern, Matrix linkMatrix) {
        double totalProb = 0.0;
        for (Double p : pattern.getProbabilities()) totalProb += p;
        if (Math.abs(totalProb - 1.0) > 0.01) return false;
        if (pattern.getForkNode() == null) return false;
        Map<Integer, List<Pair<Integer, Double>>> adjacency = buildAdjacencyMap(linkMatrix);
        List<Pair<Integer, Double>> forkNeighbors = adjacency.get(pattern.getForkNode());
        if (forkNeighbors == null) return false;
        Set<Integer> forkTargets = new HashSet<Integer>();
        for (Pair<Integer, Double> p : forkNeighbors) forkTargets.add(p.getLeft());
        for (Integer n : pattern.getBranchNodes()) {
            if (!forkTargets.contains(n)) return false;
        }
        return true;
    }

    public static Map<String, Double> calculateBranchDiversity(BranchPattern pattern) {
        Map<String, Double> metrics = new HashMap<String, Double>();
        List<Double> probs = pattern.getProbabilities();
        int n = probs.size();

        double entropy = 0.0;
        for (Double p : probs) {
            if (p > 0) entropy -= p * Math.log(p);
        }
        metrics.put("entropy", entropy);
        metrics.put("normalizedEntropy", n > 1 ? entropy / Math.log((double) n) : 0.0);

        List<Double> sortedProbs = new ArrayList<Double>(probs);
        java.util.Collections.sort(sortedProbs);
        double sumProbs = 0.0;
        for (Double p : probs) sumProbs += p;
        double gini = 0.0;
        for (int i = 0; i < sortedProbs.size(); i++) {
            gini += (2 * (i + 1) - n - 1) * sortedProbs.get(i);
        }
        gini /= (n - 1) * sumProbs;
        metrics.put("gini", Math.abs(gini));

        double maxProb = 0.0;
        for (Double p : probs) if (p > maxProb) maxProb = p;
        metrics.put("balance", 1.0 / (maxProb == 0.0 ? 1.0 : maxProb));
        return metrics;
    }

    public static Map<String, Object> getBranchStats(List<BranchPattern> patterns) {
        Map<String, Object> stats = new HashMap<String, Object>();
        stats.put("numPatterns", patterns.size());
        int total = 0;
        for (BranchPattern p : patterns) total += p.getBranchNodes().size();
        stats.put("totalBranchNodes", total);
        if (patterns.isEmpty()) {
            stats.put("avgBranches", 0.0);
            stats.put("maxBranches", 0);
            stats.put("minBranches", 0);
            stats.put("avgEntropy", 0.0);
            stats.put("avgBalance", 0.0);
            return stats;
        }
        int max = 0, min = Integer.MAX_VALUE;
        double avg = 0.0;
        for (BranchPattern p : patterns) {
            int sz = p.getBranchNodes().size();
            if (sz > max) max = sz;
            if (sz < min) min = sz;
            avg += sz;
        }
        avg /= patterns.size();
        stats.put("avgBranches", avg);
        stats.put("maxBranches", max);
        stats.put("minBranches", min);
        double avgE = 0.0, avgB = 0.0;
        int countE = 0, countB = 0;
        for (BranchPattern p : patterns) {
            Map<String, Double> div = calculateBranchDiversity(p);
            Double e = div.get("entropy");
            Double b = div.get("balance");
            if (e != null) { avgE += e; countE++; }
            if (b != null) { avgB += b; countB++; }
        }
        stats.put("avgEntropy", countE > 0 ? avgE / countE : 0.0);
        stats.put("avgBalance", countB > 0 ? avgB / countB : 0.0);
        return stats;
    }

    public static Pair<Integer, Double> findMostProbableBranch(BranchPattern pattern) {
        if (pattern.getBranchNodes().isEmpty()) return null;
        int maxIndex = -1;
        double maxV = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < pattern.getProbabilities().size(); i++) {
            double v = pattern.getProbabilities().get(i);
            if (v > maxV) { maxV = v; maxIndex = i; }
        }
        if (maxIndex < 0) return null;
        return new Pair<Integer, Double>(pattern.getBranchNodes().get(maxIndex), maxV);
    }

    public static Pair<Integer, Double> findLeastProbableBranch(BranchPattern pattern) {
        if (pattern.getBranchNodes().isEmpty()) return null;
        int minIndex = -1;
        double minV = Double.POSITIVE_INFINITY;
        for (int i = 0; i < pattern.getProbabilities().size(); i++) {
            double v = pattern.getProbabilities().get(i);
            if (v < minV) { minV = v; minIndex = i; }
        }
        if (minIndex < 0) return null;
        return new Pair<Integer, Double>(pattern.getBranchNodes().get(minIndex), minV);
    }
}
