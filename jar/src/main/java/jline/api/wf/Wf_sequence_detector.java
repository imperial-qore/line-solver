/**
 * @file Workflow Sequence Pattern Detection
 *
 * @since LINE 3.0
 */
package jline.api.wf;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import jline.util.Pair;
import jline.util.matrix.Matrix;

/**
 * Automatic sequence structure detector for workflow networks.
 */
public final class Wf_sequence_detector {
    private Wf_sequence_detector() {}

    /**
     * Detect sequence patterns in a workflow network.
     */
    public static List<List<Integer>> detectSequences(Matrix linkMatrix, List<Integer> serviceNodes) {
        List<List<Integer>> chains = new ArrayList<List<Integer>>();

        List<Pair<Integer, Integer>> serviceConnections = findServiceConnections(linkMatrix, serviceNodes);
        if (serviceConnections.isEmpty()) {
            return chains;
        }

        Map<Integer, Integer> nodeCounts = countNodeOccurrences(serviceConnections);

        int countOnce = 0;
        for (int v : nodeCounts.values()) {
            if (v == 1) countOnce++;
        }
        int numSequences = countOnce / 2;

        List<Pair<Integer, Integer>> remainingConnections = new ArrayList<Pair<Integer, Integer>>(serviceConnections);

        for (int seqIndex = 0; seqIndex < numSequences; seqIndex++) {
            if (remainingConnections.isEmpty()) break;

            List<Integer> sequence = buildSequenceChain(remainingConnections);
            if (!sequence.isEmpty()) {
                chains.add(sequence);
            }
        }

        return chains;
    }

    /**
     * Find all connections between service nodes in the workflow.
     */
    private static List<Pair<Integer, Integer>> findServiceConnections(Matrix linkMatrix, List<Integer> serviceNodes) {
        List<Pair<Integer, Integer>> connections = new ArrayList<Pair<Integer, Integer>>();
        Set<Integer> serviceSet = new HashSet<Integer>(serviceNodes);

        for (int i = 0; i < linkMatrix.getNumRows(); i++) {
            int startNode = (int) linkMatrix.get(i, 0);
            int endNode = (int) linkMatrix.get(i, 1);

            if (serviceSet.contains(startNode) && serviceSet.contains(endNode)) {
                connections.add(new Pair<Integer, Integer>(startNode, endNode));
            }
        }

        return connections;
    }

    /**
     * Count how many times each service node appears in connections.
     */
    private static Map<Integer, Integer> countNodeOccurrences(List<Pair<Integer, Integer>> connections) {
        Map<Integer, Integer> counts = new HashMap<Integer, Integer>();

        for (Pair<Integer, Integer> p : connections) {
            int start = p.getLeft();
            int end = p.getRight();
            counts.put(start, counts.containsKey(start) ? counts.get(start) + 1 : 1);
            counts.put(end, counts.containsKey(end) ? counts.get(end) + 1 : 1);
        }

        return counts;
    }

    /**
     * Build a sequence chain starting from the first available connection.
     */
    private static List<Integer> buildSequenceChain(List<Pair<Integer, Integer>> connections) {
        if (connections.isEmpty()) return new ArrayList<Integer>();

        List<Integer> sequence = new ArrayList<Integer>();
        List<Integer> usedConnections = new ArrayList<Integer>();

        Pair<Integer, Integer> firstConn = connections.get(0);
        int first = firstConn.getLeft();
        int last = firstConn.getRight();
        sequence.add(first);
        sequence.add(last);
        usedConnections.add(0);

        boolean foundExtension = true;
        while (foundExtension) {
            foundExtension = false;
            int currentSize = sequence.size();

            for (int i = 1; i < connections.size(); i++) {
                if (usedConnections.contains(i)) continue;

                Pair<Integer, Integer> p = connections.get(i);
                int start = p.getLeft();
                int end = p.getRight();

                if (start == last) {
                    last = end;
                    sequence.add(end);
                    usedConnections.add(i);
                    foundExtension = true;
                } else if (end == first) {
                    first = start;
                    sequence.add(0, start);
                    usedConnections.add(i);
                    foundExtension = true;
                }
            }

            foundExtension = foundExtension && sequence.size() > currentSize;
        }

        // Remove used connections (from highest to lowest index)
        java.util.Collections.sort(usedConnections, java.util.Collections.<Integer>reverseOrder());
        for (int index : usedConnections) {
            connections.remove(index);
        }

        return sequence;
    }

    /**
     * Validate that a sequence chain is properly connected.
     */
    public static boolean validateSequence(List<Integer> sequence, Matrix linkMatrix) {
        if (sequence.size() < 2) return false;

        Set<Pair<Integer, Integer>> connections = new HashSet<Pair<Integer, Integer>>();
        for (int i = 0; i < linkMatrix.getNumRows(); i++) {
            int start = (int) linkMatrix.get(i, 0);
            int end = (int) linkMatrix.get(i, 1);
            connections.add(new Pair<Integer, Integer>(start, end));
        }

        for (int i = 0; i < sequence.size() - 1; i++) {
            Pair<Integer, Integer> connection = new Pair<Integer, Integer>(sequence.get(i), sequence.get(i + 1));
            if (!connections.contains(connection)) {
                return false;
            }
        }

        return true;
    }

    /**
     * Get sequence statistics for analysis.
     */
    public static Map<String, Object> getSequenceStats(List<List<Integer>> sequences) {
        Map<String, Object> stats = new HashMap<String, Object>();

        stats.put("numSequences", sequences.size());

        int totalNodes = 0;
        for (List<Integer> seq : sequences) {
            totalNodes += seq.size();
        }
        stats.put("totalNodes", totalNodes);

        double avgLength = 0.0;
        if (!sequences.isEmpty()) {
            double sum = 0.0;
            for (List<Integer> seq : sequences) {
                sum += seq.size();
            }
            avgLength = sum / sequences.size();
        }
        stats.put("avgLength", avgLength);

        int maxLength = 0;
        int minLength = 0;
        if (!sequences.isEmpty()) {
            maxLength = Integer.MIN_VALUE;
            minLength = Integer.MAX_VALUE;
            for (List<Integer> seq : sequences) {
                if (seq.size() > maxLength) maxLength = seq.size();
                if (seq.size() < minLength) minLength = seq.size();
            }
        }
        stats.put("maxLength", maxLength);
        stats.put("minLength", minLength);

        return stats;
    }
}
