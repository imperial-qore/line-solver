/**
 * @file Workflow Pattern Updater
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

import jline.GlobalConstants;
import jline.util.matrix.Matrix;

/**
 * Workflow pattern updater for simplifying and optimizing workflow networks.
 */
public final class Wf_pattern_updater {

    private Wf_pattern_updater() {}

    /** Service parameters (phase-type distributions). */
    public static final class ServiceParameters {
        private final Matrix alpha;
        private final Matrix T;

        public ServiceParameters(Matrix alpha, Matrix T) {
            this.alpha = alpha;
            this.T = T;
        }
        public Matrix getAlpha() { return alpha; }
        public Matrix getT() { return T; }
    }

    /** Updated workflow representation. */
    public static final class UpdatedWorkflow {
        private final Matrix linkMatrix;
        private final Map<Integer, ServiceParameters> serviceParameters;

        public UpdatedWorkflow(Matrix linkMatrix, Map<Integer, ServiceParameters> serviceParameters) {
            this.linkMatrix = linkMatrix;
            this.serviceParameters = serviceParameters;
        }
        public Matrix getLinkMatrix() { return linkMatrix; }
        public Map<Integer, ServiceParameters> getServiceParameters() { return serviceParameters; }
    }

    public static UpdatedWorkflow updatePatterns(
            Matrix linkMatrix,
            List<Integer> serviceNodes,
            List<Integer> forkNodes,
            List<Integer> joinNodes,
            List<Integer> routerNodes,
            Map<Integer, ServiceParameters> serviceParams) {
        Matrix currentMatrix = linkMatrix.copy();
        Map<Integer, ServiceParameters> currentParams = new HashMap<Integer, ServiceParameters>(serviceParams);

        List<List<Integer>> sequences = Wf_sequence_detector.detectSequences(currentMatrix, serviceNodes);
        currentMatrix = updateSequencePatterns(currentMatrix, sequences, serviceNodes, currentParams);

        List<List<Integer>> parallels = Wf_parallel_detector.detectParallel(currentMatrix, serviceNodes, forkNodes, joinNodes);
        currentMatrix = updateParallelPatterns(currentMatrix, parallels, serviceNodes, currentParams);

        List<Integer> loops = Wf_loop_detector.detectLoops(currentMatrix, serviceNodes, routerNodes, joinNodes);
        currentMatrix = updateLoopPatterns(currentMatrix, loops, serviceNodes, routerNodes, currentParams);

        List<Wf_branch_detector.BranchPattern> branches = Wf_branch_detector.detectBranches(currentMatrix, serviceNodes, joinNodes);
        currentMatrix = updateBranchPatterns(currentMatrix, branches, serviceNodes, currentParams);

        return new UpdatedWorkflow(currentMatrix, currentParams);
    }

    private static Matrix updateSequencePatterns(
            Matrix linkMatrix,
            List<List<Integer>> sequences,
            List<Integer> serviceNodes,
            Map<Integer, ServiceParameters> serviceParams) {
        if (sequences.isEmpty()) return linkMatrix;
        Set<Integer> serviceSet = new HashSet<Integer>(serviceNodes);
        Matrix updatedMatrix = linkMatrix.copy();
        List<Integer> rowsToRemove = new ArrayList<Integer>();

        for (int i = 0; i < updatedMatrix.getNumRows(); i++) {
            int start = (int) updatedMatrix.get(i, 0);
            int end = (int) updatedMatrix.get(i, 1);
            if (serviceSet.contains(start) && serviceSet.contains(end)) {
                for (List<Integer> sequence : sequences) {
                    for (int j = 0; j < sequence.size() - 1; j++) {
                        if (sequence.get(j) == start && sequence.get(j + 1) == end) {
                            rowsToRemove.add(i);
                            break;
                        }
                    }
                }
            }
        }

        Matrix finalMatrix = removeMatrixRows(updatedMatrix, rowsToRemove);

        for (List<Integer> sequence : sequences) {
            if (sequence.size() >= 2) {
                int firstNode = sequence.get(0);
                int lastNode = sequence.get(sequence.size() - 1);
                replaceNodeReferences(finalMatrix, lastNode, firstNode);

                List<ServiceParameters> sequenceParams = new ArrayList<ServiceParameters>();
                for (Integer n : sequence) sequenceParams.add(serviceParams.get(n));
                ServiceParameters convolvedParams = convolveSequence(sequenceParams);
                serviceParams.put(firstNode, convolvedParams);

                for (int i = 1; i < sequence.size(); i++) {
                    serviceParams.remove(sequence.get(i));
                }
            }
        }

        return finalMatrix;
    }

    private static Matrix updateParallelPatterns(
            Matrix linkMatrix,
            List<List<Integer>> parallels,
            List<Integer> serviceNodes,
            Map<Integer, ServiceParameters> serviceParams) {
        if (parallels.isEmpty()) return linkMatrix;
        Matrix updatedMatrix = linkMatrix.copy();

        for (List<Integer> parallel : parallels) {
            if (parallel.size() >= 2) {
                int firstNode = parallel.get(0);
                int[] forkJoin = findForkJoinForParallel(updatedMatrix, parallel);
                Integer forkNode = forkJoin[0] == -1 ? null : forkJoin[0];
                Integer joinNode = forkJoin[1] == -1 ? null : forkJoin[1];

                if (forkNode != null && joinNode != null) {
                    List<Integer> rowsToRemove = findRowsInvolvingNodes(updatedMatrix, parallel);
                    removeMatrixRows(updatedMatrix, rowsToRemove);

                    replaceNodeReferences(updatedMatrix, forkNode, firstNode);
                    replaceNodeReferences(updatedMatrix, joinNode, firstNode);

                    List<ServiceParameters> parallelParams = new ArrayList<ServiceParameters>();
                    for (Integer n : parallel) parallelParams.add(serviceParams.get(n));
                    ServiceParameters convolvedParams = convolveParallel(parallelParams);
                    serviceParams.put(firstNode, convolvedParams);

                    for (int i = 1; i < parallel.size(); i++) {
                        serviceParams.remove(parallel.get(i));
                    }
                }
            }
        }

        return updatedMatrix;
    }

    private static Matrix updateLoopPatterns(
            Matrix linkMatrix,
            List<Integer> loops,
            List<Integer> serviceNodes,
            List<Integer> routerNodes,
            Map<Integer, ServiceParameters> serviceParams) {
        if (loops.isEmpty()) return linkMatrix;
        Matrix updatedMatrix = linkMatrix.copy();
        Set<Integer> routerSet = new HashSet<Integer>(routerNodes);

        for (Integer loopNode : loops) {
            double loopProb = Wf_loop_detector.getLoopProbability(loopNode, updatedMatrix, routerNodes);
            if (loopProb > GlobalConstants.Zero) {
                List<Integer> involvedRouters = findRoutersForLoop(updatedMatrix, loopNode, routerSet);
                List<Integer> rowsToRemove = findLoopConnections(updatedMatrix, loopNode, involvedRouters);
                removeMatrixRows(updatedMatrix, rowsToRemove);
                for (Integer router : involvedRouters) {
                    replaceNodeReferences(updatedMatrix, router, loopNode);
                }
                updateLoopTransitionProbabilities(updatedMatrix, loopNode);

                ServiceParameters originalParams = serviceParams.get(loopNode);
                ServiceParameters loopedParams = convolveLoop(originalParams, loopProb);
                serviceParams.put(loopNode, loopedParams);
            }
        }
        return updatedMatrix;
    }

    private static Matrix updateBranchPatterns(
            Matrix linkMatrix,
            List<Wf_branch_detector.BranchPattern> branches,
            List<Integer> serviceNodes,
            Map<Integer, ServiceParameters> serviceParams) {
        if (branches.isEmpty()) return linkMatrix;
        Matrix updatedMatrix = linkMatrix.copy();

        for (Wf_branch_detector.BranchPattern branch : branches) {
            List<Integer> branchNodes = branch.getBranchNodes();
            if (branchNodes.size() >= 2 && branch.getForkNode() != null) {
                int firstNode = branchNodes.get(0);
                List<Integer> rowsToRemove = findRowsInvolvingNodes(updatedMatrix, branchNodes);
                removeMatrixRows(updatedMatrix, rowsToRemove);

                replaceNodeReferences(updatedMatrix, branch.getForkNode(), firstNode);
                if (branch.getJoinNode() != null) {
                    replaceNodeReferences(updatedMatrix, branch.getJoinNode(), firstNode);
                }

                List<ServiceParameters> branchParams = new ArrayList<ServiceParameters>();
                for (Integer n : branchNodes) branchParams.add(serviceParams.get(n));
                ServiceParameters convolvedParams = convolveBranches(branchParams, branch.getProbabilities());
                serviceParams.put(firstNode, convolvedParams);

                for (int i = 1; i < branchNodes.size(); i++) {
                    serviceParams.remove(branchNodes.get(i));
                }
            }
        }

        return updatedMatrix;
    }

    private static Matrix removeMatrixRows(Matrix matrix, List<Integer> rowsToRemove) {
        if (rowsToRemove.isEmpty()) return matrix;
        List<Integer> sorted = new ArrayList<Integer>(rowsToRemove);
        java.util.Collections.sort(sorted, java.util.Collections.<Integer>reverseOrder());
        Matrix result = matrix.copy();

        for (Integer row : sorted) {
            if (row >= 0 && row < result.getNumRows()) {
                int newRows = result.getNumRows() - 1;
                Matrix newMatrix = Matrix.zeros(newRows, result.getNumCols());
                int newRowIndex = 0;
                for (int i = 0; i < result.getNumRows(); i++) {
                    if (i != row) {
                        for (int j = 0; j < result.getNumCols(); j++) {
                            newMatrix.set(newRowIndex, j, result.get(i, j));
                        }
                        newRowIndex++;
                    }
                }
                return newMatrix;
            }
        }
        return result;
    }

    private static void replaceNodeReferences(Matrix matrix, int oldNode, int newNode) {
        for (int i = 0; i < matrix.getNumRows(); i++) {
            if ((int) matrix.get(i, 0) == oldNode) {
                matrix.set(i, 0, (double) newNode);
            }
            if ((int) matrix.get(i, 1) == oldNode) {
                matrix.set(i, 1, (double) newNode);
            }
        }
    }

    private static List<Integer> findRowsInvolvingNodes(Matrix matrix, List<Integer> nodes) {
        Set<Integer> nodeSet = new HashSet<Integer>(nodes);
        List<Integer> rowsToRemove = new ArrayList<Integer>();
        for (int i = 0; i < matrix.getNumRows(); i++) {
            int start = (int) matrix.get(i, 0);
            int end = (int) matrix.get(i, 1);
            if (nodeSet.contains(start) || nodeSet.contains(end)) {
                rowsToRemove.add(i);
            }
        }
        return rowsToRemove;
    }

    private static int[] findForkJoinForParallel(Matrix matrix, List<Integer> parallel) {
        return new int[]{-1, -1};
    }

    private static List<Integer> findRoutersForLoop(Matrix matrix, int loopNode, Set<Integer> routerSet) {
        List<Integer> routers = new ArrayList<Integer>();
        for (int i = 0; i < matrix.getNumRows(); i++) {
            int start = (int) matrix.get(i, 0);
            int end = (int) matrix.get(i, 1);
            if ((start == loopNode && routerSet.contains(end))
                    || (end == loopNode && routerSet.contains(start))) {
                if (routerSet.contains(start)) routers.add(start);
                if (routerSet.contains(end)) routers.add(end);
            }
        }
        Set<Integer> dedup = new java.util.LinkedHashSet<Integer>(routers);
        return new ArrayList<Integer>(dedup);
    }

    private static List<Integer> findLoopConnections(Matrix matrix, int loopNode, List<Integer> routers) {
        List<Integer> connections = new ArrayList<Integer>();
        Set<Integer> routerSet = new HashSet<Integer>(routers);
        for (int i = 0; i < matrix.getNumRows(); i++) {
            int start = (int) matrix.get(i, 0);
            int end = (int) matrix.get(i, 1);
            if ((start == loopNode && routerSet.contains(end))
                    || (end == loopNode && routerSet.contains(start))
                    || (routerSet.contains(start) && end == loopNode)) {
                connections.add(i);
            }
        }
        return connections;
    }

    private static void updateLoopTransitionProbabilities(Matrix matrix, int loopNode) {
        for (int i = 0; i < matrix.getNumRows(); i++) {
            int start = (int) matrix.get(i, 0);
            if (start == loopNode) {
                matrix.set(i, 2, 1.0);
            }
        }
    }

    private static ServiceParameters convolveSequence(List<ServiceParameters> params) {
        if (!params.isEmpty()) return params.get(0);
        return new ServiceParameters(Matrix.ones(1, 1), Matrix.zeros(1, 1));
    }

    private static ServiceParameters convolveParallel(List<ServiceParameters> params) {
        if (!params.isEmpty()) return params.get(0);
        return new ServiceParameters(Matrix.ones(1, 1), Matrix.zeros(1, 1));
    }

    private static ServiceParameters convolveLoop(ServiceParameters params, double loopProb) {
        return params;
    }

    private static ServiceParameters convolveBranches(List<ServiceParameters> params, List<Double> probs) {
        if (!params.isEmpty()) return params.get(0);
        return new ServiceParameters(Matrix.ones(1, 1), Matrix.zeros(1, 1));
    }

    public static boolean validateUpdatedWorkflow(UpdatedWorkflow workflow) {
        // see _kb/03-api-layer.md for rationale (wf/ section)
        Set<Integer> referencedNodes = new HashSet<Integer>();
        for (int i = 0; i < workflow.linkMatrix.getNumRows(); i++) {
            referencedNodes.add((int) workflow.linkMatrix.get(i, 0));
            referencedNodes.add((int) workflow.linkMatrix.get(i, 1));
        }
        for (Integer node : workflow.serviceParameters.keySet()) {
            if (!referencedNodes.contains(node)) {
                return false;
            }
        }
        return true;
    }

    public static Map<String, Object> getUpdateStats(Matrix originalMatrix, UpdatedWorkflow updatedWorkflow) {
        Map<String, Object> stats = new HashMap<String, Object>();
        stats.put("originalLinks", originalMatrix.getNumRows());
        stats.put("updatedLinks", updatedWorkflow.linkMatrix.getNumRows());
        stats.put("linksReduced", originalMatrix.getNumRows() - updatedWorkflow.linkMatrix.getNumRows());
        stats.put("serviceNodes", updatedWorkflow.serviceParameters.size());

        double reductionRatio = 0.0;
        if (originalMatrix.getNumRows() > 0) {
            reductionRatio = (double) (originalMatrix.getNumRows() - updatedWorkflow.linkMatrix.getNumRows())
                    / originalMatrix.getNumRows();
        }
        stats.put("reductionRatio", reductionRatio);
        return stats;
    }
}
