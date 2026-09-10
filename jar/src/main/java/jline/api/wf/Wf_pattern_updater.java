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

                List<ServiceParameters> sequenceParams = collectParams(sequence, serviceParams);
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
                    updatedMatrix = removeMatrixRows(updatedMatrix, rowsToRemove);

                    replaceNodeReferences(updatedMatrix, forkNode, firstNode);
                    replaceNodeReferences(updatedMatrix, joinNode, firstNode);

                    List<ServiceParameters> parallelParams = collectParams(parallel, serviceParams);
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
                updatedMatrix = removeMatrixRows(updatedMatrix, rowsToRemove);
                for (Integer router : involvedRouters) {
                    replaceNodeReferences(updatedMatrix, router, loopNode);
                }
                updateLoopTransitionProbabilities(updatedMatrix, loopNode);

                ServiceParameters originalParams = serviceParams.get(loopNode);
                if (originalParams != null) {
                    serviceParams.put(loopNode, convolveLoop(originalParams, loopProb));
                }
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
                updatedMatrix = removeMatrixRows(updatedMatrix, rowsToRemove);

                replaceNodeReferences(updatedMatrix, branch.getForkNode(), firstNode);
                if (branch.getJoinNode() != null) {
                    replaceNodeReferences(updatedMatrix, branch.getJoinNode(), firstNode);
                }

                List<ServiceParameters> branchParams = collectParams(branchNodes, serviceParams);
                ServiceParameters convolvedParams = convolveBranches(branchParams, branch.getProbabilities());
                serviceParams.put(firstNode, convolvedParams);

                for (int i = 1; i < branchNodes.size(); i++) {
                    serviceParams.remove(branchNodes.get(i));
                }
            }
        }

        return updatedMatrix;
    }

    /** The laws of the listed nodes, skipping the ones with none, as in the reference. */
    private static List<ServiceParameters> collectParams(
            List<Integer> nodes, Map<Integer, ServiceParameters> serviceParams) {
        List<ServiceParameters> out = new ArrayList<ServiceParameters>();
        for (Integer n : nodes) {
            ServiceParameters p = serviceParams.get(n);
            if (p != null) out.add(p);
        }
        return out;
    }

    static Matrix removeMatrixRows(Matrix matrix, List<Integer> rowsToRemove) {
        if (rowsToRemove.isEmpty() || matrix.getNumRows() == 0) return matrix;
        boolean[] keep = new boolean[matrix.getNumRows()];
        java.util.Arrays.fill(keep, true);
        for (Integer row : rowsToRemove) {
            if (row != null && row >= 0 && row < matrix.getNumRows()) keep[row] = false;
        }
        int newRows = 0;
        for (int i = 0; i < keep.length; i++) {
            if (keep[i]) newRows++;
        }
        Matrix result = Matrix.zeros(newRows, matrix.getNumCols());
        int newRowIndex = 0;
        for (int i = 0; i < matrix.getNumRows(); i++) {
            if (!keep[i]) continue;
            for (int j = 0; j < matrix.getNumCols(); j++) {
                result.set(newRowIndex, j, matrix.get(i, j));
            }
            newRowIndex++;
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

    /** Flatten an entry law held as a row or a column matrix. */
    private static double[] alphaVector(Matrix alpha) {
        double[] v = new double[alpha.getNumRows() * alpha.getNumCols()];
        int k = 0;
        for (int i = 0; i < alpha.getNumRows(); i++) {
            for (int j = 0; j < alpha.getNumCols(); j++) {
                v[k++] = alpha.get(i, j);
            }
        }
        return v;
    }

    private static Matrix rowVector(double[] v) {
        Matrix m = Matrix.zeros(1, v.length);
        for (int i = 0; i < v.length; i++) {
            m.set(0, i, v[i]);
        }
        return m;
    }

    /** -T e, the exit rate out of each phase. */
    private static double[] exitRate(Matrix T) {
        double[] r = new double[T.getNumRows()];
        for (int i = 0; i < T.getNumRows(); i++) {
            double s = 0.0;
            for (int j = 0; j < T.getNumCols(); j++) {
                s += T.get(i, j);
            }
            r[i] = -s;
        }
        return r;
    }

    /** 1 - alpha e, the mass the entry law leaves for instantaneous completion. */
    private static double exitProb(double[] alpha) {
        double s = 0.0;
        for (int i = 0; i < alpha.length; i++) {
            s += alpha[i];
        }
        return 1.0 - s;
    }

    /** Copy src into dst with its top-left corner at (r0, c0). */
    private static void place(Matrix dst, Matrix src, int r0, int c0) {
        for (int i = 0; i < src.getNumRows(); i++) {
            for (int j = 0; j < src.getNumCols(); j++) {
                dst.set(r0 + i, c0 + j, src.get(i, j));
            }
        }
    }

    private static ServiceParameters unitParams() {
        return new ServiceParameters(Matrix.ones(1, 1), Matrix.zeros(1, 1));
    }

    /** Convolution of the durations, i.e. the service laws run one after another. */
    static ServiceParameters convolveSequence(List<ServiceParameters> params) {
        if (params.isEmpty()) return unitParams();
        if (params.size() == 1) return params.get(0);

        double[] alpha = alphaVector(params.get(0).getAlpha());
        Matrix T = params.get(0).getT();
        for (int k = 1; k < params.size(); k++) {
            double[] a2 = alphaVector(params.get(k).getAlpha());
            Matrix T2 = params.get(k).getT();
            int n1 = alpha.length;
            int n2 = a2.length;
            double ex = exitProb(alpha);
            double[] er = exitRate(T);

            double[] na = new double[n1 + n2];
            System.arraycopy(alpha, 0, na, 0, n1);
            for (int j = 0; j < n2; j++) {
                na[n1 + j] = ex * a2[j];
            }

            Matrix nt = Matrix.zeros(n1 + n2, n1 + n2);
            place(nt, T, 0, 0);
            place(nt, T2, n1, n1);
            for (int i = 0; i < n1; i++) {
                for (int j = 0; j < n2; j++) {
                    nt.set(i, n1 + j, er[i] * a2[j]);
                }
            }

            alpha = na;
            T = nt;
        }
        return new ServiceParameters(rowVector(alpha), T);
    }

    /** Maximum of the durations, i.e. a fork whose join waits for every branch. */
    static ServiceParameters convolveParallel(List<ServiceParameters> params) {
        if (params.isEmpty()) return unitParams();
        if (params.size() == 1) return params.get(0);

        double[] alpha = alphaVector(params.get(0).getAlpha());
        Matrix T = params.get(0).getT();
        for (int k = 1; k < params.size(); k++) {
            double[] a2 = alphaVector(params.get(k).getAlpha());
            Matrix T2 = params.get(k).getT();
            int n1 = alpha.length;
            int n2 = a2.length;
            int np = n1 * n2;
            double e1 = exitProb(alpha);
            double e2 = exitProb(a2);
            double[] r1 = exitRate(T);
            double[] r2 = exitRate(T2);

            double[] na = new double[np + n1 + n2];
            for (int i = 0; i < n1; i++) {
                for (int j = 0; j < n2; j++) {
                    na[i * n2 + j] = alpha[i] * a2[j];
                }
            }
            for (int i = 0; i < n1; i++) {
                na[np + i] = e2 * alpha[i];
            }
            for (int j = 0; j < n2; j++) {
                na[np + n1 + j] = e1 * a2[j];
            }

            Matrix nt = Matrix.zeros(np + n1 + n2, np + n1 + n2);
            for (int i = 0; i < n1; i++) {
                for (int j = 0; j < n2; j++) {
                    int r = i * n2 + j;
                    // T1 (x) I + I (x) T2 on the block where both branches are alive
                    for (int ii = 0; ii < n1; ii++) {
                        nt.set(r, ii * n2 + j, nt.get(r, ii * n2 + j) + T.get(i, ii));
                    }
                    for (int jj = 0; jj < n2; jj++) {
                        nt.set(r, i * n2 + jj, nt.get(r, i * n2 + jj) + T2.get(j, jj));
                    }
                    // branch 2 finishes first -> only branch 1 is left, in phase i
                    nt.set(r, np + i, r2[j]);
                    // branch 1 finishes first -> only branch 2 is left, in phase j
                    nt.set(r, np + n1 + j, r1[i]);
                }
            }
            place(nt, T, np, np);
            place(nt, T2, np + n1, np + n1);

            alpha = na;
            T = nt;
        }
        return new ServiceParameters(rowVector(alpha), T);
    }

    /** Geometric repetition: the exit flow re-enters through alpha with probability loopProb. */
    static ServiceParameters convolveLoop(ServiceParameters params, double loopProb) {
        if (loopProb <= 0.0 || loopProb >= 1.0) return params;
        double[] alpha = alphaVector(params.getAlpha());
        double[] er = exitRate(params.getT());
        Matrix nt = params.getT().copy();
        for (int i = 0; i < nt.getNumRows(); i++) {
            for (int j = 0; j < nt.getNumCols(); j++) {
                nt.set(i, j, nt.get(i, j) + loopProb * er[i] * alpha[j]);
            }
        }
        return new ServiceParameters(rowVector(alpha), nt);
    }

    /** Probabilistic choice among the alternatives, on a block-diagonal generator. */
    static ServiceParameters convolveBranches(List<ServiceParameters> params, List<Double> probs) {
        if (params.isEmpty()) return unitParams();
        if (params.size() == 1) return params.get(0);

        double[] w = new double[params.size()];
        double total = 0.0;
        for (int i = 0; i < params.size(); i++) {
            Double pi = (probs != null && i < probs.size()) ? probs.get(i) : null;
            w[i] = (pi == null) ? 0.0 : pi.doubleValue();
            total += w[i];
        }
        if (total > 0.0) {
            for (int i = 0; i < w.length; i++) {
                w[i] /= total;
            }
        } else {
            for (int i = 0; i < w.length; i++) {
                w[i] = 1.0 / params.size();
            }
        }

        int n = 0;
        for (int i = 0; i < params.size(); i++) {
            n += params.get(i).getAlpha().getNumRows() * params.get(i).getAlpha().getNumCols();
        }
        double[] alpha = new double[n];
        Matrix T = Matrix.zeros(n, n);
        int off = 0;
        for (int i = 0; i < params.size(); i++) {
            double[] ai = alphaVector(params.get(i).getAlpha());
            for (int j = 0; j < ai.length; j++) {
                alpha[off + j] = w[i] * ai[j];
            }
            place(T, params.get(i).getT(), off, off);
            off += ai.length;
        }
        return new ServiceParameters(rowVector(alpha), T);
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
