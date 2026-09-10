/**
 * @file Test if a CTMC has product form using Kolmogorov's criteria
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import java.util.ArrayList;
import java.util.List;

import jline.util.matrix.Matrix;

public final class Ctmc_testpf_kolmogorov {
    private Ctmc_testpf_kolmogorov() {}

    /**
     * Test if a CTMC has product form using Kolmogorov's criteria.
     *
     * Both the forward and the reverse product are taken on the same generator.
     * Taking the reverse product on the time-reversed generator instead, as this
     * method did before 2026-08-01, makes the test a tautology: the stationary
     * ratios of Qrev telescope to 1 around any cycle, so the two products were
     * identically equal and every chain passed.
     */
    public static boolean ctmc_testpf_kolmogorov(Matrix Q) {
        Matrix Q_norm = ctmc_makeStochastic(Q);

        int n = Q.length();

        Matrix A = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j && Q_norm.get(i, j) > 0) {
                    A.set(i, j, 1.0);
                }
            }
        }

        for (int start = 0; start < n; start++) {
            for (int target = 0; target < n; target++) {
                if (target != start && A.get(target, start) > 0) {
                    boolean[] unusedNodes = new boolean[n];
                    for (int i = 0; i < n; i++) unusedNodes[i] = true;
                    List<Integer> emptyPath = new ArrayList<Integer>();
                    List<List<Integer>> cycles = findPaths(A, unusedNodes, emptyPath, start, target);

                    for (List<Integer> cycle : cycles) {
                        List<Integer> completeCycle = new ArrayList<Integer>(cycle);
                        completeCycle.add(start);

                        double q = 1.0;
                        double qr = 1.0;
                        for (int i = 0; i < completeCycle.size() - 1; i++) {
                            q *= Q_norm.get(completeCycle.get(i), completeCycle.get(i + 1));
                            qr *= Q_norm.get(completeCycle.get(i + 1), completeCycle.get(i));
                        }

                        if (q == 0.0) {
                            continue;
                        }
                        if (Math.abs(q - qr) / Math.abs(q) > 1e-6) {
                            return false;
                        }
                    }
                }
            }
        }

        return true;
    }

    private static Matrix ctmc_makeStochastic(Matrix Q) {
        int n = Q.length();
        Matrix result = Q.copy();

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j && result.get(i, j) < 0) {
                    result.set(i, j, 0.0);
                }
            }
        }

        for (int i = 0; i < n; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    rowSum += result.get(i, j);
                }
            }
            result.set(i, i, -rowSum);
        }

        return result;
    }

    private static List<List<Integer>> findPaths(Matrix adj, boolean[] nodes,
                                                 List<Integer> currentPath, int start, int target) {
        List<List<Integer>> paths = new ArrayList<List<Integer>>();

        boolean[] newNodes = nodes.clone();
        newNodes[start] = false;

        List<Integer> newPath = new ArrayList<Integer>(currentPath);
        newPath.add(start);

        if (start == target) {
            paths.add(newPath);
            return paths;
        }

        List<Integer> childList = new ArrayList<Integer>();
        for (int j = 0; j < adj.getNumCols(); j++) {
            if (adj.get(start, j) > 0 && newNodes[j]) {
                childList.add(j);
            }
        }

        if (childList.isEmpty()) {
            return paths;
        }

        for (int child : childList) {
            List<List<Integer>> childPaths = findPaths(adj, newNodes, newPath, child, target);
            paths.addAll(childPaths);
        }

        return paths;
    }
}
