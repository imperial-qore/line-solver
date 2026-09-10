/**
 * @file Mean hitting times of a discrete-time Markov chain
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.List;

public final class Dtmc_hitting_time {
    private Dtmc_hitting_time() {}

    /**
     * Mean number of steps to reach any target state from each state of a DTMC.
     * Target states have zero hitting time; the others solve
     * (I - P_NT) h_NT = 1 over the non-target block. Twin of the MATLAB
     * dtmc_hitting_time and of the Python api.mc.dtmc_hitting_time; state
     * indices are 0-based, as in the rest of the JAR.
     *
     * @param P            transition matrix
     * @param targetStates indices of the target states
     * @return column vector of mean hitting times, one entry per state
     */
    public static Matrix dtmc_hitting_time(Matrix P, int[] targetStates) {
        int n = P.getNumRows();
        boolean[] isTarget = new boolean[n];
        for (int i = 0; i < targetStates.length; i++) {
            isTarget[targetStates[i]] = true;
        }
        List<Integer> nonTarget = new ArrayList<Integer>();
        for (int i = 0; i < n; i++) {
            if (!isTarget[i]) {
                nonTarget.add(i);
            }
        }

        Matrix h = new Matrix(n, 1);
        int m = nonTarget.size();
        if (m == 0) {
            return h;
        }

        Matrix A = new Matrix(m, m);
        Matrix b = new Matrix(m, 1);
        for (int i = 0; i < m; i++) {
            b.set(i, 0, 1.0);
            for (int j = 0; j < m; j++) {
                double p = P.get(nonTarget.get(i), nonTarget.get(j));
                A.set(i, j, (i == j ? 1.0 : 0.0) - p);
            }
        }
        Matrix hNT = new Matrix(m, 1);
        if (!Matrix.solveSafe(A, b, hNT)) {
            throw new RuntimeException("The hitting-time system is singular: a state cannot reach the target set.");
        }
        for (int i = 0; i < m; i++) {
            h.set(nonTarget.get(i), 0, hNT.get(i, 0));
        }
        return h;
    }
}
