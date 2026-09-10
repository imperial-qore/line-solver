/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 *
 * Reference:
 * Bini, D. A., Meini, B., Steffe, S., Van Houdt, B. (2006, October).
 * Structured Markov chains solver: software tools. In Proceeding from the
 * 2006 workshop on Tools for solving structured Markov chains (p. 14). ACM.
 */
package jline.lib.butools.mam;

import java.util.List;

import jline.lib.smc.GIM1_R.GIM1ROptions;
import jline.lib.smc.GIM1_R;
import jline.util.matrix.Matrix;

public final class GM1FundamentalMatrix {
    private GM1FundamentalMatrix() {}

    /**
     * Method for solving G/M/1 type matrix equation.
     */
    public enum GM1Method {
        CR,  // Cyclic Reduction
        RR,  // Ramaswami Reduction
        NI,  // Newton Iteration
        FI,  // Functional Iteration
        IS   // Invariant Subspace
    }

    /**
     * Returns matrix R corresponding to the G/M/1 type Markov chain given by matrices A.
     *
     * Matrix R is the minimal non-negative solution of the following matrix equation:
     * R = A_0 + R*A_1 + R^2*A_2 + R^3*A_3 + ...
     *
     * @param A List of matrix blocks of the G/M/1 type generator from 0 to M-1.
     * @param precision Matrix R is computed iteratively up to this precision
     * @param maxNumIt The maximal number of iterations
     * @param method The method used to solve the matrix-quadratic equation
     * @return The R matrix of the G/M/1 type Markov chain
     */
    public static Matrix gm1FundamentalMatrix(List<Matrix> A, double precision, int maxNumIt, GM1Method method) {
        if (A.isEmpty()) {
            throw new IllegalArgumentException("GM1FundamentalMatrix: At least one matrix block required");
        }

        int N = A.get(0).getNumCols();

        // Concatenate matrices horizontally
        Matrix Am = Matrix.zeros(N, N * A.size());
        for (int i = 0; i < A.size(); i++) {
            for (int r = 0; r < N; r++) {
                for (int c = 0; c < N; c++) {
                    Am.set(r, i * N + c, A.get(i).get(r, c));
                }
            }
        }

        // Use SMC library solver
        GIM1ROptions options = new GIM1ROptions("FI", maxNumIt, 0);
        return GIM1_R.gim1_R(Am, options);
    }

    public static Matrix gm1FundamentalMatrix(List<Matrix> A) {
        return gm1FundamentalMatrix(A, 1e-14, 50, GM1Method.CR);
    }

    public static Matrix gm1FundamentalMatrix(List<Matrix> A, double precision) {
        return gm1FundamentalMatrix(A, precision, 50, GM1Method.CR);
    }

    public static Matrix gm1FundamentalMatrix(List<Matrix> A, double precision, int maxNumIt) {
        return gm1FundamentalMatrix(A, precision, maxNumIt, GM1Method.CR);
    }
}
