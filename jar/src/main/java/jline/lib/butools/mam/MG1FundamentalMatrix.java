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

import jline.lib.smc.MG1FIOptions;
import jline.lib.smc.MG1_FI;
import jline.util.matrix.Matrix;

public final class MG1FundamentalMatrix {
    private MG1FundamentalMatrix() {}

    /**
     * Method for solving M/G/1 type matrix equation.
     */
    public enum MG1Method {
        CR,  // Cyclic Reduction
        RR,  // Ramaswami Reduction
        NI,  // Newton Iteration
        FI,  // Functional Iteration
        IS   // Invariant Subspace
    }

    /**
     * Returns matrix G corresponding to the M/G/1 type Markov chain defined by matrices A.
     *
     * Matrix G is the minimal non-negative solution of the following matrix equation:
     * G = A_0 + A_1*G + A_2*G^2 + A_3*G^3 + ...
     *
     * @param A List of matrix blocks of the M/G/1 type generator from 0 to M-1.
     * @param precision Matrix G is computed iteratively up to this precision
     * @param maxNumIt The maximal number of iterations
     * @param method The method used to solve the matrix-quadratic equation
     * @return The G matrix of the M/G/1 type Markov chain (G is stochastic)
     */
    public static Matrix mg1FundamentalMatrix(List<Matrix> A, double precision, int maxNumIt, MG1Method method) {
        if (A.isEmpty()) {
            throw new IllegalArgumentException("MG1FundamentalMatrix: At least one matrix block required");
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

        // Use SMC library solver (using Functional Iteration method)
        MG1FIOptions options = new MG1FIOptions("U-Based", maxNumIt, 0);
        return MG1_FI.mg1_fi(Am, options);
    }

    public static Matrix mg1FundamentalMatrix(List<Matrix> A) {
        return mg1FundamentalMatrix(A, 1e-14, 50, MG1Method.CR);
    }

    public static Matrix mg1FundamentalMatrix(List<Matrix> A, double precision) {
        return mg1FundamentalMatrix(A, precision, 50, MG1Method.CR);
    }

    public static Matrix mg1FundamentalMatrix(List<Matrix> A, double precision, int maxNumIt) {
        return mg1FundamentalMatrix(A, precision, maxNumIt, MG1Method.CR);
    }
}
