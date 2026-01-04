/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 *
 * Reference:
 * Bini, D. A., Meini, B., Steffé, S., Van Houdt, B. (2006, October).
 * Structured Markov chains solver: software tools. In Proceeding from the
 * 2006 workshop on Tools for solving structured Markov chains (p. 14). ACM.
 */
package jline.lib.butools.mam

import jline.lib.smc.gim1_R
import jline.lib.smc.GIM1ROptions
import jline.util.matrix.Matrix

/**
 * Method for solving G/M/1 type matrix equation.
 */
enum class GM1Method {
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
 * @param precision Matrix R is computed iteratively up to this precision (default 1e-14)
 * @param maxNumIt The maximal number of iterations (default 50)
 * @param method The method used to solve the matrix-quadratic equation (default CR)
 * @return The R matrix of the G/M/1 type Markov chain
 */
fun gm1FundamentalMatrix(
    A: List<Matrix>,
    precision: Double = 1e-14,
    maxNumIt: Int = 50,
    method: GM1Method = GM1Method.CR
): Matrix {
    if (A.isEmpty()) {
        throw IllegalArgumentException("GM1FundamentalMatrix: At least one matrix block required")
    }

    val N = A[0].numCols

    // Concatenate matrices horizontally
    val Am = Matrix.zeros(N, N * A.size)
    for (i in A.indices) {
        for (r in 0 until N) {
            for (c in 0 until N) {
                Am[r, i * N + c] = A[i][r, c]
            }
        }
    }

    // Use SMC library solver
    val options = GIM1ROptions(mode = "FI", maxNumIt = maxNumIt, verbose = 0)
    return gim1_R(Am, options)
}
