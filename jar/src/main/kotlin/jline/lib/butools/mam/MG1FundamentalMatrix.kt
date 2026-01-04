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

import jline.lib.smc.mg1_fi
import jline.lib.smc.MG1FIOptions
import jline.util.matrix.Matrix

/**
 * Method for solving M/G/1 type matrix equation.
 */
enum class MG1Method {
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
 * @param precision Matrix G is computed iteratively up to this precision (default 1e-14)
 * @param maxNumIt The maximal number of iterations (default 50)
 * @param method The method used to solve the matrix-quadratic equation (default CR)
 * @return The G matrix of the M/G/1 type Markov chain (G is stochastic)
 */
fun mg1FundamentalMatrix(
    A: List<Matrix>,
    precision: Double = 1e-14,
    maxNumIt: Int = 50,
    method: MG1Method = MG1Method.CR
): Matrix {
    if (A.isEmpty()) {
        throw IllegalArgumentException("MG1FundamentalMatrix: At least one matrix block required")
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

    // Use SMC library solver (using Functional Iteration method)
    val options = MG1FIOptions(mode = "U-Based", maxNumIt = maxNumIt, verbose = 0)
    return mg1_fi(Am, options)
}
