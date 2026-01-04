/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import java.util.Random
import kotlin.math.floor
import kotlin.math.ln

/**
 * Generates random samples from a discrete marked Markovian arrival process.
 *
 * @param D The D0...DN matrices of the DMMAP (as MatrixCell)
 * @param K The number of samples to generate
 * @param initial Optional initial state (1-indexed). If not provided, drawn from stationary distribution.
 * @param prec Numerical precision for validation
 * @param random Random number generator
 * @return For single type (2 matrices): IntArray of inter-arrival times
 *         For multiple types: Array of IntArray where each element is [time, type]
 */
fun samplesFromDMMAP(
    D: MatrixCell,
    K: Int,
    initial: Int? = null,
    prec: Double = 1e-14,
    random: Random = Random()
): Any {
    if (!checkDMMAPRepresentation(D, prec)) {
        throw IllegalArgumentException("SamplesFromDMMAP: Input isn't a valid DMMAP representation!")
    }

    val N = D[0].numRows
    val numTypes = D.size()

    // Determine initial state
    var state = if (initial != null) {
        initial - 1 // Convert to 0-indexed
    } else {
        // Draw initial state according to stationary distribution
        val margDist = marginalDistributionFromDMMAP(D, prec)
        val stst = margDist.alpha

        // Compute cumulative distribution
        val cummInitial = DoubleArray(N)
        cummInitial[0] = stst[0, 0]
        for (i in 1 until N) {
            cummInitial[i] = cummInitial[i - 1] + stst[0, i]
        }

        // Draw initial state
        val r = random.nextDouble()
        var s = 0
        while (s < N - 1 && cummInitial[s] <= r) {
            s++
        }
        s
    }

    // Auxiliary variables
    val diag0 = DoubleArray(N) { i -> D[0][i, i] }
    val sojourn = DoubleArray(N) { i -> 1.0 / (1.0 - diag0[i]) }
    val logp = DoubleArray(N) { i -> ln(diag0[i]) }

    // Build next state probability matrix
    // nextpr = cumsum([diag(sojourn)*(D0 - diag), diag(sojourn)*D1, ...], 2)
    val totalCols = N * numTypes
    val nextpr = Matrix(N, totalCols)

    for (i in 0 until N) {
        var colOffset = 0
        for (dIdx in 0 until numTypes) {
            for (j in 0 until N) {
                val value = if (dIdx == 0 && i == j) {
                    0.0 // Remove diagonal from D0
                } else {
                    sojourn[i] * D[dIdx][i, j]
                }
                nextpr[i, colOffset + j] = value
            }
            colOffset += N
        }
    }

    // Compute cumulative sum along rows
    for (i in 0 until N) {
        for (j in 1 until totalCols) {
            nextpr[i, j] = nextpr[i, j] + nextpr[i, j - 1]
        }
    }

    // Generate samples
    return if (numTypes > 2) {
        // Multiple types: return array of [time, type] pairs
        val samples = Array(K) { IntArray(2) }
        for (n in 0 until K) {
            var time = 0

            // Play state transitions
            while (state < N) {
                // Geometric sojourn time
                time += 1 + floor(ln(random.nextDouble()) / logp[state]).toInt()

                // Choose next state
                val r = random.nextDouble()
                var nstate = 0
                while (nstate < totalCols - 1 && nextpr[state, nstate] <= r) {
                    nstate++
                }
                state = nstate
            }

            samples[n][0] = time
            samples[n][1] = state / N // Type (0-indexed)
            state = state % N
        }
        samples
    } else {
        // Single type: return just inter-arrival times
        val samples = IntArray(K)
        for (n in 0 until K) {
            var time = 0

            // Play state transitions
            while (state < N) {
                // Geometric sojourn time
                time += 1 + floor(ln(random.nextDouble()) / logp[state]).toInt()

                // Choose next state
                val r = random.nextDouble()
                var nstate = 0
                while (nstate < totalCols - 1 && nextpr[state, nstate] <= r) {
                    nstate++
                }
                state = nstate
            }

            samples[n] = time
            state = state % N
        }
        samples
    }
}

/**
 * Overload for Array<Matrix>.
 */
fun samplesFromDMMAP(
    D: Array<Matrix>,
    K: Int,
    initial: Int? = null,
    prec: Double = 1e-14,
    random: Random = Random()
): Any {
    val cell = MatrixCell(D.size)
    for (i in D.indices) {
        cell[i] = D[i]
    }
    return samplesFromDMMAP(cell, K, initial, prec, random)
}
