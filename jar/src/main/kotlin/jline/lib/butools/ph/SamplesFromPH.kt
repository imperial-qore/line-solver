/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph

import jline.util.matrix.Matrix
import java.util.Random
import kotlin.math.ln

/**
 * Generates random samples from a phase-type distribution.
 *
 * @param alpha The initial probability vector of the phase-type distribution.
 * @param A The transient generator matrix of the phase-type distribution.
 * @param K The number of samples to generate.
 * @param random Random number generator (optional).
 * @return The vector of random samples.
 */
fun samplesFromPH(alpha: Matrix, A: Matrix, K: Int, random: Random = Random()): DoubleArray {
    val N = alpha.length()

    // Cumulative initial distribution
    val cummInitial = DoubleArray(N)
    var cumSum = 0.0
    for (i in 0 until N) {
        cumSum += alpha[i]
        cummInitial[i] = cumSum
    }

    // Sojourn times
    val sojourn = DoubleArray(N)
    for (i in 0 until N) {
        sojourn[i] = -1.0 / A[i, i]
    }

    // Next state probabilities (cumulative)
    val nextpr = Array(N) { DoubleArray(N + 1) }
    for (i in 0 until N) {
        var rowSum = 0.0
        for (j in 0 until N) {
            if (i != j) {
                val prob = sojourn[i] * A[i, j]
                rowSum += prob
                nextpr[i][j] = rowSum
            } else {
                nextpr[i][j] = rowSum
            }
        }
        // Probability of absorption
        nextpr[i][N] = 1.0
    }

    val samples = DoubleArray(K)
    for (n in 0 until K) {
        var time = 0.0

        // Draw initial state from initial distribution
        val r1 = random.nextDouble()
        var state = 0
        while (state < N - 1 && cummInitial[state] <= r1) {
            state++
        }

        // Play state transitions until absorption
        while (state < N) {
            time -= ln(random.nextDouble()) * sojourn[state]
            val r2 = random.nextDouble()
            var nstate = 0
            while (nstate < N && nextpr[state][nstate] <= r2) {
                nstate++
            }
            state = nstate
        }
        samples[n] = time
    }

    return samples
}

/**
 * Overload for DoubleArray alpha.
 */
fun samplesFromPH(alpha: DoubleArray, A: Matrix, K: Int, random: Random = Random()): DoubleArray {
    return samplesFromPH(Matrix(alpha), A, K, random)
}
