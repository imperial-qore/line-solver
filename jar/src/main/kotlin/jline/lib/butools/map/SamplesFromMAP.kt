/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map

import jline.lib.butools.mc.drpSolve
import jline.util.matrix.Matrix
import java.util.Random
import kotlin.math.ln

/**
 * Generates random samples from a Markovian arrival process.
 *
 * @param D0 The D0 matrix of the Markovian arrival process
 * @param D1 The D1 matrix of the Markovian arrival process
 * @param K The number of samples to generate
 * @param initial Optional initial state (null means steady-state)
 * @param random Random number generator (optional)
 * @return The vector of random samples (inter-arrival times)
 */
fun samplesFromMAP(
    D0: Matrix,
    D1: Matrix,
    K: Int,
    initial: Int? = null,
    random: Random = Random()
): DoubleArray {
    val N = D0.numRows

    // Compute stationary distribution
    val P = D0.neg().inv().mult(D1)
    val pi = drpSolve(P)

    // Cumulative initial distribution
    val cummInitial = DoubleArray(N)
    var cumSum = 0.0
    for (i in 0 until N) {
        cumSum += pi[0, i]
        cummInitial[i] = cumSum
    }

    // Sojourn times
    val sojourn = DoubleArray(N)
    for (i in 0 until N) {
        sojourn[i] = -1.0 / D0[i, i]
    }

    // Transition probabilities within D0 (cumulative)
    val nextprD0 = Array(N) { DoubleArray(N) }
    for (i in 0 until N) {
        var rowSum = 0.0
        for (j in 0 until N) {
            if (i != j) {
                val prob = sojourn[i] * D0[i, j]
                rowSum += prob
                nextprD0[i][j] = rowSum
            } else {
                nextprD0[i][j] = rowSum
            }
        }
    }

    // Probability of transitioning via D1 (arrival)
    val arrivalProb = DoubleArray(N)
    for (i in 0 until N) {
        var d1sum = 0.0
        for (j in 0 until N) {
            d1sum += D1[i, j]
        }
        arrivalProb[i] = sojourn[i] * d1sum
    }

    // Transition probabilities within D1 (cumulative, normalized)
    val nextprD1 = Array(N) { DoubleArray(N) }
    for (i in 0 until N) {
        var rowSum = 0.0
        for (j in 0 until N) {
            val prob = D1[i, j]
            rowSum += prob
            nextprD1[i][j] = rowSum
        }
        // Normalize
        if (rowSum > 0) {
            for (j in 0 until N) {
                nextprD1[i][j] /= rowSum
            }
        }
    }

    val samples = DoubleArray(K)

    // Draw initial state
    var state = if (initial != null) {
        initial
    } else {
        val r = random.nextDouble()
        var s = 0
        while (s < N - 1 && cummInitial[s] <= r) {
            s++
        }
        s
    }

    for (n in 0 until K) {
        var time = 0.0

        // Transitions until arrival
        while (true) {
            time -= ln(random.nextDouble()) * sojourn[state]
            val r = random.nextDouble()

            // Check if this is an arrival (transition via D1)
            if (r < arrivalProb[state]) {
                // Transition via D1 - arrival occurs
                val r2 = random.nextDouble()
                var nstate = 0
                while (nstate < N - 1 && nextprD1[state][nstate] <= r2) {
                    nstate++
                }
                state = nstate
                break
            } else {
                // Transition via D0 - no arrival
                val r2 = random.nextDouble()
                var nstate = 0
                while (nstate < N - 1 && nextprD0[state][nstate] <= r2) {
                    nstate++
                }
                state = nstate
            }
        }
        samples[n] = time
    }

    return samples
}
