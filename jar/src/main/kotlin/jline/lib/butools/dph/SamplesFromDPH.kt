/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph

import jline.util.matrix.Matrix
import java.util.Random
import kotlin.math.floor
import kotlin.math.ln

/**
 * Generates random samples from a discrete phase-type distribution.
 *
 * @param alpha The initial probability vector of the discrete phase-type distribution.
 * @param A The transition probability matrix of the discrete phase-type distribution.
 * @param K The number of samples to generate.
 * @param random Random number generator (optional, defaults to a new Random instance)
 * @return The vector of random samples
 */
fun samplesFromDPH(alpha: Matrix, A: Matrix, K: Int, random: Random = Random()): IntArray {
    val N = alpha.length()

    // Compute cumulative initial distribution
    val cummInitial = DoubleArray(N)
    var sum = 0.0
    for (i in 0 until N) {
        sum += alpha[i]
        cummInitial[i] = sum
    }

    // logp = log(diag(A))
    val logp = DoubleArray(N)
    for (i in 0 until N) {
        logp[i] = if (A[i, i] > 0) ln(A[i, i]) else Double.NEGATIVE_INFINITY
    }

    // sojourn = 1 / (1 - diag(A))
    val sojourn = DoubleArray(N)
    for (i in 0 until N) {
        sojourn[i] = 1.0 / (1.0 - A[i, i])
    }

    // nextpr = diag(sojourn) * A - diag(diag(nextpr))
    // Then add column for absorption: [nextpr, 1 - sum(nextpr, 2)]
    // Then cumsum along rows
    val nextpr = Array(N) { DoubleArray(N + 1) }
    for (i in 0 until N) {
        var rowSum = 0.0
        for (j in 0 until N) {
            if (i != j) {
                nextpr[i][j] = sojourn[i] * A[i, j]
                rowSum += nextpr[i][j]
            }
        }
        nextpr[i][N] = 1.0 - rowSum
    }

    // Cumsum along rows
    val nextprCum = Array(N) { DoubleArray(N + 1) }
    for (i in 0 until N) {
        var cumSum = 0.0
        for (j in 0..N) {
            cumSum += nextpr[i][j]
            nextprCum[i][j] = cumSum
        }
    }

    // Generate samples
    val x = IntArray(K)
    for (n in 0 until K) {
        var time = 0

        // Draw initial state
        val r = random.nextDouble()
        var state = 0
        while (state < N && cummInitial[state] <= r) {
            state++
        }

        // Play state transitions until absorption
        while (state < N) {
            // Sojourn time in current state (geometric distribution)
            if (logp[state] > Double.NEGATIVE_INFINITY) {
                val sojournTime = 1 + floor(ln(random.nextDouble()) / logp[state]).toInt()
                time += sojournTime
            } else {
                time += 1
            }

            // Next state transition
            val r2 = random.nextDouble()
            var nstate = 0
            while (nstate <= N && nextprCum[state][nstate] <= r2) {
                nstate++
            }
            state = nstate
        }

        x[n] = time
    }

    return x
}

/**
 * Overload for DoubleArray alpha.
 */
fun samplesFromDPH(alpha: DoubleArray, A: Matrix, K: Int, random: Random = Random()): IntArray {
    return samplesFromDPH(Matrix(alpha), A, K, random)
}
