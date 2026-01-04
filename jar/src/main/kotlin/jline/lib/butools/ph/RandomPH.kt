/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph

import jline.lib.butools.mc.ctmcSolve
import jline.util.matrix.Matrix
import java.util.Random
import kotlin.math.abs
import kotlin.math.min

/**
 * Result class for RandomPH.
 */
data class PHRepresentation(val alpha: Matrix, val A: Matrix)

/**
 * Returns a random phase-type distribution with a given order.
 *
 * @param order The size of the phase-type distribution.
 * @param mean The mean of the phase-type distribution (default 1.0).
 * @param zeroEntries The number of zero entries in the initial vector,
 *        generator matrix and closing vector (default 0).
 * @param maxTrials The maximum number of trials to find a proper PH
 *        (that has an irreducible phase process and none of its parameters
 *        is all-zero). The default value is 1000.
 * @param prec Numerical precision for checking the irreducibility.
 *        The default value is 1e-7.
 * @param random Random number generator (optional).
 * @return The PHRepresentation containing alpha and A.
 * @throws IllegalArgumentException if no feasible PH is found.
 */
fun randomPH(
    order: Int,
    mean: Double = 1.0,
    zeroEntries: Int = 0,
    maxTrials: Int = 1000,
    prec: Double = 1e-7,
    random: Random = Random()
): PHRepresentation {
    if (zeroEntries > (order + 1) * (order - 1)) {
        throw IllegalArgumentException("RandomPH: Too many zero entries requested!")
    }

    // Generate all possible zero distributions
    val zeroDistr = allZeroDistr(order, zeroEntries)

    var trials = 0
    while (trials < maxTrials) {
        // Shuffle the zero distributions
        val zdixList = zeroDistr.shuffled(random)

        for (zDistr in zdixList) {
            val B = Matrix.zeros(order, order + 2)
            for (i in 0 until order) {
                val rp = (0 until order + 1).shuffled(random)
                val a = DoubleArray(order + 1)
                for (j in 0 until order + 1 - zDistr[i]) {
                    a[rp[j]] = random.nextDouble()
                }
                // Fill B matrix
                for (j in 0 until i) {
                    B[i, j] = a[j]
                }
                for (j in i until order + 1) {
                    B[i, j + 1] = a[j]
                }
            }

            // Construct PH parameters
            val A = Matrix.zeros(order, order)
            for (i in 0 until order) {
                for (j in 0 until order) {
                    A[i, j] = B[i, j]
                }
            }
            val aVec = DoubleArray(order) { i -> B[i, order + 1] }
            val alphaVec = DoubleArray(order) { i -> B[i, order] }

            // A = A - diag(sum(A,2) + a)
            for (i in 0 until order) {
                var rowSum = 0.0
                for (j in 0 until order) {
                    rowSum += A[i, j]
                }
                A[i, i] = A[i, i] - rowSum - aVec[i]
            }

            // Check if all zero
            var allAZero = true
            var allAlphaZero = true
            var allAVecZero = true
            for (i in 0 until order) {
                for (j in 0 until order) {
                    if (A[i, j] != 0.0) allAZero = false
                }
                if (alphaVec[i] != 0.0) allAlphaZero = false
                if (aVec[i] != 0.0) allAVecZero = false
            }

            if (allAZero || allAlphaZero || allAVecZero) {
                continue
            }

            // Normalize alpha
            var alphaSum = 0.0
            for (i in 0 until order) {
                alphaSum += alphaVec[i]
            }
            val alpha = Matrix(1, order)
            for (i in 0 until order) {
                alpha[0, i] = alphaVec[i] / alphaSum
            }

            // D = A + a * alpha
            val D = Matrix.zeros(order, order)
            for (i in 0 until order) {
                for (j in 0 until order) {
                    D[i, j] = A[i, j] + aVec[i] * alpha[0, j]
                }
            }

            // Check if rank is order - 1 (irreducible)
            if (D.rank() == order - 1) {
                val pi = ctmcSolve(D)
                var minPi = Double.MAX_VALUE
                for (i in 0 until pi.length()) {
                    if (abs(pi[i]) < minPi) minPi = abs(pi[i])
                }

                if (minPi > prec) {
                    // Scale to the mean value
                    val moms = momentsFromPH(alpha, A, 1)
                    val scaleFactor = moms[0] / mean
                    val scaledA = A.scale(scaleFactor)
                    return PHRepresentation(alpha, scaledA)
                }
            }
            trials++
        }
    }

    throw IllegalArgumentException("No feasible random PH found with such many zero entries!")
}

/**
 * Helper function to generate all zero distributions.
 */
private fun allZeroDistr(states: Int, zeros: Int): List<IntArray> {
    if (states == 1) {
        return listOf(intArrayOf(zeros))
    }

    val result = mutableListOf<IntArray>()
    for (iz in 0..zeros) {
        val subResults = allZeroDistr(states - 1, zeros - iz)
        for (subResult in subResults) {
            val xt = (subResult.toList() + iz).sorted().toIntArray()
            // Check if we have it already
            var found = false
            for (existing in result) {
                if (existing.contentEquals(xt)) {
                    found = true
                    break
                }
            }
            if (!found) {
                result.add(xt)
            }
        }
    }
    return result
}
