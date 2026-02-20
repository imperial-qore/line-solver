/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map

import jline.lib.butools.mc.ctmcSolve
import jline.util.matrix.Matrix
import java.util.Random
import kotlin.math.abs

/**
 * Result class for RandomMAP.
 */
data class MAPRepresentation(val D0: Matrix, val D1: Matrix)

/**
 * Returns a random Markovian arrival process with given mean value.
 *
 * @param order The size of the MAP
 * @param mean The mean inter-arrival times of the MAP (default 1.0)
 * @param zeroEntries The number of zero entries in the D0 and D1 matrices (default 0)
 * @param maxTrials The maximum number of trials to find a proper MAP (default 1000)
 * @param prec Numerical precision for checking the irreducibility (default 1e-7)
 * @param random Random number generator (optional)
 * @return The MAPRepresentation containing D0 and D1 matrices
 * @throws IllegalArgumentException if no feasible MAP is found
 */
@JvmOverloads
fun randomMAP(
    order: Int,
    mean: Double = 1.0,
    zeroEntries: Int = 0,
    maxTrials: Int = 1000,
    prec: Double = 1e-7,
    random: Random = Random()
): MAPRepresentation {
    if (zeroEntries > 2 * order * (order - 1)) {
        throw IllegalArgumentException("RandomMAP: Too many zero entries requested!")
    }

    // Generate all possible zero distributions
    val zeroDistr = allZeroDistr(order, zeroEntries)

    var trials = 0
    while (trials < maxTrials) {
        val zdixList = zeroDistr.shuffled(random)

        for (zDistr in zdixList) {
            // Create random D0 and D1 matrices
            val D0 = Matrix.zeros(order, order)
            val D1 = Matrix.zeros(order, order)

            // Fill with random values respecting zero entries
            val totalEntries = 2 * order * order
            val nonZeroEntries = totalEntries - zeroEntries

            // Randomly distribute entries
            val allPositions = mutableListOf<Pair<Int, Int>>()
            for (i in 0 until order) {
                for (j in 0 until order) {
                    allPositions.add(Pair(i, j))
                    allPositions.add(Pair(i, j + order)) // Second matrix
                }
            }

            val shuffledPositions = allPositions.shuffled(random)
            for (idx in 0 until nonZeroEntries.coerceAtMost(shuffledPositions.size)) {
                val (i, j) = shuffledPositions[idx]
                val value = random.nextDouble()
                if (j < order) {
                    if (i != j) { // Off-diagonal for D0
                        D0[i, j] = value
                    }
                } else {
                    D1[i, j - order] = value
                }
            }

            // Set diagonal of D0 to make row sums of D0+D1 equal to 0
            for (i in 0 until order) {
                var rowSum = 0.0
                for (j in 0 until order) {
                    if (i != j) rowSum += D0[i, j]
                    rowSum += D1[i, j]
                }
                D0[i, i] = -rowSum
            }

            // Check if valid MAP
            var allD0Zero = true
            var allD1Zero = true
            for (i in 0 until order) {
                for (j in 0 until order) {
                    if (D0[i, j] != 0.0) allD0Zero = false
                    if (D1[i, j] != 0.0) allD1Zero = false
                }
            }

            if (allD0Zero || allD1Zero) {
                continue
            }

            // Check if irreducible: D = D0 + D1 should have rank n-1
            val D = D0.add(D1)
            if (D.rank() == order - 1) {
                val pi = ctmcSolve(D)
                var minPi = Double.MAX_VALUE
                for (i in 0 until pi.length()) {
                    if (abs(pi[i]) < minPi) minPi = abs(pi[i])
                }

                if (minPi > prec) {
                    // Scale to the mean value
                    val moms = marginalMomentsFromMAP(D0, D1, 1)
                    val scaleFactor = moms[0] / mean
                    return MAPRepresentation(D0.scale(scaleFactor), D1.scale(scaleFactor))
                }
            }
            trials++
        }
    }

    throw IllegalArgumentException("No feasible random MAP found!")
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
