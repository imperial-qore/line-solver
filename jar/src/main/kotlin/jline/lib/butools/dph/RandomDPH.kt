/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph

import jline.util.matrix.Matrix
import java.util.Random

/**
 * Returns a random discrete phase-type distribution with a given mean value.
 *
 * @param order The size of the discrete phase-type distribution
 * @param mean The mean of the discrete phase-type distribution (default 10.0)
 * @param zeroEntries The number of zero entries in the initial vector, generator matrix and closing vector (default 0)
 * @param maxTrials The maximum number of trials to find a proper DPH (default 1000)
 * @param prec Numerical precision for checking the irreducibility (default 1e-7)
 * @param random Random number generator (optional, defaults to a new Random instance)
 * @return The DPH2Representation or DPH3Representation containing alpha (initial probability vector) and A (transition probability matrix)
 *
 * Note: If the procedure fails, try to increase the 'maxTrials' parameter, or increase the mean value.
 */
fun randomDPH(
    order: Int,
    mean: Double = 10.0,
    zeroEntries: Int = 0,
    maxTrials: Int = 1000,
    prec: Double = 1e-7,
    random: Random = Random()
): MGRepresentation {
    if (zeroEntries > (order + 1) * (order - 1)) {
        throw IllegalArgumentException("RandomDPH: Too many zero entries requested!")
    }

    // Generate all possible distributions of zero entries among rows
    val zeroDistr = allZeroDistr(order, zeroEntries)

    var trials = 0
    while (trials < maxTrials) {
        // Randomly select a configuration
        val zdixList = (0 until zeroDistr.size).shuffled(random)

        for (zdix in zdixList) {
            val zDistr = zeroDistr[zdix]

            // Create B matrix (order x (order+2))
            val B = Matrix.zeros(order, order + 2)

            for (i in 0 until order) {
                val rp = (0 until order + 1).shuffled(random)
                val a = DoubleArray(order + 1)

                for (j in 0 until order + 1 - zDistr[i]) {
                    a[rp[j]] = random.nextDouble()
                }

                // B(i, 1:i-1) = a(1:i-1)
                for (j in 0 until i) {
                    B[i, j] = a[j]
                }
                // B(i, i+1:end) = a(i:end)
                for (j in i until order + 1) {
                    B[i, j + 1] = a[j]
                }
            }

            // Construct DPH parameters
            var A = Matrix.zeros(order, order)
            for (i in 0 until order) {
                for (j in 0 until order) {
                    A[i, j] = B[i, j]
                }
            }

            val aVec = Matrix(order, 1)
            for (i in 0 until order) {
                aVec[i, 0] = B[i, order + 1]
            }

            // Scale rows
            val sc = DoubleArray(order)
            var hasZeroRow = false
            for (i in 0 until order) {
                var rowSum = 0.0
                for (j in 0 until order) {
                    rowSum += A[i, j]
                }
                rowSum += aVec[i, 0]
                sc[i] = rowSum
                if (rowSum == 0.0) {
                    hasZeroRow = true
                    break
                }
            }

            if (hasZeroRow) continue

            for (i in 0 until order) {
                for (j in 0 until order) {
                    A[i, j] /= sc[i]
                }
                aVec[i, 0] /= sc[i]
            }

            // Extract alpha
            var alpha = Matrix(1, order)
            for (i in 0 until order) {
                alpha[0, i] = B[i, order]
            }

            // Check for all-zero matrices
            if (A.elementMax() == 0.0 || alpha.elementMax() == 0.0 || aVec.elementMax() == 0.0) {
                continue
            }

            // Normalize alpha
            val alphaSum = alpha.elementSum()
            if (alphaSum == 0.0) continue
            alpha = alpha.scale(1.0 / alphaSum)

            // Check irreducibility: rank(I - A) should be order
            val I = Matrix.eye(order)
            val IminusA = I.sub(A)
            val rank = IminusA.rank()

            if (rank == order) {
                // Check if alpha * inv(I - A) has all positive elements
                val invIminusA = IminusA.inv()
                val alphaInv = alpha.mult(invIminusA)
                var allPositive = true
                for (i in 0 until order) {
                    if (kotlin.math.abs(alphaInv[0, i]) <= prec) {
                        allPositive = false
                        break
                    }
                }

                if (allPositive) {
                    // Scale diagonals to achieve target mean
                    val d = DoubleArray(order) { random.nextDouble() }

                    // Compute current mean with diagonal modification
                    var scaledA = Matrix.zeros(order, order)
                    for (i in 0 until order) {
                        for (j in 0 until order) {
                            if (i == j) {
                                scaledA[i, j] = d[i]
                            } else {
                                scaledA[i, j] = (1 - d[i]) * A[i, j]
                            }
                        }
                    }

                    val testMoms = momentsFromDPH(alpha, scaledA, 1)
                    val currentMean = testMoms[0]

                    // Scale to target mean
                    for (i in 0 until order) {
                        d[i] = 1 - (1 - d[i]) * currentMean / mean
                    }

                    // Rebuild A with scaled diagonals
                    for (i in 0 until order) {
                        for (j in 0 until order) {
                            if (i == j) {
                                A[i, j] = d[i]
                            } else {
                                A[i, j] = (1 - d[i]) * A[i, j]
                            }
                        }
                    }

                    if (checkDPHRepresentation(alpha, A, prec)) {
                        return MGRepresentation(alpha, A)
                    }
                }
            }

            trials++
        }
    }

    throw IllegalArgumentException("RandomDPH: No feasible random DPH found! Try increasing mean or maxTrials.")
}

/**
 * Generate all possible distributions of zero entries among rows.
 */
private fun allZeroDistr(states: Int, zeros: Int): List<IntArray> {
    if (states == 1) {
        return listOf(intArrayOf(zeros))
    }

    val result = mutableListOf<IntArray>()
    for (iz in 0..zeros) {
        val subDistrs = allZeroDistr(states - 1, zeros - iz)
        for (subDistr in subDistrs) {
            val combined = (subDistr.toList() + iz).sorted().toIntArray()
            // Check if we already have this combination
            val found = result.any { existing ->
                existing.size == combined.size && existing.indices.all { existing[it] == combined[it] }
            }
            if (!found) {
                result.add(combined)
            }
        }
    }
    return result
}
