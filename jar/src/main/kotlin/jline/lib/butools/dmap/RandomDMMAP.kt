/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.lib.butools.mc.dtmcSolve
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import java.util.Random
import kotlin.math.abs

/**
 * Returns a random discrete marked Markovian arrival process.
 *
 * @param order The size of the DMMAP
 * @param types The number of different arrival types
 * @param mean The mean inter-arrival time of the DMMAP (default: 10.0)
 * @param zeroEntries The number of zero entries in the D matrices (default: 0)
 * @param maxTrials Maximum number of trials to find a proper DMMAP (default: 1000)
 * @param prec Numerical precision for checking irreducibility (default: 1e-7)
 * @param random Random number generator
 * @return The D0...Dtypes matrices of the DMMAP
 */
fun randomDMMAP(
    order: Int,
    types: Int,
    mean: Double = 10.0,
    zeroEntries: Int = 0,
    maxTrials: Int = 1000,
    prec: Double = 1e-7,
    random: Random = Random()
): MatrixCell {
    if (types < 1) {
        throw IllegalArgumentException("RandomDMMAP: 'types' must be positive integer!")
    }

    if (zeroEntries > (order + 1) * (order - 1) + types * (order * order - 1)) {
        throw IllegalArgumentException("RandomDMMAP: Too many zero entries requested!")
    }

    // Generate all possible zero distributions among rows
    val zeroDistr = allZeroDistributions(order, zeroEntries)

    var trials = 0
    while (trials < maxTrials) {
        val indices = (0 until zeroDistr.size).shuffled(random)

        for (zdix in indices) {
            val zDistr = zeroDistr[zdix]

            var bad = false
            for (z in zDistr) {
                if (z >= (types + 1) * order - 1) {
                    bad = true
                    break
                }
            }
            if (bad) {
                trials++
                continue
            }

            val B = Matrix(order, (types + 1) * order)
            for (i in 0 until order) {
                val rp = (0 until (types + 1) * order - 1).shuffled(random)
                val a = DoubleArray((types + 1) * order - 1)
                for (j in 0 until (types + 1) * order - 1 - zDistr[i]) {
                    a[rp[j]] = random.nextDouble()
                }
                for (j in 0 until i) {
                    B[i, j] = a[j]
                }
                for (j in i + 1 until (types + 1) * order) {
                    B[i, j] = a[j - 1]
                }
            }

            val D = MatrixCell(types + 1)
            val sc = DoubleArray(order)
            for (k in 0..types) {
                D[k] = Matrix(order, order)
                for (i in 0 until order) {
                    for (j in 0 until order) {
                        D[k][i, j] = B[i, k * order + j]
                        sc[i] += D[k][i, j]
                    }
                }
            }

            if (sc.any { it == 0.0 }) continue

            for (k in 0..types) {
                for (i in 0 until order) {
                    for (j in 0 until order) {
                        D[k][i, j] = D[k][i, j] / sc[i]
                    }
                }
            }

            var sumD = Matrix(order, order)
            for (k in 0..types) {
                sumD = sumD.add(D[k])
            }

            val I = Matrix.eye(order)
            if (D[0].rank() == order && I.sub(sumD).rank() == order - 1) {
                val alpha = dtmcSolve(sumD)
                var minAlpha = Double.POSITIVE_INFINITY
                for (j in 0 until order) {
                    minAlpha = minOf(minAlpha, abs(alpha[0, j]))
                }

                if (minAlpha > prec) {
                    var fullZero = false
                    for (k in 0..types) {
                        var allZero = true
                        outer@ for (i in 0 until order) {
                            for (j in 0 until order) {
                                if (D[k][i, j] != 0.0) {
                                    allZero = false
                                    break@outer
                                }
                            }
                        }
                        if (allZero) {
                            fullZero = true
                            break
                        }
                    }

                    if (!fullZero) {
                        val d = DoubleArray(order) { random.nextDouble() }

                        val Dv = MatrixCell(types + 1)
                        for (k in 0..types) {
                            Dv[k] = Matrix(order, order)
                            for (i in 0 until order) {
                                for (j in 0 until order) {
                                    Dv[k][i, j] = (1.0 - d[i]) * D[k][i, j]
                                }
                            }
                        }
                        for (i in 0 until order) {
                            Dv[0][i, i] += d[i]
                        }

                        val m = marginalMomentsFromDMMAP(Dv, 1, prec)

                        for (i in 0 until order) {
                            d[i] = 1.0 - (1.0 - d[i]) * m[0] / mean
                        }

                        for (k in 0..types) {
                            for (i in 0 until order) {
                                for (j in 0 until order) {
                                    D[k][i, j] = (1.0 - d[i]) * D[k][i, j]
                                }
                            }
                        }
                        for (i in 0 until order) {
                            D[0][i, i] += d[i]
                        }

                        if (checkDMMAPRepresentation(D, prec)) {
                            return D
                        }
                    }
                }
            }
            trials++
        }
    }

    throw IllegalStateException("No feasible random DMMAP found!")
}

private fun allZeroDistributions(states: Int, zeros: Int): List<IntArray> {
    if (states == 1) {
        return listOf(intArrayOf(zeros))
    }
    val result = mutableListOf<IntArray>()
    for (iz in 0..zeros) {
        val subDistributions = allZeroDistributions(states - 1, zeros - iz)
        for (sub in subDistributions) {
            val combined = (sub.toList() + iz).sorted().toIntArray()
            if (result.none { it.contentEquals(combined) }) {
                result.add(combined)
            }
        }
    }
    return result
}
