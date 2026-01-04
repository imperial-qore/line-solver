/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 *
 * Reference:
 * A. van de Liefvoort. The moment problem for continuous distributions.
 * Technical report, University of Missouri, WP-CM-1990-02, Kansas City, 1990.
 */
package jline.lib.butools.ph

import jline.lib.butools.ReducedMomsFromMoms
import jline.util.matrix.Matrix
import kotlin.math.ceil

/**
 * Result class for MEFromMoments containing both alpha and A.
 */
data class MERepresentation(val alpha: Matrix, val A: Matrix)

/**
 * Creates a matrix-exponential distribution that has the
 * same moments as given.
 *
 * @param moms The list of moments. The order of the resulting
 *        matrix-exponential distribution is determined based on
 *        the number of moments given. To obtain a matrix exponential
 *        distribution of order M, 2*M-1 moments are required.
 * @return The MERepresentation containing alpha (initial vector) and A (matrix parameter).
 */
fun meFromMoments(moms: DoubleArray): MERepresentation {
    val K = appie(ReducedMomsFromMoms(moms))
    val N = ceil(moms.size / 2.0).toInt()

    // Create T matrix (lower triangular with ones)
    val T = Matrix.zeros(N, N)
    for (i in 0 until N) {
        for (j in 0..i) {
            T[i, j] = 1.0
        }
    }

    // Create U matrix (upper triangular with 1/(N-i+1) values)
    val U = Matrix.zeros(N, N)
    for (i in 0 until N) {
        for (j in i until N) {
            U[i, j] = 1.0 / (N - i)
        }
    }

    // Compute alpha = [1, 0, 0, ...] * inv(T) * U
    val alphaInit = Matrix.zeros(1, N)
    alphaInit[0, 0] = 1.0
    val alpha = alphaInit.mult(T.inv()).mult(U)

    // Compute A = inv(-inv(U)*T*K*inv(T)*U)
    val Uinv = U.inv()
    val Tinv = T.inv()
    val inner = Uinv.neg().mult(T).mult(K).mult(Tinv).mult(U)
    val A = inner.inv()

    return MERepresentation(alpha, A)
}

/**
 * Implements the Appie algorithm for computing the K matrix from reduced moments.
 */
private fun appie(rmom: DoubleArray): Matrix {
    val m = rmom.size
    val actualM: Int
    val rm: DoubleArray

    if (m % 2 == 0) {
        rm = DoubleArray(m) { i -> if (i == 0) 1.0 else rmom[i - 1] }
        actualM = m / 2
    } else {
        rm = DoubleArray(m + 1) { i -> if (i == 0) 1.0 else rmom[i - 1] }
        actualM = ceil(m / 2.0).toInt()
    }

    val twoM = 2 * actualM
    val f = DoubleArray(twoM)
    f[0] = 1.0
    val y = DoubleArray(twoM)
    val dd = Matrix.zeros(twoM, twoM)

    var n = 0
    var k = 0
    var q = 1.0
    val d = IntArray(actualM)
    val alphaArr = Array(actualM) { DoubleArray(actualM) }
    val beta = DoubleArray(actualM)

    for (i in 1 until twoM) {
        dd[i, i - 1] = 1.0
    }

    for (i in 0 until twoM) {
        // Compute ro = q * rm * f
        var ro = 0.0
        for (j in 0 until rm.size.coerceAtMost(f.size)) {
            ro += q * rm[j] * f[j]
        }

        val nold = n
        n = nold + 1

        val yold = y.copyOf()

        if (n > 0 && ro != 0.0) {
            if (k > 0) {
                val power = d[k - 1] + n - 1
                beta[k - 1] = ro / Math.pow(rm[1], power.toDouble())
            }
            k++
            d[k - 1] = n
            n = -n
            q = q / ro

            // y = dd * f
            for (j in 0 until twoM) {
                y[j] = 0.0
                for (l in 0 until twoM) {
                    y[j] += dd[j, l] * f[l]
                }
            }
        } else if (n <= 0) {
            if (k > 0) {
                val j = nold + d[k - 1] + 1
                if (j > 0 && j <= actualM) {
                    alphaArr[k - 1][j - 1] = ro / Math.pow(rm[1], (j - 1).toDouble())
                }
            }
        }

        // f = dd * f - ro * yold
        val newF = DoubleArray(twoM)
        for (j in 0 until twoM) {
            for (l in 0 until twoM) {
                newF[j] += dd[j, l] * f[l]
            }
            newF[j] -= ro * yold[j]
        }
        for (j in 0 until twoM) {
            f[j] = newF[j]
        }
    }

    // Check sum of d
    var sumD = 0
    for (i in 0 until actualM) {
        sumD += d[i]
    }
    if (sumD != actualM) {
        throw IllegalArgumentException("MEFromMoments: Insufficient matrix order!")
    }

    // Construct K matrix
    val K = Matrix.zeros(actualM, actualM)
    K[0, 0] = rm[1]
    for (i in 0 until actualM - 1) {
        K[i, i + 1] = rm[1]
    }

    var ind = d[0]
    for (i in 1 until actualM) {
        if (ind < actualM) {
            val inc = d[i]
            ind += inc
            if (ind <= actualM) {
                K[ind - 1, ind - inc - d[i - 1]] = beta[i - 1]
                for (j in 0 until inc) {
                    K[ind - 1, ind - j - 1] = alphaArr[i][j]
                }
            }
        }
    }

    return K
}
