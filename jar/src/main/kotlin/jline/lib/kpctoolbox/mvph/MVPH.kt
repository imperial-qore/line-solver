package jline.lib.kpctoolbox.mvph

import jline.util.matrix.Matrix
import org.apache.commons.math3.linear.LUDecomposition
import org.apache.commons.math3.linear.MatrixUtils

/**
 * Multivariate Phase-Type Distribution (MVPH) functions.
 *
 * Ported from MATLAB: matlab/lib/kpctoolbox/mvph/
 */

/**
 * Computes the joint moment E[X^n1 * Y^n2] of a bivariate phase-type distribution.
 *
 * The bivariate phase-type distribution has:
 * - Initial vector alpha
 * - First phase transition matrix S (for X)
 * - Transition matrix D between phases (from X completion to Y start)
 * - Second phase transition matrix T (for Y)
 *
 * Formula: E[X^n1 * Y^n2] = n1! * n2! * alpha * inv(-S)^(n1+1) * D * inv(-T)^(n2+1) * (-T) * e
 *
 * @param alpha Initial probability vector
 * @param S First phase generator matrix (for X)
 * @param T Second phase generator matrix (for Y)
 * @param D Transition matrix between phases
 * @param n1 Power for first variable X
 * @param n2 Power for second variable Y
 * @return Joint moment E[X^n1 * Y^n2]
 */
fun mvph_joint(
    alpha: DoubleArray,
    S: Matrix,
    T: Matrix,
    D: Matrix,
    n1: Int,
    n2: Int
): Double {
    val sSize = S.numRows
    val tSize = T.numRows

    // Convert S to Apache Commons Math format and compute inv(-S)
    val negS = MatrixUtils.createRealMatrix(sSize, sSize)
    for (i in 0 until sSize) {
        for (j in 0 until sSize) {
            negS.setEntry(i, j, -S.get(i, j))
        }
    }

    // Compute inv(-S)
    val invNegS = LUDecomposition(negS).solver.inverse

    // Compute inv(-S)^(n1+1)
    var invNegSPow = MatrixUtils.createRealIdentityMatrix(sSize)
    for (p in 0 until n1 + 1) {
        invNegSPow = invNegSPow.multiply(invNegS)
    }

    // Convert T to Apache Commons Math format and compute inv(-T)
    val negT = MatrixUtils.createRealMatrix(tSize, tSize)
    for (i in 0 until tSize) {
        for (j in 0 until tSize) {
            negT.setEntry(i, j, -T.get(i, j))
        }
    }

    // Compute inv(-T)
    val invNegT = LUDecomposition(negT).solver.inverse

    // Compute inv(-T)^(n2+1)
    var invNegTPow = MatrixUtils.createRealIdentityMatrix(tSize)
    for (p in 0 until n2 + 1) {
        invNegTPow = invNegTPow.multiply(invNegT)
    }

    // Convert D to Apache Commons Math format
    val dMatrix = MatrixUtils.createRealMatrix(sSize, tSize)
    for (i in 0 until sSize) {
        for (j in 0 until tSize) {
            dMatrix.setEntry(i, j, D.get(i, j))
        }
    }

    // Compute (-T) - note: negT is already -T
    // But we need the actual -T matrix for the final multiplication
    val minusT = MatrixUtils.createRealMatrix(tSize, tSize)
    for (i in 0 until tSize) {
        for (j in 0 until tSize) {
            minusT.setEntry(i, j, -T.get(i, j))
        }
    }

    // Compute the chain: alpha * inv(-S)^(n1+1) * D * inv(-T)^(n2+1) * (-T) * e

    // Step 1: alpha * inv(-S)^(n1+1)
    val alphaVec = MatrixUtils.createRowRealMatrix(alpha)
    val step1 = alphaVec.multiply(invNegSPow)

    // Step 2: step1 * D
    val step2 = step1.multiply(dMatrix)

    // Step 3: step2 * inv(-T)^(n2+1)
    val step3 = step2.multiply(invNegTPow)

    // Step 4: step3 * (-T)
    val step4 = step3.multiply(minusT)

    // Step 5: step4 * e (column vector of ones)
    val ones = MatrixUtils.createColumnRealMatrix(DoubleArray(tSize) { 1.0 })
    val step5 = step4.multiply(ones)

    // Get the scalar result
    val result = step5.getEntry(0, 0)

    // Multiply by factorials
    return factorial(n1) * factorial(n2) * result
}

/**
 * Computes the mean of the first variable in a bivariate PH distribution.
 * E[X] = mvph_joint(alpha, S, T, D, 1, 0)
 */
fun mvph_mean_x(alpha: DoubleArray, S: Matrix, T: Matrix, D: Matrix): Double {
    return mvph_joint(alpha, S, T, D, 1, 0)
}

/**
 * Computes the mean of the second variable in a bivariate PH distribution.
 * E[Y] = mvph_joint(alpha, S, T, D, 0, 1)
 */
fun mvph_mean_y(alpha: DoubleArray, S: Matrix, T: Matrix, D: Matrix): Double {
    return mvph_joint(alpha, S, T, D, 0, 1)
}

/**
 * Computes the covariance of a bivariate PH distribution.
 * Cov(X, Y) = E[XY] - E[X]*E[Y]
 */
fun mvph_cov(alpha: DoubleArray, S: Matrix, T: Matrix, D: Matrix): Double {
    val exy = mvph_joint(alpha, S, T, D, 1, 1)
    val ex = mvph_joint(alpha, S, T, D, 1, 0)
    val ey = mvph_joint(alpha, S, T, D, 0, 1)
    return exy - ex * ey
}

/**
 * Computes the correlation of a bivariate PH distribution.
 * Corr(X, Y) = Cov(X, Y) / (StdDev(X) * StdDev(Y))
 */
fun mvph_corr(alpha: DoubleArray, S: Matrix, T: Matrix, D: Matrix): Double {
    val ex = mvph_joint(alpha, S, T, D, 1, 0)
    val ey = mvph_joint(alpha, S, T, D, 0, 1)
    val ex2 = mvph_joint(alpha, S, T, D, 2, 0)
    val ey2 = mvph_joint(alpha, S, T, D, 0, 2)
    val exy = mvph_joint(alpha, S, T, D, 1, 1)

    val varX = ex2 - ex * ex
    val varY = ey2 - ey * ey
    val covXY = exy - ex * ey

    if (varX <= 0 || varY <= 0) {
        return 0.0
    }

    return covXY / (kotlin.math.sqrt(varX) * kotlin.math.sqrt(varY))
}

/**
 * Computes factorial.
 */
private fun factorial(n: Int): Double {
    if (n <= 1) return 1.0
    var result = 1.0
    for (i in 2..n) {
        result *= i.toDouble()
    }
    return result
}
