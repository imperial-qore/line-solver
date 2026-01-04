/**
 * @file Single Queue Network Interpolation (SQNI) approximate MVA
 * 
 * Implements the Single Queue Network Interpolation method for analyzing multi-class closed
 * queueing networks. Provides efficient approximation by reducing multi-class networks to
 * single-queue representations with interpolation-based corrections.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.mva

import jline.util.matrix.Matrix
import kotlin.math.sqrt

data class PfqnSqniResult(val Q: Matrix, val U: Matrix, val X: Matrix)

fun pfqn_sqni(N: Matrix, L: Matrix, Z: Matrix): PfqnSqniResult {
    val queueIdx = 0
    val C = L.length()
    val Nt = N.elementSum()
    val Q = Matrix(2, C)
    val U = Matrix(2, C)
    val X = Matrix(1, C)

    if (Nt <= 0) {
        return PfqnSqniResult(Q, U, X)
    }

    if (Nt == 1.0) {
        for (r in 0 until C) {
            val Xr = N.get(0, r) / (Z.get(r) + L.get(r))
            X.set(0, r, Xr)
            U.set(queueIdx, r, Xr * L.get(r))
            Q.set(queueIdx, r, Xr * L.get(r))
        }
    } else {
        for (r in 0 until C) {
            val Nr = N.get(0, r)
            val Lr = L.get(r)
            val Zr = Z.get(r)

            val Nvec_1r = N.copy()
            Nvec_1r.set(0, r, Nvec_1r.get(0, r) - 1)

            val sumN = N.elementSum()
            var sumBrPart = 0.0
            for (i in 0 until C) {
                if (i != r) {
                    val Zi = Z.get(i)
                    val Li = L.get(i)
                    val Ni = Nvec_1r.get(0, i)
                    sumBrPart += Zi * Ni / (Zi + Li + Li * (sumN - 2))
                }
            }

            val BrVec = Matrix(1, C)
            for (i in 0 until C) {
                val Zi = Z.get(i)
                val Li = L.get(i)
                val Ni = N.get(0, i)
                val denom = Zi + Li + Li * (sumN - 1 - sumBrPart)
                BrVec.set(0, i, Ni / denom * Zi)
            }

            var BrSum = 0.0
            for (i in 0 until C) {
                if (i != r) {
                    BrSum += BrVec.get(0, i)
                }
            }

            val Br = Lr * BrSum
            val Xr = if (Lr == 0.0) {
                Nr / Zr
            } else {
                val sqrtTerm =
                    sqrt(Br * Br - 2 * Br * Lr * Nt - 2 * Br * Zr + Lr * Lr * Nt * Nt + 2 * Lr * Nt * Zr - 4 * Nr * Lr * Zr + Zr * Zr)
                (Zr - sqrtTerm - Br + Lr * Nt) / (2 * Lr * Zr)
            }

            X.set(0, r, Xr)
            U.set(queueIdx, r, Xr * Lr)
            Q.set(queueIdx, r, Nr - Xr * Zr)
        }
    }

    for (r in 0 until C) {
        if (Z.get(r) == 0.0) {
            val denom = L.get(r) * (1 + Q.sumRows(queueIdx))
            val Xr = N.get(0, r) / denom
            X.set(0, r, Xr)
            U.set(queueIdx, r, Xr * L.get(r))
            Q.set(queueIdx, r, N.get(0, r) - Xr * Z.get(r))
        }
    }

    return PfqnSqniResult(Q, U, X)
}
/**
 * PFQN sqni algorithms
 */
@Suppress("unused")
class PfqnSqniAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}