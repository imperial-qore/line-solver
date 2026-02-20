/**
 * @file Knessl-Tier asymptotic ray method for normalizing constant computation
 * 
 * Implements the Knessl-Tier asymptotic expansion using the ray method for computing
 * normalizing constants in closed product-form queueing networks. Provides high-accuracy
 * asymptotic approximation for large population systems with self-loop preprocessing.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.io.Ret
import jline.GlobalConstants
import jline.api.pfqn.mva.pfqn_aql
import jline.api.pfqn.mva.pfqn_bs
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Knessl-Tier asymptotic expansion of the normalizing constant using the ray method.
 *
 * @param L - demands at all stations
 * @param N - number of jobs for each class
 * @param Z - think time for each class
 * @return normalizing constant (G) and its logarithm (lG)
 */

fun pfqn_kt(L: Matrix, N: Matrix, Z: Matrix): Ret.pfqnNc {/* ---------- trivial / degenerate cases ---------- */

    var L = L
    var N = N
    var Z = Z
    val Ntot0 = N.sumRows(0)
    val out = Ret.pfqnNc(0.0, 0.0)
    if (L.isEmpty || N.isEmpty || Ntot0 <= GlobalConstants.Zero) {
        out.G = 1.0
        return out
    }
    if (Z.isEmpty) {
        Z = Matrix(1, L.numCols) // 1×R vector of zeros
    }

    /* ---------- pre-processing: remove self-loop classes ---------- */
    val Morig = L.numRows
    val Rorig = L.numCols
    var slcDemandFactor = 0.0

    if (Rorig > 1) {
        val isSLC = BooleanArray(Rorig)

        for (r in 0..<Rorig) {
            var nnz = 0
            var firstRow = -1
            for (k in 0..<Morig) {
                if (L[k, r] > GlobalConstants.Zero) {
                    nnz++
                    firstRow = k
                }
            }

            /* single-station class with zero think-time → self-loop customer */
            if (nnz == 1 && Z[0, r] == 0.0) {/* replicate that single row N(r) times and append at bottom   */

                val rep = Math.round(N[0, r]).toInt()
                if (rep > 0) {
                    val block = Matrix(rep, Rorig)
                    for (i in 0..<rep) {
                        for (c in 0..<Rorig) {
                            block[i, c] = L[firstRow, c]
                        }
                    }
                    L = Matrix.concatRows(L, block, null) // <-- Matrix API
                }

                isSLC[r] = true
                slcDemandFactor = N[0, r] * FastMath.log(L[firstRow, r]) // see MATLAB code
            }
        }

        /* actually drop those self-loop columns from L, N, Z */
        var remaining = 0
        for (b in isSLC) if (!b) remaining++

        val newL = Matrix(L.numRows, remaining)
        val newN = Matrix(1, remaining)
        val newZ = Matrix(1, remaining)

        var c = 0
        for (r in 0..<Rorig) {
            if (isSLC[r]) continue
            for (k in 0..<L.numRows) {
                newL[k, c] = L[k, r]
            }
            newN[0, c] = N[0, r]
            newZ[0, c] = Z[0, r]
            c++
        }
        L = newL
        N = newN
        Z = newZ
    }

    /* ---------- core ray-method computation ---------- */
    val M = L.numRows
    val R = L.numCols
    val Ntot = N.sumRows(0)

    val beta = DoubleArray(R)
    for (r in 0..<R) {
        beta[r] = N[0, r] / Ntot
    }

    /* choose exact (BS) or asymptotic (AQL) solver for X,Q */
    val XQ = if (Ntot <= 4.0) pfqn_bs(L, N, Z) else pfqn_aql(L, N, Z)
    val X = XQ.X // 1 × R visit ratios

    /* build covariance-like matrix C */
    val delta = Matrix.eye(R)
    val C = Matrix(R, R)
    for (s in 0..<R) {
        for (r in 0..<R) {
            var SK = 0.0
            for (k in 0..<M) {
                var Uk = 0.0
                for (t in 0..<R) {
                    Uk += L[k, t] * X[0, t]
                }
                val denom = FastMath.max(GlobalConstants.FineTol, 1.0 - Uk)
                SK += X[0, s] * X[0, r] * L[k, s] * L[k, r] / (denom * denom)
            }
            C[s, r] = delta[s, r] * beta[s] + SK / Ntot
        }
    }

    /* denominator ∏ₖ max(ε,1−Uₖ) */
    var Den = 1.0
    for (k in 0..<M) {
        var Uk = 0.0
        for (r in 0..<R) {
            Uk += L[k, r] * X[0, r]
        }
        Den *= FastMath.max(GlobalConstants.FineTol, 1.0 - Uk)
    }

    /* log(G) */
    var lG = (FastMath.log(FastMath.pow(2.0 * FastMath.PI, -R / 2.0) / FastMath.sqrt(FastMath.pow(Ntot,
        R) * C.det())) - FastMath.log(Den) + slcDemandFactor)

    for (r in 0..<R) {
        lG += -Ntot * beta[r] * FastMath.log(X[0, r])
    }

    out.lG = lG
    out.G = FastMath.exp(lG)
    return out
}
/**
 * PFQN kt algorithms
 */
@Suppress("unused")
class PfqnKtAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}