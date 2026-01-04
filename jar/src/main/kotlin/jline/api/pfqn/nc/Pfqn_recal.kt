/**
 * @file RECAL normalizing constant method for closed queueing networks
 *
 * Implements the RECAL (Recursive Calculation) algorithm for computing normalizing constants
 * in closed product-form queueing networks. Provides exact computation using recursive
 * enumeration techniques with support for multiple job classes and station multiplicities.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.api.pfqn.pfqn_unique
import jline.api.pfqn.pfqn_combine_mi
import jline.io.Ret
import jline.util.Maths
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * RECAL method to compute the normalizing constant of a load-independent closed queueing network model
 *
 * @param L - demands at all stations
 * @param N - number of jobs for each class
 * @return normalizing constant (G), its logarithm (lG))
 */

fun pfqn_recal(L: Matrix, N: Matrix): Ret.pfqnNc {
    val R = L.numCols
    val Z = Matrix(1, R)
    return pfqn_recal(L, N, Z)
}

/**
 * RECAL method to compute the normalizing constant of a load-independent closed queueing network model
 *
 * @param L - demands at all stations
 * @param N - number of jobs for each class
 * @param Z - think time for each class
 * @return normalizing constant (G), its logarithm (lG))
 */

fun pfqn_recal(L: Matrix, N: Matrix, Z: Matrix): Ret.pfqnNc {
    val M = L.numRows
    val m0 = Matrix(1, M)
    m0.fill(1.0)
    return pfqn_recal(L, N, Z, m0)
}

/**
 * RECAL method to compute the normalizing constant of a load-independent closed queueing network model
 *
 * @param L  - demands at all stations
 * @param N  - number of jobs for each class
 * @param Z  - think time for each class
 * @param m0 - vector with multiplicities of each station
 * @return normalizing constant (G), its logarithm (lG))
 */

fun pfqn_recal(L: Matrix, N: Matrix, Z: Matrix, m0: Matrix): Ret.pfqnNc {
    val M_original = L.numRows
    val R = L.numCols

    // Normalize m0 to row vector
    var m0_local = m0
    if (m0_local.numRows > 1) {
        m0_local = m0_local.transpose()
    }

    // Detect and consolidate replicated stations
    val uniqueResult = pfqn_unique(L)
    val L_reduced = uniqueResult.L_unique
    val mapping = uniqueResult.mapping
    val M = L_reduced.numRows

    // Combine user-provided m0 with detected multiplicity
    m0_local = pfqn_combine_mi(m0_local, mapping, M)

    val Ntot = N.elementSum().toInt()

    // Fix floating point precision issue by rounding
    val G_1_size = (Maths.nCk((Ntot + (M + 1) - 1).toDouble(), Ntot.toDouble()) + 0.5).toInt()

    var G_1 = Matrix(1, G_1_size)
    G_1.fill(1.0)
    var G = G_1.copy()

    // Element M+1 is for the delay server
    var I_1 = Maths.multichoose((M + 1).toDouble(), Ntot.toDouble())
    var n = 0

    for (r in 0..<R) {
        for (nr in 1..N[0, r].toInt()) {
            n++
            val I = Maths.multichoose((M + 1).toDouble(), (Ntot + 1 - (n + 1)).toDouble())

            // G should be sized to match I.numRows
            G = Matrix(1, I.numRows)

            for (i in 0..<I.numRows) {
                val m = I.getRow(i)
                val mZ = Matrix.extractColumns(m, 0, M)  // extract first M elements like m(1:M) in MATLAB
                val I_1_subset = Matrix.extractColumns(I_1, 0, M)  // extract first M columns like I_1(:,1:M) in MATLAB
                val mzIndex = Matrix.matchrow(I_1_subset, mZ)

                G[i] = Z[0, r] * G_1[mzIndex] / nr

                for (j in 0..<M) {
                    // here m is modified
                    m[j] = m[j] + 1
                    val fullIndex = Matrix.matchrow(I_1, m)
                    G[i] = G[i] + (m[j] + m0_local[0, j] - 1) * L_reduced[j, r] * G_1[fullIndex] / nr
                    m[j] = m[j] - 1
                }
            }
            I_1 = I
            G_1 = G
        }
    }

    return Ret.pfqnNc(G[0], FastMath.log(G[0]))
}
/**
 * PFQN recal algorithms
 */
@Suppress("unused")
class PfqnRecalAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}