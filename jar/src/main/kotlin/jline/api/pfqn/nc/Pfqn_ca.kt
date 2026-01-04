/**
 * Convolution Algorithm for Product-Form Networks
 * 
 * Implements the classic Buzen convolution algorithm for computing normalizing constants
 * in closed product-form queueing networks. Provides exact results with computational
 * complexity exponential in the number of job classes.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.io.Ret
import jline.GlobalConstants
import jline.GlobalConstants.NegInf
import jline.util.Maths
import jline.util.PopulationLattice
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Compute the normalizing constant using the convolution algorithm
 *
 * @param L - demands at all stations
 * @param N - number of jobs for each class
 * @return normalizing constant and its logarithm
 */
fun pfqn_ca(L: Matrix, N: Matrix): Ret.pfqnNc {
    val Z = N.copy()
    Z.zero()
    return pfqn_ca(L, N, Z)
}

/**
 * Compute the normalizing constant using the convolution algorithm
 *
 * @param L - demands at all stations
 * @param N - number of jobs for each class
 * @param Z - think times
 * @return normalizing constant and its logarithm
 */
fun pfqn_ca(L: Matrix, N: Matrix, Z: Matrix = Matrix(1, L.numCols)): Ret.pfqnNc {
    var Z = Z
    val M = L.numRows
    val R = L.numCols

    if (M == 0) {
        // Calculate the logarithm of the factorial of every element in N
        val tmp = Matrix(1, N.length())
        for (i in 0..<N.length()) {
            tmp[0, i] = -Maths.factln(N[i])
        }
        var lGn = tmp.sumRows().sumCols().value()

        val tmp2 = Z!!.sumCols()
        for (i in 0..<tmp2.length()) {
            tmp2[i] = FastMath.log(tmp2[i])
        }
        if (N.length() == 1) {
            lGn += (N[0] * tmp2.sumRows()[0])
        } else if (tmp2.length() == 1) {
            lGn += (tmp2[0] * N.sumRows().sumCols().value())
        } else {
            val tmp3 = Matrix(1, N.length())
            for (i in 0..<N.length()) {
                tmp3[i] = N[i]
            }
            lGn += tmp3.elementMult(tmp2, null).sumRows()[0]
        }
        val Gn = FastMath.exp(lGn)

        return Ret.pfqnNc(Gn, lGn)
    }

    if (N.elementMin() < 0) {
        return Ret.pfqnNc(0.0, NegInf)
    }

    if (N.sumRows().sumCols()[0] == 0.0) {
        return Ret.pfqnNc(1.0, 0.0)
    }

    Z = if (Z.isEmpty) {
        val temp = Matrix(1, R)
        temp.fill(0.0)
        temp
    } else Z

    var product_N_plus_one = 1
    for (i in 0..<N.length()) {
        product_N_plus_one = (product_N_plus_one * (N[i] + 1)).toInt()
    }
    val G = Matrix(M + 1, product_N_plus_one)
    G.fill(1.0)
    var n = PopulationLattice.pprod(N)

    while (FastMath.abs(n.sumRows().sumCols()[0] + 1) > GlobalConstants.FineTol) {
        val idxn = PopulationLattice.hashpop(n, N)
        G[0, idxn] = pfqn_pff_delay(Z, n)
        for (m in 1..<M + 1) {
            G[m, idxn] = G[m - 1, idxn]
            for (r in 0..<R) {
                if (n[r] >= 1) {
                    n[r] = n[r] - 1
                    val idxn_1r = PopulationLattice.hashpop(n, N)
                    n[r] = n[r] + 1
                    val tmp_res = G[m, idxn] + L[m - 1, r] * G[m, idxn_1r]
                    G[m, idxn] = tmp_res
                }
            }
        }
        n = PopulationLattice.pprod(n, N)
    }

    val Gn = G[M, G.numCols - 1]
    val lGn = FastMath.log(Gn)
    return Ret.pfqnNc(Gn, lGn)
}
/**
 * PFQN ca algorithms
 */
@Suppress("unused")
class PfqnCaAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}