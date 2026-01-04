/**
 * @file Grundmann-Möller cubature method for normalizing constant computation
 * 
 * Implements the cubature (multi-dimensional integration) approach for computing normalizing
 * constants in closed product-form queueing networks. Uses Grundmann-Möller simplex quadrature
 * with configurable order for exact or approximate computation.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.io.Ret
import jline.util.Maths
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Cubature method to compute the normalizing constant of a load-independent closed queueing network model
 *
 * @param L - demands at all stations
 * @param N - number of jobs for each class
 * @param Z - think time for each class
 * @return normalizing constant, its logarithm, and mean performance metrics computed as a by-product
 */

fun pfqn_cub(L: Matrix, N: Matrix, Z: Matrix): Ret.pfqnNc {
    val order = FastMath.ceil((N.elementSum() - 1) / 2.0).toInt()
    val atol = 1e-8
    return pfqn_cub(L, N, Z, order, atol)
}

/**
 * Cubature method to compute the normalizing constant of a load-independent closed queueing network model
 *
 * @param L     - demands at all stations
 * @param N     - number of jobs for each class
 * @param Z     - think time for each class
 * @param order - cubature order ( ceil((|N|-1)/2) or above is exact)
 * @param atol  - numerical tolerance for simplex quadrature
 * @return normalizing constant, its logarithm, and mean performance metrics computed as a by-product
 */

fun pfqn_cub(L: Matrix, N: Matrix, Z: Matrix, order: Int, atol: Double): Ret.pfqnNc {
    var order = order
    var atol = atol
    val M = L.numRows
    L.numCols

    if (L.isEmpty || N.isEmpty || N.elementSum() == 0.0) {
        return Ret.pfqnNc(1.0, 0.0)
    }

    if (order < 0) {
        order = FastMath.ceil((N.elementSum() - 1) / 2.0).toInt()
    }

    if (atol <= 0) {
        atol = 1e-8
    }

    if (Z.isEmpty || Z.elementSum() < atol) {
        val Nt = N.elementSum()
        val beta = N.scale(1.0 / Nt)

        val simplexResult = Maths.simplexquad(Ret.pfqnCUB(L, Nt, beta), M - 1, order, atol)
        val Q = simplexResult.Q
        val Gn = Q[Q.size - 1] * FastMath.exp(Maths.factln(N.elementSum() + M - 1) - N.factln().elementSum())
        return Ret.pfqnNc(Gn, FastMath.log(Gn))
    } else {
        val steps = 10000
        val Nt = N.elementSum()
        val beta = N.scale(1 / Nt)
        var Gn = 0.0
        val vmax = Nt * 10
        val dv = vmax / steps

        var v = 0.0
        while (v <= vmax) {
            val Lv = L.scale(v).add(Z.repmat(M, 1))
            val simplexResult = Maths.simplexquad(Ret.pfqnCUB(Lv, Nt, beta), M - 1, order, atol)
            val Q = simplexResult.Q
            val dG = FastMath.exp(-v) * FastMath.pow(v, M - 1) * Q[Q.size - 1] * dv
            Gn += dG

            if (v > 0 && dG / Gn < atol) {
                break
            }
            v += dv
        }

        Gn *= FastMath.exp(-N.factln().elementSum())
        return Ret.pfqnNc(Gn, FastMath.log(Gn))
    }
}
/**
 * PFQN cub algorithms
 */
@Suppress("unused")
class PfqnCubAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}