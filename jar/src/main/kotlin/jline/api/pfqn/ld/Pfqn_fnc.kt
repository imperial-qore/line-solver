/**
 * @file Functional server scaling factor computation for load-dependent systems
 * 
 * Computes scaling factors for load-dependent functional servers in product-form queueing networks.
 * Handles the mathematical transformation of load-dependent service rates into functional scaling
 * parameters, supporting both automatic parameter selection and user-specified scaling constants.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.ld

import jline.GlobalConstants.Inf
import jline.io.Ret
import jline.util.Utils
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Compute scaling factor of a load-dependent functional server use to calculate the mean
 *
 * @param alpha - load-dependent scalings
 * @return functional server scalings
 */

fun pfqn_fnc(alpha: Matrix): Ret.pfqnFnc {
    val M = alpha.numRows
    var c = Matrix(1, M)
    c.zero()
    var mu = pfqn_fnc(alpha, c).mu
    if (!mu.isFinite) {
        c = Matrix.ones(1, M)
        c.scaleEq(-0.5)
        mu = pfqn_fnc(alpha, c).mu
    }
    var dt = 0.0
    while (!mu.isFinite) {
        dt += 0.05
        val c_scalar = -0.5 + dt
        c = Matrix.singleton(c_scalar)
        mu = pfqn_fnc(alpha, c).mu
        if (c_scalar >= 2) {
            break
        }
    }
    return Ret.pfqnFnc(mu, c)
}

/**
 * Compute scaling factor of a load-dependent functional server use to calculate the mean instantiated
 * with scaling constant c.
 *
 * @param alpha - load-dependent scalings
 * @return functional server scalings
 */

fun pfqn_fnc(alpha: Matrix, c: Matrix): Ret.pfqnFnc {
    val M = alpha.numRows
    val N = alpha.numCols
    val mu = Matrix(M, N)
    mu.zero()
    for (i in 0..<M) {
        mu[i, 0] = alpha[i, 0] / (1 + c[i])
        val alphanum = Matrix(N, N)
        alphanum.zero()
        val alphaden = alphanum.copy()
        for (n in 1..<N) {
            alphanum[n, 0] = alpha[i, n]
            alphaden[n, 0] = alpha[i, n - 1]
            for (k in 1..<n) {
                alphanum[n, k] = alphanum[n, k - 1] * alpha[i, n - k]
                alphaden[n, k] = alphaden[n, k - 1] * alpha[i, n - k - 1]
            }
        }
        for (n in 1..<N) {
            var rho = 0.0
            var muden = 1.0
            for (k in 0..<n) {
                muden *= mu[i, k]
                rho += (alphanum[n, k] - alphaden[n, k]) / muden
            }
            mu[i, n] = (alphanum[n, n - 1] * alpha[i, 0] / muden) / (1 - rho)
        }
    }
    for (i in 0..<M) {
        for (j in 0..<N) {
            if (java.lang.Double.isNaN(mu[i, j]) || FastMath.abs(mu[i, j]) > 1e15) {
                mu[i, j] = Inf
            }
        }
    }
    for (i in 0..<M) {
        if (Matrix.extractRows(mu, i, i + 1, null).isFinite) {
            continue
        }
        var replaceWithInf = false
        for (j in 0..<N) {
            if (replaceWithInf) {
                mu[i, j] = Inf
            } else if (Utils.isInf(mu[i, j])) {
                replaceWithInf = true
            }
        }
    }
    return Ret.pfqnFnc(mu, c)
}
/**
 * PFQN fnc algorithms
 */
@Suppress("unused")
class PfqnFncAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}