/**
 * @file Cache Xi Terms via Fixed Point Iteration
 * 
 * Estimates cache xi terms using fixed-point algorithms. Xi terms represent
 * throughput duals in queueing network analysis adapted for cache systems,
 * essential for performance evaluation of hierarchical cache architectures.
 * 
 * @since LINE 3.0
 */
package jline.api.cache

import jline.io.Ret
import jline.util.Maths
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Estimate cache xi terms using a fixed-point algorithm. The xi terms represent duals of the throughput of a queueing
 * network but in the caching setting.
 *
 * @param gamma - Matrix representing the cache access factors.
 * @param m     - Matrix representing the cache capacity vector.
 * @param xi    - Optional matrix representing the initial value of the fixed-point iteration.
 * @return cacheXiFpReturn - An object containing the calculated xi terms, initial state probabilities (pi0),
 * the cache state probabilities (pij), and the number of iterations (it).
 */

fun cache_xi_fp(gamma: Matrix, m: Matrix, xi: Matrix?): Ret.cacheXiFp {
    var xi = xi
    val n = gamma.numRows
    val h = gamma.numCols
    val tol = 1.0e-14
    val pi0 = Matrix(1, n)
    pi0.fill(1.0 / (h + 1))
    val pij = Matrix(n, h)
    if (xi == null || xi.isEmpty) {
        xi = Matrix(1, h)
        for (l in 0..<h) {
            val mean = Matrix.extractColumn(gamma, l, null).elementSum() / gamma.numRows
            xi[l] = m[l] / mean / (n + m.elementSum() - 1)
        }
    }
    var it = 1
    while (it < 1.0e4) {
        val pi0_1 = Matrix(pi0)
        var mul = pi0_1.mult(gamma, null)
        xi = Matrix(m.numRows, m.numCols)
        for (i in 0..<xi.numRows) {
            for (j in 0..<xi.numCols) {
                xi[i, j] = m[i, j] / mul[i, j]
            }
        }
        val intermediate = gamma.elementMult(xi.repmat(n, 1), null)
        mul = gamma.mult(xi.repmat(n, 1).transpose())
        for (i in 0..<pij.numRows) {
            for (j in 0..<pij.numCols) {
                pij[i, j] = FastMath.abs(intermediate[i, j]) / FastMath.abs(1 + mul[i, j])
            }
        }
        for (i in 0..<n) {
            pi0[0, i] = Maths.max(tol, 1 - pij.sumRows(i))
        }
        var DELTA = 0.0
        for (i in 0..<n) {
            DELTA += FastMath.abs(1 - pi0[i] / pi0_1[i])
        }
        if (DELTA < tol) {
            break
        }
        it++
    }
    for (i in 0..<xi!!.numRows) {
        for (j in 0..<xi.numCols) {
            if (xi[i, j] < 0) {
                xi[i, j] = tol
            }
        }
    }
    return Ret.cacheXiFp(xi, pi0, pij, it)
}
/**
 * Cache xi fp algorithms
 */
@Suppress("unused")
class CacheXiFpAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}