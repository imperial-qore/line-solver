/**
 * @file Single-class load-dependent normalizing constant auxiliary computation
 * 
 * Provides specialized auxiliary function for computing normalizing constants in single-class
 * load-dependent closed queueing networks. Implements efficient recursive computation optimized
 * for single-class systems with state-dependent service rates.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.ld

import jline.io.Ret
import jline.solvers.SolverOptions
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Auxiliary function used by pfqn_gld to computer the normalizing constant in a single-class load-dependent model.
 *
 * @param L       - demands at all stations
 * @param N       - number of jobs for each class
 * @param mu      - load-dependent scaling factors
 * @param options - solver options
 * @return normalizing constant (G) and its logarithm (lG)
 */

fun pfqn_gldsingle(L: Matrix, N: Matrix, mu: Matrix, options: SolverOptions?): Ret.pfqnNc {
    val M = L.numRows
    val R = L.numCols

    if (R > 1) {
        throw RuntimeException("pfqn_gldsingle: multiclass model detected. pfqn_gldsingle is for single class models.")
    }

    val g: MutableMap<Ret.pfqnGldIndex, Double> = HashMap()
    
    // Initialize boundary conditions: g(0+1, n+1, 1+1) = 0 for n=1:N
    var n = 1
    while (n <= N[0]) {
        g[Ret.pfqnGldIndex(0 + 1, n + 1, 1 + 1)] = 0.0
        n++
    }

    for (m in 1..M) {
        // Initialize boundary conditions: g(m+1, 0+1, tm+1) = 1 for tm=1:(N+1)
        var tm = 1
        while (tm <= N[0] + 1) {
            g[Ret.pfqnGldIndex(m + 1, 0 + 1, tm + 1)] = 1.0
            tm++
        }
        var n = 1
        while (n <= N[0]) {
            var tm = 1
            while (tm <= N[0] - n + 1) {
                val gPrev = g[Ret.pfqnGldIndex(m - 1 + 1, n + 1, 1 + 1)] ?: 0.0
                val gCurr = g[Ret.pfqnGldIndex(m + 1, n - 1 + 1, tm + 1 + 1)] ?: 0.0
                g[Ret.pfqnGldIndex(m + 1, n + 1, tm + 1)] = 
                    gPrev + L[m - 1] * gCurr / mu[m - 1, tm - 1]
                tm++
            }
            n++
        }
    }
    val G = g[Ret.pfqnGldIndex(M + 1, N[0].toInt() + 1, 1 + 1)]!!
    val lG = FastMath.log(G)
    return Ret.pfqnNc(G, lG)
}
/**
 * PFQN gldsingle algorithms
 */
@Suppress("unused")
class PfqnGldsingleAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}