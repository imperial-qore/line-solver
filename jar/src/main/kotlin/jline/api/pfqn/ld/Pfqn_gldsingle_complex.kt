/**
 * @file Complex-valued single-class load-dependent auxiliary normalizing constant computation
 * 
 * Provides specialized auxiliary function for computing normalizing constants in single-class
 * load-dependent closed queueing networks with complex-valued service demands. Implements
 * efficient recursive computation optimized for complex parameter systems.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.ld

import jline.io.Ret
import jline.solvers.SolverOptions
import jline.util.matrix.ComplexMatrix
import jline.util.matrix.Matrix
import org.apache.commons.math3.complex.Complex

/**
 * Auxiliary function used by pfqn_gld to computer the normalizing constant in a single-class load-dependent model
 * with complex demands.
 *
 * @param L       - demands at all stations
 * @param N       - number of jobs for each class
 * @param mu      - load-dependent scaling factors
 * @param options - solver options
 * @return normalizing constant (G) and its logarithm (lG)
 */

fun pfqn_gldsingle_complex(L: ComplexMatrix, N: Matrix, mu: Matrix, options: SolverOptions?): Ret.pfqnNcComplex {
    val M = L.numRows
    val R = L.numCols

    if (R > 1) {
        throw RuntimeException("pfqn_gldsingle_complex: multiclass model detected. pfqn_gldsingle_complex is for single class models.")
    }
    val g: MutableMap<Ret.pfqnGldIndex, Complex> = HashMap()
    g[Ret.pfqnGldIndex(1, 1, 1)] = Complex.valueOf(0.0)
    var n = 1
    while (n <= N[0]) {
        g[Ret.pfqnGldIndex(1, n + 1, 2)] = Complex.valueOf(0.0)
        n++
    }
    for (m in 1..M) {
        var tm = 1
        while (tm <= N[0] + 1) {
            g[Ret.pfqnGldIndex(m + 1, 1, tm + 1)] = Complex.valueOf(1.0)
            tm++
        }
        var n = 1
        while (n <= N[0]) {
            var tm = 1
            while (tm <= N[0] - n + 1) {
                g[Ret.pfqnGldIndex(m + 1, n + 1, tm + 1)] =
                    g[Ret.pfqnGldIndex(m, n + 1, 2)]!!.add(L[m - 1].multiply(g[Ret.pfqnGldIndex(m + 1, n, tm + 2)])
                        .divide(mu[m - 1, tm - 1]))
                tm++
            }
            n++
        }
    }
    val G = g[Ret.pfqnGldIndex(M + 1, N[0].toInt() + 1, 2)]
    val lG = G!!.log()
    return Ret.pfqnNcComplex(G, lG)
}
/**
 * PFQN gldsingle complex algorithms
 */
@Suppress("unused")
class PfqnGldsingleComplexAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}