/**
 * @file Multi-server load-dependent scaling factor computation
 * 
 * Computes load-dependent scaling factors for multi-server queueing stations with finite
 * server capacity. Implements recursive algorithms for computing state-dependent service rates
 * in multi-server environments with load-dependent service mechanisms.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.ld

import jline.util.Maths
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

fun pfqn_mu_ms(N: Int, m: Int, c: Int): Matrix {
    val mu = Matrix.zeros(1, N)
    val g = Matrix.zeros(m, N + 1) // table

    for (n in 0..N) {
        for (i in 0..<m) {
            g[i, n] = pfqn_mu_ms_gnaux(n, i, c, g)
        }
    }

    for (n in 1..N) {
        mu[0, n - 1] = g[m - 1, n - 1] / g[m - 1, n]
    }

    return mu
}

fun pfqn_mu_ms_gnaux(n: Int, m: Int, c: Int, g: Matrix): Double {
    if (n == 0) {
        return 1.0
    } else {
        if (m == 0) {
            var prodmin = 1.0
            for (i in 1..n) {
                prodmin *= FastMath.min(i, c).toDouble()
            }
            return 1.0 / prodmin
        } else {
            var gn = 0.0
            for (k in 0..n) {
                val a: Double
                if (k > 0) {
                    val avec = Matrix(1, k)
                    for (t in 0..<k) {
                        avec[0, t] = Maths.min((1 + t).toDouble(), c.toDouble())
                    }
                    a = avec.elementMult()
                } else {
                    a = 1.0
                }
                val b = 1 / g[m - 1, n - k]
                gn += 1.0 / (a * b)
            }
            return gn
        }
    }
}
/**
 * PFQN mu ms algorithms
 */
@Suppress("unused")
class PfqnMuMsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}