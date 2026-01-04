/**
 * CTMC Transient Analysis via Uniformization
 * 
 * Computes the transient probability distribution of CTMCs using the uniformization
 * method, which transforms the continuous-time problem into a weighted sum of DTMC
 * powers. Provides numerically stable computation of time-dependent state probabilities.
 *
 * @since LINE 3.0
 */
package jline.api.mc

import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import kotlin.math.exp

/**
 * Return the transient probability distribution of the CTMC via the
 * uniformization method.
 *
 * @param pi0 Initial state the CTMC
 * @param Q   Infinitesimal generator of the CTMC
 * @param t   Transient analysis period boundary [0,t]
 * @return Transient probability vector at time t
 */

fun ctmc_uniformization(pi0: Matrix, Q: Matrix, t: Double): Matrix {
    val tol = 1e-12
    val maxiter = 100
    var q = 0.0
    val n = Q.numCols
    for (i in 0..<n) {
        q = FastMath.max(q, 1.1 * FastMath.abs(Q[i, i]))
    }
    var Qs = Matrix.eye(n)
    Qs = Qs.add(1.0 / q, Q)
    var k = 0
    var s = 1.0
    var r = 1.0
    var iter = 0
    var kmax = 1
    while (iter < maxiter) {
        iter++
        k++
        r = r * (q * t) / k
        s = s + r
        if (1 - FastMath.exp(-q * t) * s <= tol) {
            kmax = k
            break
        }
    }

    var pi = Matrix(1, n)
    Matrix(1, n)
    Matrix(1, n)
    pi0.scaleEq(exp(-q * t), pi)
    val P = Matrix(pi0)
    var ri = FastMath.exp(-q * t)
    for (j in 0..<kmax) {
        P.multEq(Qs)
        ri = ri * (q * t / (j + 1))
        pi = pi.add(ri, P)
    }
    return pi
}
/**
 * CTMC uniformization algorithms
 */
@Suppress("unused")
class CtmcUniformizationAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}