/**
 * @file Markovian Arrival Process hyperexponential distribution fitting
 * 
 * Constructs MAP representations of two-phase hyperexponential renewal processes.
 * Used for modeling high-variability arrival processes with coefficient of variation greater than 1.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.FastMath

/**
 * Fit a two-phase Hyper-exponential renewal process as a MAP
 *
 * @param mean: mean inter-arrival time of the process
 * @param scv:  squared coefficient of variation of inter-arrival times
 * @param p:    p: probability of being served in phase 1 (DEFAULT: p=0.99)
 * @return Fitted hyper-exponential process
 */

fun map_hyperexp(mean: Double, scv: Double, p: Double): MatrixCell? {
    var p = p
    val D = MatrixCell()
    if (p == 0.0) {
        p = 0.99
    }

    val e2 = (1.0 + scv) * mean * mean
    val delta = -4.0 * p * mean * mean + 4.0 * p * p * mean * mean + 2.0 * e2 * p - 2.0 * e2 * p * p
    var mu2 = (-2.0 * mean + 2.0 * p * mean + FastMath.sqrt(delta)) / (e2 * p - 2.0 * mean * mean)
    var mu1 = mu2 * p / (p - 1.0 + mean * mu2)
    val D0 = Matrix(2, 2)
    val D1 = Matrix(2, 2)
    D0[0, 0] = -mu1
    D0[1, 1] = -mu2
    D1[0, 0] = mu1 * p
    D1[0, 1] = mu1 * (1.0 - p)
    D1[1, 0] = mu2 * p
    D1[1, 1] = mu2 * (1.0 - p)
    D[0] = D0
    D[1] = D1

    if (map_isfeasible(D)) {
        return D
    } else {
        mu2 = (-2 * mean + 2 * p * mean - FastMath.sqrt(delta)) / (e2 * p - 2 * mean * mean)
        mu1 = mu2 * p / (p - 1 + mean * mu2)
        D0.zero()
        D1.zero()
        D0[0, 0] = -mu1
        D0[1, 1] = -mu2
        D1[0, 0] = mu1 * p
        D1[0, 1] = mu1 * (1.0 - p)
        D1[1, 0] = mu2 * p
        D1[1, 1] = mu2 * (1.0 - p)
        D[0] = D0
        D[1] = D1
        return if (p > 1e-6 && !map_isfeasible(D)) {
            map_hyperexp(mean, scv, p / 10.0)
        } else {
            null
        }
    }
}
/**
 * MAP hyperexp algorithms
 */
@Suppress("unused")
class MapHyperexpAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}