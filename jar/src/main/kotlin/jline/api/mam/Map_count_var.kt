/**
 * @file Markovian Arrival Process counting process variance analysis
 * 
 * Computes variance of MAP counting processes over specified time intervals using matrix
 * exponential methods. Essential for analyzing variability in arrival patterns and burstiness.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.Maths
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.FastMath

/**
 * Computes the variance of the counting process over a specified interval length for a given Markovian Arrival Process (MAP).
 *
 *
 * This function calculates the variance in the number of events occurring in a MAP over an interval of length `t`.
 * It uses the mean counting process function for a single interval length and retrieves the first element.
 *
 *
 * @param MAP The Markovian Arrival Process stored in a MatrixCell.
 * @param t   The length of the interval over which to compute the variance.
 * @return The variance of the counting process over the interval `t`.
 */
fun map_count_var(MAP: MatrixCell, t: Double): Double {
    return map_count_mean(MAP, doubleArrayOf(t))[0]
}

/**
 * Computes the variance of the counting process over multiple specified interval lengths for a given Markovian Arrival Process (MAP).
 *
 *
 * This function calculates the variance in the number of events occurring in a MAP over each interval length provided in the array `t`.
 * It employs the MAP's infinitesimal generator, stationary distribution, and other parameters to determine the variance.
 * The results are returned in an array corresponding to the input intervals.
 *
 *
 * @param MAP The Markovian Arrival Process stored in a MatrixCell.
 * @param t   An array of interval lengths over which to compute the variance.
 * @return An array of doubles, where each element represents the variance of the counting process over the corresponding interval length in `t`.
 */
fun map_count_var(MAP: MatrixCell, t: DoubleArray): DoubleArray {
    val ret = DoubleArray(t.size)
    val l = map_lambda(MAP)
    val D = map_infgen(MAP)
    val theta = map_piq(MAP)
    val D1 = MAP[1]
    val n = MAP[0].length()
    val e = Matrix.ones(n, 1)
    val tmp = e.mult(theta).sub(D).inv()
    val c = theta.mult(D1).mult(tmp)
    val d = tmp.mult(D1).mult(e)
    val ll = theta.mult(D1).mult(e)
    val l2 = 2 * FastMath.pow(l, 2)
    for (i in t.indices) {
        val t_i = t[i]
        val expmDt = Maths.matrixExp(D.scale(t_i))
        val term1 = ll.sub(Matrix.createLike(ll).fill(l2)).add(c.mult(MAP[1]).mult(e)).scale(t_i)
        val term2 = c.mult(Matrix.eye(n).sub(expmDt)).mult(d).scale(2.0)
        ret[i] = term1.sub(term2).value()
    }
    return ret
}
/**
 * MAP count var algorithms
 */
@Suppress("unused")
class MapCountVar {
    companion object {
        // Class documentation marker for Dokka
    }
}