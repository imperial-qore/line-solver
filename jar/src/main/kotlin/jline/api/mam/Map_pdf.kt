/**
 * @file Markovian Arrival Process probability density function computation
 * 
 * Computes PDF values for MAP inter-arrival times using matrix exponential techniques.
 * Essential for probability analysis and statistical modeling of arrival processes.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.Maths
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the probability density function (PDF) of a MAP at specified time points.
 *
 * This function evaluates the PDF of the inter-arrival times of a Markovian Arrival Process
 * at the given time points. The PDF is computed using the formula:
 * f(t) = π * exp(D0 * t) * (-D0) * e
 *
 * where π is the steady-state probability of the embedded DTMC and e is a column vector of ones.
 *
 * @param MAP The Markovian Arrival Process stored in a MatrixCell, containing the (D0, D1) matrices
 * @param tset Array of time points at which to evaluate the PDF
 * @return Array of PDF values corresponding to the time points
 */
fun map_pdf(MAP: MatrixCell, tset: DoubleArray): DoubleArray {
    return map_pdf(MAP[0], MAP[1], tset)
}

/**
 * Computes the probability density function (PDF) of a MAP at specified time points.
 *
 * @param D0 The hidden transition matrix of the MAP
 * @param D1 The visible transition matrix of the MAP
 * @param tset Array of time points at which to evaluate the PDF
 * @return Array of PDF values corresponding to the time points
 */
fun map_pdf(D0: Matrix, D1: Matrix, tset: DoubleArray): DoubleArray {
    val pi = map_pie(D0, D1)
    val e = Matrix.ones(D1.numRows, 1)
    val minusD0 = D0.scale(-1.0)
    
    val result = DoubleArray(tset.size)
    
    for (i in tset.indices) {
        val t = tset[i]
        if (t < 0) {
            result[i] = 0.0
        } else if (t == 0.0) {
            // At t=0, the PDF is typically 0 for continuous distributions
            result[i] = 0.0
        } else {
            // Compute π * exp(D0 * t) * (-D0) * e
            val D0t = D0.scale(t)
            val expD0t = Maths.matrixExp(D0t)
            val temp = pi.mult(expD0t).mult(minusD0).mult(e)
            result[i] = temp.get(0, 0)
        }
    }
    
    return result
}

/**
 * Computes the probability density function (PDF) of a MAP at a single time point.
 *
 * @param MAP The Markovian Arrival Process stored in a MatrixCell
 * @param t Time point at which to evaluate the PDF
 * @return The PDF value at time t
 */
fun map_pdf(MAP: MatrixCell, t: Double): Double {
    return map_pdf(MAP, doubleArrayOf(t))[0]
}

/**
 * Computes the probability density function (PDF) of a MAP at a single time point.
 *
 * @param D0 The hidden transition matrix of the MAP
 * @param D1 The visible transition matrix of the MAP
 * @param t Time point at which to evaluate the PDF
 * @return The PDF value at time t
 */
fun map_pdf(D0: Matrix, D1: Matrix, t: Double): Double {
    return map_pdf(D0, D1, doubleArrayOf(t))[0]
}
/**
 * MAP pdf algorithms
 */
@Suppress("unused")
class MapPdfAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}