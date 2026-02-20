/**
 * @file Markovian Arrival Process autocorrelation decay rate analysis
 * 
 * Computes gamma parameter measuring autocorrelation decay rates in MAP processes.
 * Essential for characterizing temporal correlation patterns and burstiness in arrival streams.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.FastMath

/**
 * Computes the gamma parameter for a MAP, which is the autocorrelation decay rate.
 *
 * For MAPs of order 1 (Poisson process), gamma is 0.
 * For MAPs of order 2, gamma is the ratio of ACF at lag 2 to ACF at lag 1.
 * For higher-order MAPs, this implementation uses a simplified approximation.
 *
 * @param D0 the hidden transition matrix of the MAP
 * @param D1 the visible transition matrix of the MAP
 * @return the gamma parameter (autocorrelation decay rate)
 */
fun map_gamma(D0: Matrix, D1: Matrix): Double {
    val n = D0.numRows
    
    return when {
        n == 1 -> {
            // Poisson process: no correlation
            0.0
        }
        n == 2 -> {
            // Second-order MAP: geometric ACF
            val acf1 = map_acf(D0, D1, 1).value()
            if (FastMath.abs(acf1) < 1e-8) {
                // Phase-type distribution
                0.0
            } else {
                map_acf(D0, D1, 2).value() / acf1
            }
        }
        else -> {
            // Higher-order MAP: simplified approximation
            // Use the basic ratio approach as an approximation
            val acf1 = map_acf(D0, D1, 1).value()
            if (FastMath.abs(acf1) < 1e-8) {
                0.0
            } else {
                val acf2 = map_acf(D0, D1, 2).value()
                val gamma = acf2 / acf1
                // Bound gamma to reasonable range
                FastMath.max(-0.999, FastMath.min(0.999, gamma))
            }
        }
    }
}

/**
 * Computes the gamma parameter for a MAP using a MatrixCell.
 *
 *
 * This method is a convenience overload that extracts the D0 and D1 matrices from the provided MatrixCell
 * and calculates the gamma parameter.
 *
 * @param MAP the MatrixCell representing the MAP, containing the D0 and D1 matrices
 * @return the gamma parameter
 */

fun map_gamma(MAP: MatrixCell?): Double {
    return if (MAP != null) {
        map_gamma(MAP[0], MAP[1])
    } else {
        0.0
    }
}
/**
 * MAP gamma algorithms
 */
@Suppress("unused")
class MapGammaAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}