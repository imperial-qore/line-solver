/**
 * @file Markov Modulated Poisson Process single-parameter fitting
 * 
 * Fits MMPP(2) models using simplified single-parameter approach for specific scenarios.
 * Provides efficient fitting for cases with limited statistical information.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.FastMath

/**
 * Fits a 2-phase Markov Modulated Poisson Process (MMPP2) based on the specified parameters.
 *
 * @param mean the mean inter-arrival time
 * @param scv  the squared coefficient of variation (SCV)
 * @param skew the skewness of the inter-arrival times
 * @param idc  the index of dispersion for counts (IDC)
 * @return a MatrixCell representing the fitted MMPP2 transition matrices
 */
fun mmpp2_fit1(mean: Double, scv: Double, skew: Double, idc: Double): MatrixCell {
    val E1 = mean
    val E2 = (1 + scv) * FastMath.pow(E1, 2)
    val g2 = -(scv - idc) / (-1 + idc)
    val E3 = if (skew == -1.0) {
        -1.0
    } else {
        -(2 * FastMath.pow(E1, 3) - 3 * E1 * E2 - skew * FastMath.pow(E2 - FastMath.pow(E1, 2), 1.5))
    }

    return map2_fit(E1, E2, E3, g2).MAP
}
/**
 * Mmpp2 Fit1 algorithms
 */
@Suppress("unused")
class Mmpp2Fit1Algo {
    companion object {
        // Class documentation marker for Dokka
    }
}