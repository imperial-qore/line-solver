/**
 * @file Absorbing Phase-type distribution two-phase fitting
 * 
 * Fits APH(2) distributions to match given moments with automatic feasibility adjustment.
 * Essential for modeling service time distributions with bounded support and low variability.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.io.Ret

/**
 * Fits an acyclic phase-type (APH) distribution with two phases (APH(2)) to match the given moments of a random variable.
 *
 *
 * This method attempts to fit an APH(2) distribution based on the first three moments (M1, M2, M3) of the distribution.
 * If a valid APH(2) distribution cannot be directly obtained, the method adjusts the moments and retries the fitting process.
 * If the fitting still fails, it throws an exception indicating that feasibility could not be restored.
 *
 * @param M1 the first moment (mean) of the distribution
 * @param M2 the second moment of the distribution
 * @param M3 the third moment of the distribution
 * @return an object of type `aph2_fit_return_type` containing the fitted APH(2) distribution and all possible fits
 * @throws RuntimeException if the APH(2) fitting is not feasible with the given or adjusted moments
 */
fun aph2_fit(M1: Double, M2: Double, M3: Double): Ret.mamAPH2Fit {
    val result = Ret.mamAPH2Fit()
    result.APHS = aph2_fitall(M1, M2, M3)

    if (result.APHS.isEmpty()) {
        result.APHS = aph2_fitall(M1, aph2_adjust(M1, M2, M3)[0]!!, aph2_adjust(M1, M2, M3)[1]!!)
        if (result.APHS.isEmpty()) {
            throw RuntimeException("Fitting APH(2): feasibility could not be restored")
        }
    }

    result.APH = result.APHS[0]

    return result
}
/**
 * APH 2 fit algorithms
 */
@Suppress("unused")
class Aph2FitAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}