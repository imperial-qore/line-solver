/**
 * @file Markovian Arrival Process feasible block matrix construction
 * 
 * Constructs feasible MAP representations when exact moment matching fails by adjusting
 * parameters and ensuring mathematical constraints are satisfied. Used for robust MAP fitting.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.io.line_warning
import jline.util.matrix.Matrix

/**
 * Fits the most similar feasible MAP when exact moment matching fails.
 * Ensures feasibility constraints are met.
 * 
 * @param E1 First moment (mean)
 * @param E2 Second moment 
 * @param E3 Third moment
 * @param G2 Autocorrelation decay ratio ρ(i)/ρ(i-1)
 * @param OPT Optional parameter ('scv' means E2 is squared coefficient of variation)
 * @return MAP as Array<Matrix> where result[0] = D0 and result[1] = D1
 */
fun map_feasblock(E1: Double, E2: Double, E3: Double, G2: Double, OPT: String? = null): Array<Matrix> {
    var actualE2 = E2
    var actualE3 = E3
    
    // Handle exponential case
    if (actualE2 == 2 * E1 * E1) {
        val D0 = Matrix(arrayOf(doubleArrayOf(-1.0, 0.0), doubleArrayOf(0.0, -1.0)))
        val D1 = Matrix(arrayOf(doubleArrayOf(0.5, 0.5), doubleArrayOf(0.5, 0.5)))
        val scaledMAP = map_scale(D0, D1, E1)
        return arrayOf(scaledMAP[0], scaledMAP[1])
    }
    
    // Handle OPT parameter
    if (OPT != null && OPT.equals("scv", ignoreCase = true)) {
        actualE2 = (1 + E2) * E1 * E1
    }
    
    val tolerance = 1e-10  // kpcfit_tol equivalent
    
    // Check feasibility constraints and adjust if necessary
    if (actualE2 <= 2 * E1 * E1) {
        line_warning("map_feasblock", "E2 failure (SCV≤1), setting SCV=1.001")
        actualE2 = (1 + tolerance) * E1 * E1
    }

    val minE3 = (3.0/2.0) * actualE2 * actualE2 / E1
    if (actualE3 <= minE3) {
        line_warning("map_feasblock", "E3 failure, setting E3=(3/2+1e-6)*E2^2/E1")
        actualE3 = (3.0/2.0 + tolerance) * actualE2 * actualE2 / E1
    }
    
    // Call map_block with adjusted parameters
    return map_block(E1, actualE2, actualE3, G2)
}
/**
 * MAP feasblock algorithms
 */
@Suppress("unused")
class MapFeasblockAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}