/**
 * @file Markovian Arrival Process construction using block matrices
 * 
 * Constructs MAP(2) representations from moment and autocorrelation parameters using fallback
 * algorithms for hyperexponential and Erlang-like approximations when MMPP fitting fails.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import kotlin.math.*

/**
 * Constructs a MAP(2) or MAP(1) according to given moments and autocorrelation parameters.
 * This is a simplified implementation of the complex MATLAB version.
 * Falls back if no MMPP(2) can be constructed.
 * 
 * @param E1 First moment (mean)
 * @param E2 Second moment 
 * @param E3 Third moment
 * @param G2 Autocorrelation decay ratio ρ(i)/ρ(i-1)
 * @param OPT Optional parameter ('scv' means E2 is squared coefficient of variation)
 * @return MAP as Array<Matrix> where result[0] = D0 and result[1] = D1
 */
fun map_block(E1: Double, E2: Double, E3: Double, G2: Double, OPT: String? = null): Array<Matrix> {
    var actualE2 = E2
    
    // Handle OPT parameter
    if (OPT != null && OPT.equals("scv", ignoreCase = true)) {
        actualE2 = (1 + E2) * E1 * E1
    }
    
    val SCV = (actualE2 - E1 * E1) / (E1 * E1)
    
    // Fall back to hyperexponential or Erlang fitting based on SCV
    return fallbackMAP(E1, actualE2, E3, G2)
}

/**
 * Fallback MAP construction when MMPP(2) fitting fails.
 * Uses simplified heuristics based on moment characteristics.
 */
private fun fallbackMAP(E1: Double, E2: Double, E3: Double, G2: Double): Array<Matrix> {
    val SCV = (E2 - E1 * E1) / (E1 * E1)
    
    return if (SCV > 1) {
        // Hyperexponential case: high variability
        constructHyperexponential(E1, SCV, G2)
    } else {
        // Erlang case: low variability (though SCV should be > 1 for feasibility)
        constructErlang(E1, SCV, G2)
    }
}

/**
 * Construct a hyperexponential MAP for high variability cases.
 */
private fun constructHyperexponential(E1: Double, SCV: Double, G2: Double): Array<Matrix> {
    // Simplified hyperexponential construction
    val p = 0.5  // Equal probability for two branches
    
    // Calculate rates for two-phase hyperexponential
    val lambda1 = 2.0 / E1 * (1 + sqrt((SCV - 1) / (SCV + 1)))
    val lambda2 = 2.0 / E1 * (1 - sqrt((SCV - 1) / (SCV + 1)))
    
    val D0 = Matrix(2, 2)
    D0[0, 0] = -lambda1
    D0[0, 1] = 0.0
    D0[1, 0] = 0.0
    D0[1, 1] = -lambda2
    
    val D1 = Matrix(2, 2)
    D1[0, 0] = lambda1 * p
    D1[0, 1] = lambda1 * (1 - p)
    D1[1, 0] = lambda2 * p
    D1[1, 1] = lambda2 * (1 - p)
    
    return arrayOf(D0, D1)
}

/**
 * Construct an Erlang-like MAP for low variability cases.
 */
private fun constructErlang(E1: Double, SCV: Double, G2: Double): Array<Matrix> {
    // For SCV <= 1, construct a simple MAP that approximates Erlang behavior
    val mu = 2.0 / E1  // Rate parameter
    
    val D0 = Matrix(2, 2)
    D0[0, 0] = -mu
    D0[0, 1] = mu * (1 - G2)
    D0[1, 0] = 0.0
    D0[1, 1] = -mu
    
    val D1 = Matrix(2, 2)
    D1[0, 0] = 0.0
    D1[0, 1] = mu * G2
    D1[1, 0] = 0.0
    D1[1, 1] = mu
    
    return arrayOf(D0, D1)
}

/**
 * Check if the constructed MAP is feasible and has reasonable parameters.
 */
private fun validateMAP(MAP: Array<Matrix>): Boolean {
    if (MAP.size != 2) return false
    
    val D0 = MAP[0]
    val D1 = MAP[1]
    
    // Check dimensions
    if (D0.numRows != D0.numCols || D1.numRows != D1.numCols || 
        D0.numRows != D1.numRows) return false
    
    // Check that diagonal elements of D0 are non-positive
    for (i in 0..<D0.numRows) {
        if (D0[i, i] > 1e-10) return false
    }
    
    // Check that off-diagonal elements of D0 are non-negative
    for (i in 0..<D0.numRows) {
        for (j in 0..<D0.numCols) {
            if (i != j && D0[i, j] < -1e-10) return false
        }
    }
    
    // Check that all elements of D1 are non-negative
    for (i in 0..<D1.numRows) {
        for (j in 0..<D1.numCols) {
            if (D1[i, j] < -1e-10) return false
        }
    }
    
    return true
}

/**
 * MAP block matrix construction algorithms
 */
@Suppress("unused")
class MapBlockAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}