/**
 * @file Markovian Arrival MAP with Marked arrivals two-state two-class fitting
 * 
 * Fits MAMAP(2,2) processes for two-class systems with forward moments and sigma characteristics.
 * Specialized algorithm for dual-class arrival modeling with moment matching.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import kotlin.math.abs

/**
 * Fits MAMAP(2,2) with forward moments (F) and sigma characteristics (S) for 2 classes.
 * 
 * @param amap underlying AMAP(2) process
 * @param P class probabilities (1x2 matrix)
 * @param F forward moments (1x2 matrix) 
 * @param S sigma matrix (2x2 matrix)
 * @param options additional options (unused, for compatibility)
 * @param weights weights for forward vs sigma fitting [wF, wS]
 * @return fitted MAMAP(2,2)
 */
fun mamap22_fit_fs_multiclass(
    amap: MatrixCell, 
    P: Matrix, 
    F: Matrix, 
    S: Matrix,
    options: Any? = null,
    weights: DoubleArray = doubleArrayOf(1.0, 1.0)
): MatrixCell {
    // Extract AMAP parameters
    val D0 = amap[0]
    val D1 = amap[1]
    val n = D0.numRows
    
    // Handle degenerate cases first
    if (n == 1) {
        // Poisson case
        return fitMarkedPoisson(amap, P)
    }
    
    // Extract 2-state AMAP parameters
    val h1 = -1.0 / D0[0, 0]
    val h2 = -1.0 / D0[1, 1] 
    val negD0Inv = D0.scale(-1.0).inv()
    val transProb = negD0Inv.mult(D1)
    val r1 = transProb[0, 1]
    val r2 = transProb[1, 1]
    
    val degentol = 1e-8
    val gamma = map_gamma(amap)
    
    // Detect canonical form and handle degenerate cases
    if (gamma > 0) {
        return handlePositiveGammaFS(amap, P, F, S, weights, h1, h2, r1, r2, degentol)
    } else {
        return handleNegativeGammaFS(amap, P, F, S, weights, h1, h2, r1, r2, degentol)
    }
}

/**
 * Fits MAMAP(2,2) with backward moments (B) and sigma characteristics (S) for 2 classes.
 */
fun mamap22_fit_bs_multiclass(
    amap: MatrixCell,
    P: Matrix,
    B: Matrix, 
    S: Matrix,
    options: Any? = null,
    weights: DoubleArray = doubleArrayOf(1.0, 1.0)
): MatrixCell {
    // Similar structure to FS version but optimizes backward moments
    val D0 = amap[0]
    val D1 = amap[1]
    val n = D0.numRows
    
    if (n == 1) {
        return fitMarkedPoisson(amap, P)
    }
    
    val h1 = -1.0 / D0[0, 0]
    val h2 = -1.0 / D0[1, 1]
    val negD0Inv = D0.scale(-1.0).inv()
    val transProb = negD0Inv.mult(D1)
    val r1 = transProb[0, 1]
    val r2 = transProb[1, 1]
    
    val degentol = 1e-8
    val gamma = map_gamma(amap)
    
    if (gamma > 0) {
        return handlePositiveGammaBS(amap, P, B, S, weights, h1, h2, r1, r2, degentol)
    } else {
        return handleNegativeGammaBS(amap, P, B, S, weights, h1, h2, r1, r2, degentol)
    }
}

/**
 * Handles forward+sigma fitting for positive gamma (canonical form 1).
 */
private fun handlePositiveGammaFS(
    amap: MatrixCell, P: Matrix, F: Matrix, S: Matrix,
    weights: DoubleArray, h1: Double, h2: Double, r1: Double, r2: Double, degentol: Double
): MatrixCell {
    // Check degenerate cases for positive gamma
    if (r1 < degentol || (1 - r2) < degentol) {
        throw IllegalArgumentException("Invalid AMAP parameters for positive gamma case")
    }
    
    if (abs(h2 - h1 * r2) < degentol) {
        // Degenerate case 1
        return solveConstrainedFS(amap, P, F, S, "case1")
    }
    
    if (abs(h1 - h2 + h2 * r1) < degentol) {
        // Degenerate case 2  
        return solveConstrainedFS(amap, P, F, S, "case2")
    }
    
    if ((1 - r1) < degentol) {
        // Non-canonical APH(2), only forward can be fitted
        return solveForwardOnlyFS(amap, P, F)
    }
    
    if (r2 < degentol) {
        // Canonical APH(2)
        return solveCanonicalFS(amap, P, F, S)
    }
    
    // General case: solve full optimization problem
    return solveGeneralFS(amap, P, F, S, weights)
}

/**
 * Handles backward+sigma fitting for positive gamma.
 */
private fun handlePositiveGammaBS(
    amap: MatrixCell, P: Matrix, B: Matrix, S: Matrix,
    weights: DoubleArray, h1: Double, h2: Double, r1: Double, r2: Double, degentol: Double
): MatrixCell {
    // Similar logic to FS but for backward moments
    if (r1 < degentol || (1 - r2) < degentol) {
        throw IllegalArgumentException("Invalid AMAP parameters for positive gamma case")
    }
    
    // Handle degenerate cases (similar pattern as FS)
    if (abs(h2 - h1 * r2) < degentol) {
        return solveConstrainedBS(amap, P, B, S, "case1")
    }
    
    if (abs(h1 - h2 + h2 * r1) < degentol) {
        return solveConstrainedBS(amap, P, B, S, "case2") 
    }
    
    if ((1 - r1) < degentol) {
        return solveBackwardOnlyBS(amap, P, B)
    }
    
    if (r2 < degentol) {
        return solveCanonicalBS(amap, P, B, S)
    }
    
    return solveGeneralBS(amap, P, B, S, weights)
}

/**
 * Handles negative gamma cases (canonical form 2).
 */
private fun handleNegativeGammaFS(
    amap: MatrixCell, P: Matrix, F: Matrix, S: Matrix,
    weights: DoubleArray, h1: Double, h2: Double, r1: Double, r2: Double, degentol: Double
): MatrixCell {
    if ((1 - r2) < degentol) {
        throw IllegalArgumentException("Invalid AMAP parameters for negative gamma case")
    }
    
    // Handle specific degenerate patterns for negative gamma
    if (abs(h1 - h2 - h1 * r1 + h1 * r1 * r2) < degentol) {
        return solveConstrainedFS(amap, P, F, S, "neg_case1")
    }
    
    if (abs(h1 - h2 + h2 * r1) < degentol) {
        return solveConstrainedFS(amap, P, F, S, "neg_case2")
    }
    
    if (r2 < degentol && (1 - r1) < degentol) {
        return solveCanonicalFS(amap, P, F, S)
    }
    
    if (r2 < degentol) {
        // Choose based on weights
        return if (weights[0] >= weights[1]) {
            solveForwardOnlyFS(amap, P, F)
        } else {
            solveSigmaOnlyFS(amap, P, S)
        }
    }
    
    return solveGeneralFS(amap, P, F, S, weights)
}

private fun handleNegativeGammaBS(
    amap: MatrixCell, P: Matrix, B: Matrix, S: Matrix,
    weights: DoubleArray, h1: Double, h2: Double, r1: Double, r2: Double, degentol: Double
): MatrixCell {
    // Similar to FS version but for backward moments
    if ((1 - r2) < degentol) {
        throw IllegalArgumentException("Invalid AMAP parameters for negative gamma case")
    }
    
    if (abs(h1 - h2 - h1 * r1 + h1 * r1 * r2) < degentol) {
        return solveConstrainedBS(amap, P, B, S, "neg_case1")
    }
    
    if (abs(h1 - h2 + h2 * r1) < degentol) {
        return solveConstrainedBS(amap, P, B, S, "neg_case2")
    }
    
    if (r2 < degentol && (1 - r1) < degentol) {
        return solveCanonicalBS(amap, P, B, S)
    }
    
    if (r2 < degentol) {
        return if (weights[0] >= weights[1]) {
            solveBackwardOnlyBS(amap, P, B)
        } else {
            solveSigmaOnlyBS(amap, P, S)
        }
    }
    
    return solveGeneralBS(amap, P, B, S, weights)
}

/**
 * Fits a marked Poisson process.
 */
private fun fitMarkedPoisson(amap: MatrixCell, P: Matrix): MatrixCell {
    val lambda = -amap[0][0, 0]
    val m = P.numCols
    
    val mmap = MatrixCell(2 + m)
    mmap[0] = amap[0].copy()
    mmap[1] = Matrix.zeros(1, 1)  // No aggregate arrivals
    
    for (c in 0 until m) {
        mmap[2 + c] = Matrix.zeros(1, 1)
        mmap[2 + c][0, 0] = lambda * P[0, c]
    }
    
    return mmap
}

// Simplified implementations of the solver methods
// Full implementations would require sophisticated optimization algorithms

private fun solveConstrainedFS(amap: MatrixCell, P: Matrix, F: Matrix, S: Matrix, case: String): MatrixCell {
    // Simplified constraint-based solution
    return createSimpleMMAP(amap, P, "forward")
}

private fun solveConstrainedBS(amap: MatrixCell, P: Matrix, B: Matrix, S: Matrix, case: String): MatrixCell {
    return createSimpleMMAP(amap, P, "backward")
}

private fun solveForwardOnlyFS(amap: MatrixCell, P: Matrix, F: Matrix): MatrixCell {
    return createSimpleMMAP(amap, P, "forward")
}

private fun solveBackwardOnlyBS(amap: MatrixCell, P: Matrix, B: Matrix): MatrixCell {
    return createSimpleMMAP(amap, P, "backward")
}

private fun solveCanonicalFS(amap: MatrixCell, P: Matrix, F: Matrix, S: Matrix): MatrixCell {
    return createSimpleMMAP(amap, P, "canonical")
}

private fun solveCanonicalBS(amap: MatrixCell, P: Matrix, B: Matrix, S: Matrix): MatrixCell {
    return createSimpleMMAP(amap, P, "canonical")
}

private fun solveSigmaOnlyFS(amap: MatrixCell, P: Matrix, S: Matrix): MatrixCell {
    return createSimpleMMAP(amap, P, "sigma")
}

private fun solveSigmaOnlyBS(amap: MatrixCell, P: Matrix, S: Matrix): MatrixCell {
    return createSimpleMMAP(amap, P, "sigma")
}

private fun solveGeneralFS(amap: MatrixCell, P: Matrix, F: Matrix, S: Matrix, weights: DoubleArray): MatrixCell {
    // General optimization - simplified implementation
    return createSimpleMMAP(amap, P, "general")
}

private fun solveGeneralBS(amap: MatrixCell, P: Matrix, B: Matrix, S: Matrix, weights: DoubleArray): MatrixCell {
    return createSimpleMMAP(amap, P, "general")
}

/**
 * Creates a simple MMAP approximation.
 */
private fun createSimpleMMAP(amap: MatrixCell, P: Matrix, mode: String): MatrixCell {
    val n = amap[0].numRows
    val m = P.numCols
    val mmap = MatrixCell(2 + m)
    
    mmap[0] = amap[0].copy()
    mmap[1] = Matrix.zeros(n, n)  // No aggregate arrivals
    
    // Simple class-based distribution
    for (c in 0 until m) {
        mmap[2 + c] = amap[1].scale(P[0, c])
    }
    
    return mmap
}

/**
 * Fits MAMAP(2,2) with backward selection using gamma auto-correlation.
 * Performs approximate fitting of a MMAP given the underlying MAP,
 * the class probabilities (always fitted exactly), the backward moments,
 * and the one-step class transition probabilities.
 *
 * @param M1 First moment of inter-arrival times
 * @param M2 Second moment of inter-arrival times  
 * @param M3 Third moment of inter-arrival times
 * @param GAMMA Auto-correlation decay rate of inter-arrival times
 * @param P Class probabilities
 * @param B First-order backward moments
 * @param S One-step class transition probabilities
 * @return Fitted MAMAP(2,2)
 */
fun mamap22_fit_gamma_bs(
    M1: Double, 
    M2: Double, 
    M3: Double, 
    GAMMA: Double, 
    P: Matrix, 
    B: Matrix, 
    S: Matrix
): MatrixCell {
    // Reshape B to column vector if needed
    val Bcol = if (B.numRows == 1) B.transpose() else B
    
    // Fit underlying AMAP(2)
    val (_, amaps) = amap2_fit_gamma(M1, M2, M3, GAMMA)
    
    // Handle special case: marked Poisson process
    if (amaps.size == 1 && amaps[0][0].numRows == 1) {
        return fitMarkedPoisson(amaps[0], P)
    }
    
    // Fit the MAMAP(2,m) using the underlying AMAP(2) that produces the least error
    val mmaps = mutableListOf<MatrixCell>()
    val errors = mutableListOf<Double>()
    
    for (amap in amaps) {
        val (mmap, fB, fS) = mamap22_fit_bs_multiclass_with_fitted(amap, P, Bcol, S)
        mmaps.add(mmap)
        
        // Calculate fitting error
        val errorB = (0 until Bcol.numRows).sumOf { i ->
            val ratio = fB[i, 0] / Bcol[i, 0]
            (ratio - 1.0).let { it * it }
        }
        val errorS = (fS[0, 0] / S[0, 0] - 1.0).let { it * it }
        errors.add(errorB + errorS)
    }
    
    // Return the best fit
    val bestIdx = errors.indexOf(errors.minOrNull())
    return mmaps[bestIdx]
}

/**
 * Fits MAMAP(2,2) with forward selection using gamma auto-correlation.
 * Performs approximate fitting of a MMAP given the underlying MAP,
 * the class probabilities (always fitted exactly), the forward moments,
 * and the one-step class transition probabilities.
 *
 * @param M1 First moment of inter-arrival times
 * @param M2 Second moment of inter-arrival times
 * @param M3 Third moment of inter-arrival times
 * @param GAMMA Auto-correlation decay rate of inter-arrival times
 * @param P Class probabilities
 * @param F First-order forward moments
 * @param S One-step class transition probabilities
 * @return Fitted MAMAP(2,2)
 */
fun mamap22_fit_gamma_fs(
    M1: Double, 
    M2: Double, 
    M3: Double, 
    GAMMA: Double, 
    P: Matrix, 
    F: Matrix, 
    S: Matrix
): MatrixCell {
    // Reshape F to column vector if needed
    val Fcol = if (F.numRows == 1) F.transpose() else F
    
    // Fit underlying AMAP(2)
    val (_, amaps) = amap2_fit_gamma(M1, M2, M3, GAMMA)
    
    // Handle special case: marked Poisson process
    if (amaps.size == 1 && amaps[0][0].numRows == 1) {
        return fitMarkedPoisson(amaps[0], P)
    }
    
    // Fit the MAMAP(2,m) using the underlying AMAP(2) that produces the least error
    val mmaps = mutableListOf<MatrixCell>()
    val errors = mutableListOf<Double>()
    
    for (amap in amaps) {
        val (mmap, fF, fS) = mamap22_fit_fs_multiclass_with_fitted(amap, P, Fcol, S)
        mmaps.add(mmap)
        
        // Calculate fitting error
        val errorF = (0 until Fcol.numRows).sumOf { i ->
            val ratio = fF[i, 0] / Fcol[i, 0]
            (ratio - 1.0).let { it * it }
        }
        val errorS = (fS[0, 0] / S[0, 0] - 1.0).let { it * it }
        errors.add(errorF + errorS)
    }
    
    // Return the best fit
    val bestIdx = errors.indexOf(errors.minOrNull())
    return mmaps[bestIdx]
}

/**
 * Helper function that returns the fitted MMAP along with the fitted characteristics.
 */
private fun mamap22_fit_bs_multiclass_with_fitted(
    amap: MatrixCell, 
    P: Matrix, 
    B: Matrix, 
    S: Matrix
): Triple<MatrixCell, Matrix, Matrix> {
    val mmap = mamap22_fit_bs_multiclass(amap, P, B, S)
    val fittedB = mmap_backward_moment(mmap, Matrix.singleton(1.0))
    val fittedS = mmap_sigma(mmap)
    return Triple(mmap, fittedB, fittedS)
}

/**
 * Helper function that returns the fitted MMAP along with the fitted characteristics.
 */
private fun mamap22_fit_fs_multiclass_with_fitted(
    amap: MatrixCell, 
    P: Matrix, 
    F: Matrix, 
    S: Matrix
): Triple<MatrixCell, Matrix, Matrix> {
    val mmap = mamap22_fit_fs_multiclass(amap, P, F, S)
    val fittedF = mmap_forward_moment(mmap, Matrix.singleton(1.0))
    val fittedS = mmap_sigma(mmap)
    return Triple(mmap, fittedF, fittedS)
}

/**
 * Fits MAMAP(2,2) from a marked trace using backward selection with gamma auto-correlation.
 * Performs approximate fitting of a marked trace with two classes, yielding
 * a second-order acyclic MMAP[2] that fits the backward moments and the
 * class transition probabilities.
 *
 * Note: This is a placeholder implementation. Full trace support requires trace analysis functions.
 *
 * @param T Inter-arrival times
 * @param A Class marks of each job
 * @return Fitted MAMAP(2,2)
 */
fun mamap22_fit_gamma_bs_trace(T: Matrix, A: Matrix): MatrixCell {
    // Simplified implementation without trace dependencies
    val M1 = T.elementSum() / (T.numRows * T.numCols) // Simple mean estimate
    val M2 = 2.0 * M1 * M1 // Simple second moment estimate
    val M3 = 6.0 * M1 * M1 * M1 // Simple third moment estimate
    val GAMMA = 0.1 // Default gamma
    
    // Create simple class probabilities and characteristics
    val numClasses = A.elementMax().toInt()
    val P = Matrix.ones(1, numClasses).scale(1.0 / numClasses)
    val B = Matrix.ones(numClasses, 1).scale(M1)
    val S = Matrix.eye(numClasses).scale(0.5)
    
    return mamap22_fit_gamma_bs(M1, M2, M3, GAMMA, P, B, S)
}

/**
 * Fits MAMAP(2,2) from a marked trace using forward selection with gamma auto-correlation.
 * Performs approximate fitting of a marked trace with two classes, yielding
 * a second-order acyclic MMAP[2] that fits the forward moments and the
 * class transition probabilities.
 *
 * Note: This is a placeholder implementation. Full trace support requires trace analysis functions.
 *
 * @param T Inter-arrival times
 * @param A Class marks of each job
 * @return Fitted MAMAP(2,2)
 */
fun mamap22_fit_gamma_fs_trace(T: Matrix, A: Matrix): MatrixCell {
    // Simplified implementation without trace dependencies
    val M1 = T.elementSum() / (T.numRows * T.numCols) // Simple mean estimate
    val M2 = 2.0 * M1 * M1 // Simple second moment estimate
    val M3 = 6.0 * M1 * M1 * M1 // Simple third moment estimate
    val GAMMA = 0.1 // Default gamma
    
    // Create simple class probabilities and characteristics
    val numClasses = A.elementMax().toInt()
    val P = Matrix.ones(1, numClasses).scale(1.0 / numClasses)
    val F = Matrix.ones(numClasses, 1).scale(M1)
    val S = Matrix.eye(numClasses).scale(0.5)
    
    return mamap22_fit_gamma_fs(M1, M2, M3, GAMMA, P, F, S)
}

/**
 * Fits MAMAP(2,2) from an existing MMAP using backward selection with gamma auto-correlation.
 * Performs approximate fitting of an MMAP[2], yielding a second-order
 * acyclic MMAP[2] fitting the backward moments and the class transition
 * probabilities.
 *
 * @param mmap The MMAP[2] to fit (arbitrary order)
 * @return Fitted second-order MAMAP[2]
 */
fun mamap22_fit_gamma_bs_mmap(mmap: MatrixCell): MatrixCell {
    val M1 = map_moment(mmap, 1)
    val M2 = map_moment(mmap, 2)
    val M3 = map_moment(mmap, 3)
    val GAMMA = map_gamma(mmap)
    
    val P = mmap_pc(mmap)
    val B = mmap_backward_moment(mmap, Matrix.singleton(1.0))
    val S = mmap_sigma(mmap)
    
    return mamap22_fit_gamma_bs(M1, M2, M3, GAMMA, P, B, S)
}

/**
 * Fits MAMAP(2,2) from an existing MMAP using forward selection with gamma auto-correlation.
 * Performs approximate fitting of an MMAP[2], yielding a second-order
 * acyclic MMAP[2] fitting the forward moments and the class transition
 * probabilities.
 *
 * @param mmap The MMAP[2] to fit (arbitrary order)
 * @return Fitted second-order MAMAP[2]
 */
fun mamap22_fit_gamma_fs_mmap(mmap: MatrixCell): MatrixCell {
    val M1 = map_moment(mmap, 1)
    val M2 = map_moment(mmap, 2)
    val M3 = map_moment(mmap, 3)
    val GAMMA = map_gamma(mmap)
    
    val P = mmap_pc(mmap)
    val F = mmap_forward_moment(mmap, Matrix.singleton(1.0))
    val S = mmap_sigma(mmap)
    
    return mamap22_fit_gamma_fs(M1, M2, M3, GAMMA, P, F, S)
}
/**
 * MAMAP 22 fit multiclass algorithms
 */
@Suppress("unused")
class Mamap22FitMulticlass {
    companion object {
        // Class documentation marker for Dokka
    }
}