/**
 * @file Marked Markovian Arrival Process temporal scaling operations
 * 
 * Rescales MMAP inter-arrival distributions to achieve specified mean values.
 * Essential for model calibration and parameter adjustment in multiclass systems.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Changes the mean inter-arrival time of a Markovian Arrival Process with marked arrivals (MMAP).
 *
 *
 * This method scales the MMAP to achieve a desired mean inter-arrival time. If a single new mean (M) is provided,
 * the method scales all matrices (D0, D1, and the marking matrices) uniformly. If a matrix M is provided with different
 * desired means for each class, the method scales each class independently.
 *
 * @param MMAP the original MMAP to be scaled
 * @param M    a matrix containing the new mean inter-arrival times; can be a single value or a vector with one value per class
 * @param maxIter maximum number of iterations for the optimization algorithm (default: 30)
 * @return the scaled MMAP with adjusted mean inter-arrival times
 */

fun mmap_scale(MMAP: MatrixCell, M: Matrix, maxIter: Int = 30): MatrixCell? {
    val C = MMAP.size() - 2
    var SCALED: MatrixCell? = MatrixCell(2 + C)
    
    if (M.length() == 1) {
        // Single scaling factor case - uniform scaling
        val MOLD = map_mean(MMAP[0], MMAP[1])
        val ratio = MOLD / M[0]

        val D0 = Matrix(MMAP[0])
        D0.scaleEq(ratio)
        SCALED!![0] = D0
        val D1 = Matrix(MMAP[1])
        D1.scaleEq(ratio)
        SCALED[1] = D1

        for (c in 0..<C) {
            val a = Matrix(MMAP[2 + c])
            a.scaleEq(ratio)
            SCALED[2 + c] = a
        }
    } else {
        // Multiple scaling factors case - requires iterative approximation
        SCALED!![0] = MMAP[0].copy()
        SCALED[1] = Matrix(MMAP[0].numRows, MMAP[0].numCols, MMAP[0].numRows * MMAP[0].numCols)
        val l = mmap_count_lambda(MMAP)
        
        // Initial heuristic approximation
        for (c in 0..<C) {
            if (l[c] > 0) {
                val a = MMAP[2 + c].copy()
                a.scaleEq((1 / M[c]) / l[c])
                SCALED[2 + c] = a
                SCALED[1] = SCALED[1].add(1.0, a)
            } else {
                SCALED[2 + c] =
                    Matrix(MMAP[2 + c].numRows, MMAP[2 + c].numCols, MMAP[2 + c].numRows * MMAP[2 + c].numCols)
            }
        }

        SCALED = mmap_normalize(SCALED)
        
        // Iterative approximation to refine the solution
        if (maxIter > 0) {
            try {
                // Optimize scaling factors using a simple iterative approach
                val x = optimizeScalingFactors(MMAP, M, maxIter)
                
                // Apply optimized scaling factors
                SCALED!![0] = MMAP[0].copy()
                SCALED[1] = Matrix(MMAP[0].numRows, MMAP[0].numCols, MMAP[0].numRows * MMAP[0].numCols)
                
                for (c in 0..<C) {
                    val a = MMAP[2 + c].copy()
                    a.scaleEq(x[c])
                    SCALED[2 + c] = a
                    SCALED[1] = SCALED[1].add(1.0, a)
                }
                
                SCALED = mmap_normalize(SCALED)
            } catch (e: Exception) {
                throw RuntimeException("The input MMAP is invalid.", e)
            }
        }
    }
    return SCALED
}

/**
 * Optimize scaling factors using a simple iterative approach
 * This implements a basic optimization algorithm similar to MATLAB's fminsearchbnd
 * @param MMAP Original MMAP
 * @param M Desired mean inter-arrival times
 * @param maxIter Maximum number of iterations
 * @return Optimized scaling factors
 */
private fun optimizeScalingFactors(MMAP: MatrixCell, M: Matrix, maxIter: Int): DoubleArray {
    val C = MMAP.size() - 2
    var x = DoubleArray(C) { 1.0 } // Initial guess: all scaling factors = 1
    val tolFun = 1e-2
    val stepSize = 0.1
    val minBound = 1e-6
    
    var bestX = x.copyOf()
    var bestObjective = objectiveFunction(x, M, MMAP)
    
    for (iter in 0..<maxIter) {
        val currentObjective = objectiveFunction(x, M, MMAP)
        
        // Check for convergence
        if (currentObjective < tolFun) {
            break
        }
        
        // Simple gradient-free optimization using coordinate descent
        for (c in 0..<C) {
            val originalValue = x[c]
            
            // Try increasing the scaling factor
            x[c] = maxOf(minBound, originalValue + stepSize)
            val objPlus = objectiveFunction(x, M, MMAP)
            
            // Try decreasing the scaling factor
            x[c] = maxOf(minBound, originalValue - stepSize)
            val objMinus = objectiveFunction(x, M, MMAP)
            
            // Choose the best direction
            when {
                objPlus < currentObjective && objPlus <= objMinus -> {
                    x[c] = maxOf(minBound, originalValue + stepSize)
                }
                objMinus < currentObjective && objMinus < objPlus -> {
                    x[c] = maxOf(minBound, originalValue - stepSize)
                }
                else -> {
                    x[c] = originalValue // No improvement
                }
            }
        }
        
        val newObjective = objectiveFunction(x, M, MMAP)
        if (newObjective < bestObjective) {
            bestObjective = newObjective
            bestX = x.copyOf()
        }
    }
    
    return bestX
}

/**
 * Objective function that measures the error between desired and actual arrival rates
 * This mirrors the MATLAB objfun function
 * @param x Scaling factors for each class
 * @param M Desired mean inter-arrival times
 * @param MMAP Original MMAP
 * @return Sum of normalized errors
 */
private fun objectiveFunction(x: DoubleArray, M: Matrix, MMAP: MatrixCell): Double {
    val C = MMAP.size() - 2
    val SCALED = MatrixCell(2 + C)
    
    // Build scaled MMAP with current scaling factors
    SCALED[0] = MMAP[0].copy()
    SCALED[1] = Matrix(MMAP[0].numRows, MMAP[0].numCols, MMAP[0].numRows * MMAP[0].numCols)
    
    for (c in 0..<C) {
        val a = MMAP[2 + c].copy()
        a.scaleEq(x[c])
        SCALED[2 + c] = a
        SCALED[1] = SCALED[1].add(1.0, a)
    }
    
    val normalizedSCALED = mmap_normalize(SCALED)
    if (normalizedSCALED == null) {
        return Double.MAX_VALUE // Invalid MMAP
    }
    
    val l = mmap_count_lambda(normalizedSCALED)
    
    // Calculate sum of errors: ||1/M(c) - l(c)||
    var f = 0.0
    for (c in 0..<C) {
        val targetRate = 1.0 / M[c]
        val actualRate = l[c]
        f += kotlin.math.abs(targetRate - actualRate)
    }
    
    return f
}

/**
 * Overloaded function for backward compatibility
 * @param MMAP the original MMAP to be scaled
 * @param M a matrix containing the new mean inter-arrival times
 * @return the scaled MMAP with adjusted mean inter-arrival times
 */
fun mmap_scale(MMAP: MatrixCell, M: Matrix): MatrixCell? {
    return mmap_scale(MMAP, M, 30)
}
/**
 * MMAP scale algorithms
 */
@Suppress("unused")
class MmapScale {
    companion object {
        // Class documentation marker for Dokka
    }
}