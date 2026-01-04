/**
 * @file Absorbing Phase-type distribution moment adjustment
 * 
 * Adjusts moments to ensure feasibility bounds for APH(2) fitting procedures.
 * Essential for stabilizing parameter estimation when input moments are infeasible.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import org.apache.commons.math3.util.FastMath
import kotlin.math.*

/**
 * Adjusts the second and third moments (M2 and M3) of a distribution using the default "simple" method for fitting an APH(2) distribution.
 *
 *
 * This method calls `aph2_adjust` with the "simple" adjustment method, which is designed to ensure that the moments
 * are within feasible bounds for fitting an APH(2) distribution. This adjustment helps in stabilizing the fitting process
 * by correcting moments that may otherwise be infeasible.
 *
 * @param M1 the first moment (mean) of the distribution
 * @param M2 the second moment of the distribution
 * @param M3 the third moment of the distribution
 * @return a map with keys 0 and 1, representing the adjusted M2 and M3, respectively
 */
fun aph2_adjust(M1: Double, M2: Double, M3: Double, method: String = "simple"): Map<Int, Double> {
    val tol = 1e-4
    var M2a = 0.0
    var M3a = 0.0
    if (method == "simple") {
        val scva: Double
        val M1sq = FastMath.pow(M1, 2)
        val scv = ((M2 - M1sq) / M1sq)
        if (scv < 0.5) {
            M2a = 1.5 * FastMath.pow(M1, 2)
            scva = ((M2a - M1sq) / M1sq)
        } else {
            M2a = M2
            scva = scv
        }
        if (scva < 1) {
            val lb = 3 * FastMath.pow(M1, 3) * (3 * scva - 1 + FastMath.sqrt(2.0) * FastMath.pow(1 - scva, 1.5))
            val ub = 6 * FastMath.pow(M1, 3) * scva
            M3a = if (M3 < lb) {
                lb
            } else if (M3 > ub) {
                ub
            } else {
                M3
            }
        } else {
            val lb = 1.5 * FastMath.pow(M1, 3) * FastMath.pow(1 + scva, 2)
            M3a = if (M3 < lb) {
                lb * (1 + tol)
            } else {
                M3
            }
        }
    } else if (method == "opt_char" || method == "opt_char_gads") {
        // Optimize in moment space (characteristics)
        val results = optimizeInCharacteristicsSpace(M1, M2, M3, method.contains("gads"))
        M2a = results.first
        M3a = results.second
        
    } else if (method == "opt_param" || method == "opt_param_gads") {
        // Optimize in parameter space
        val results = optimizeInParameterSpace(M1, M2, M3, method.contains("gads"))
        M2a = results.first
        M3a = results.second
        
    } else {
        throw IllegalArgumentException("Invalid method: $method")
    }

    val result: MutableMap<Int, Double> = HashMap()
    result[0] = M2a
    result[1] = M3a
    return result
}

/**
 * Optimize in characteristics space (moment space)
 * Finds M2a and M3a values that minimize distance from original M2,M3 while ensuring feasibility
 */
private fun optimizeInCharacteristicsSpace(M1: Double, M2: Double, M3: Double, useGlobal: Boolean): Pair<Double, Double> {
    val maxIter = if (useGlobal) 100 else 50
    val tolerance = 1e-6
    
    // Try both fit types and pick the best result
    val result1 = optimizeForFitType(M1, M2, M3, 1, maxIter, tolerance)
    val result2 = optimizeForFitType(M1, M2, M3, 2, maxIter, tolerance)
    
    // Calculate distances from original values
    val dist1 = sqrt((result1.first - M2).pow(2) + (result1.second - M3).pow(2))
    val dist2 = sqrt((result2.first - M2).pow(2) + (result2.second - M3).pow(2))
    
    return if (dist1 < dist2) result1 else result2
}

/**
 * Optimize in parameter space
 * Finds optimal APH(2) parameters and computes corresponding M2a and M3a
 */
private fun optimizeInParameterSpace(M1: Double, M2: Double, M3: Double, useGlobal: Boolean): Pair<Double, Double> {
    val maxIter = if (useGlobal) 100 else 50
    val tolerance = 1e-6
    val feastol = 1e-6
    val degentol = 1e-8
    
    // Initial solution that fits M1 exactly: [l2, r1]
    var bestL2 = M1
    var bestR1 = 0.5
    var bestObj = Double.MAX_VALUE
    
    // Bounds: l2 >= feastol, degentol <= r1 <= 1-degentol
    val l2Min = feastol
    val l2Max = M1 * 10 // reasonable upper bound
    val r1Min = degentol
    val r1Max = 1.0 - degentol
    
    // Simple grid search + local optimization
    val gridSize = if (useGlobal) 20 else 10
    
    for (i in 0..gridSize) {
        for (j in 0..gridSize) {
            val l2 = l2Min + (l2Max - l2Min) * i / gridSize
            val r1 = r1Min + (r1Max - r1Min) * j / gridSize
            
            // Check constraint: l2 * r1 <= M1
            if (l2 * r1 > M1 - feastol) continue
            
            val l1 = M1 - l2 * r1
            if (l1 < feastol) continue
            
            // Calculate moments for these parameters
            val xM2 = 2 * l1.pow(2) + 2 * r1 * l1 * l2 + 2 * r1 * l2.pow(2)
            val xM3 = 6 * l1.pow(3) + 6 * r1 * l1.pow(2) * l2 + 6 * r1 * l1 * l2.pow(2) + 6 * r1 * l2.pow(3)
            
            val obj = sqrt((M2 - xM2).pow(2) + (M3 - xM3).pow(2))
            
            if (obj < bestObj) {
                bestObj = obj
                bestL2 = l2
                bestR1 = r1
            }
        }
    }
    
    // Local optimization around best point
    for (iter in 0..<maxIter) {
        val stepSize = 0.01 * (1.0 - iter.toDouble() / maxIter) // Decreasing step size
        var improved = false
        
        // Try small perturbations
        val deltas = arrayOf(-stepSize, stepSize)
        
        for (dl2 in deltas) {
            for (dr1 in deltas) {
                val newL2 = (bestL2 + dl2).coerceIn(l2Min, l2Max)
                val newR1 = (bestR1 + dr1).coerceIn(r1Min, r1Max)
                
                if (newL2 * newR1 > M1 - feastol) continue
                
                val l1 = M1 - newL2 * newR1
                if (l1 < feastol) continue
                
                val xM2 = 2 * l1.pow(2) + 2 * newR1 * l1 * newL2 + 2 * newR1 * newL2.pow(2)
                val xM3 = 6 * l1.pow(3) + 6 * newR1 * l1.pow(2) * newL2 + 6 * newR1 * l1 * newL2.pow(2) + 6 * newR1 * newL2.pow(3)
                
                val obj = sqrt((M2 - xM2).pow(2) + (M3 - xM3).pow(2))
                
                if (obj < bestObj) {
                    bestObj = obj
                    bestL2 = newL2
                    bestR1 = newR1
                    improved = true
                }
            }
        }
        
        if (!improved && iter > 10) break // Early termination if no improvement
    }
    
    // Calculate final moments
    val l1 = M1 - bestL2 * bestR1
    val M2a = 2 * l1.pow(2) + 2 * bestR1 * l1 * bestL2 + 2 * bestR1 * bestL2.pow(2)
    val M3a = 6 * l1.pow(3) + 6 * bestR1 * l1.pow(2) * bestL2 + 6 * bestR1 * l1 * bestL2.pow(2) + 6 * bestR1 * bestL2.pow(3)
    
    return Pair(M2a, M3a)
}

/**
 * Optimize for a specific fit type (1 or 2)
 */
private fun optimizeForFitType(M1: Double, M2: Double, M3: Double, fitType: Int, maxIter: Int, tolerance: Double): Pair<Double, Double> {
    var bestM2 = M2
    var bestM3 = M3
    var bestObj = 0.0 // Starting objective (distance from original)
    
    // Simple gradient-free optimization
    val stepSize = min(M2, M3) * 0.01
    
    for (iter in 0..<maxIter) {
        var improved = false
        val currentStepSize = stepSize * (1.0 - iter.toDouble() / maxIter)
        
        // Try perturbations in all directions
        val deltas = arrayOf(-currentStepSize, currentStepSize)
        
        for (dM2 in deltas) {
            for (dM3 in deltas) {
                val newM2 = max(tolerance, bestM2 + dM2)
                val newM3 = max(tolerance, bestM3 + dM3)
                
                // Check feasibility for this fit type
                if (isFeasible(M1, newM2, newM3, fitType)) {
                    val obj = sqrt((M2 - newM2).pow(2) + (M3 - newM3).pow(2))
                    
                    if (iter == 0 || obj < bestObj) {
                        bestObj = obj
                        bestM2 = newM2
                        bestM3 = newM3
                        improved = true
                    }
                }
            }
        }
        
        if (!improved && iter > 10) break
    }
    
    return Pair(bestM2, bestM3)
}

/**
 * Check feasibility of M2, M3 for APH(2) fitting
 */
private fun isFeasible(M1: Double, M2: Double, M3: Double, fitType: Int): Boolean {
    try {
        val tmp0 = (8 * M1.pow(3) * M3) / 3 - 3 * M1.pow(2) * M2.pow(2) - 2 * M1 * M2 * M3 + 2 * M2.pow(3) + M3.pow(2) / 9
        
        if (tmp0 < 0) return false // Non-negative square root argument
        
        val tmp1 = 3 * sqrt(tmp0)
        val tmp2 = M3 - 3 * M1 * M2
        val tmp3 = 6 * M2 - 12 * M1.pow(2)
        
        if (abs(tmp3) < 1e-10) return false // Avoid division by zero
        
        val (l1, l2) = if (fitType == 1) {
            Pair((tmp2 + tmp1) / tmp3, (tmp2 - tmp1) / tmp3)
        } else {
            Pair((tmp2 - tmp1) / tmp3, (tmp2 + tmp1) / tmp3)
        }
        
        if (l1 < 1e-6 || l2 < 1e-6) return false // Non-negative rates
        
        val p1 = (M1 - l1) / l2
        if (p1 < 0 || p1 > 1) return false // Valid probability
        
        return true
    } catch (e: Exception) {
        return false
    }
}
/**
 * APH 2 adjust algorithms
 */
@Suppress("unused")
class Aph2AdjustAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}