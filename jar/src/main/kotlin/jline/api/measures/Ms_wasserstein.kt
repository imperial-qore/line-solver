/**
 * @file Wasserstein (Earth Mover's) distance
 * 
 * Implements the Wasserstein distance measuring the minimum cost to transform 
 * one probability distribution into another. Also known as Earth Mover's Distance, 
 * widely used in optimal transport theory and comparing probability distributions.
 * 
 * @since LINE 3.0
 */
package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.abs
import kotlin.math.max

/**
 * Wasserstein distance (Earth Mover's Distance) between two empirical distributions.
 * Computes the 1-Wasserstein distance (p=1).
 * 
 * @param XX first sample matrix (if 2D, each column is a sample)
 * @param YY second sample matrix (if 2D, each column is a sample)
 * @return Wasserstein distance
 */
fun ms_wasserstein(XX: Matrix, YY: Matrix): Double {
    // Handle 2D matrices - compute Wasserstein for each column and return max
    if (XX.getNumCols() > 1 || YY.getNumCols() > 1) {
        require(XX.getNumCols() == YY.getNumCols()) { "Matrices must have the same number of columns" }
        
        var maxWS = 0.0
        for (col in 0 until XX.getNumCols()) {
            val xCol = XX.getColumn(col)
            val yCol = YY.getColumn(col)
            val ws = ms_wasserstein_1D(xCol, yCol)
            maxWS = max(maxWS, ws)
        }
        return maxWS
    }
    
    // Handle 1D case
    return ms_wasserstein_1D(XX, YY)
}

private fun ms_wasserstein_1D(XX: Matrix, YY: Matrix): Double {
    // Remove NaN values
    val X = mutableListOf<Double>()
    val Y = mutableListOf<Double>()
    
    for (i in 0 until XX.length()) {
        val xi = XX.get(i)
        if (!xi.isNaN()) X.add(xi)
    }
    
    for (i in 0 until YY.length()) {
        val yi = YY.get(i)
        if (!yi.isNaN()) Y.add(yi)
    }
    
    if (X.isEmpty() || Y.isEmpty()) return 0.0
    
    val nx = X.size
    val ny = Y.size
    
    // Create combined array with weights
    val combined = mutableListOf<Pair<Double, Double>>() // (value, weight difference)
    
    // Add X values with positive weight
    for (x in X) {
        combined.add(Pair(x, 1.0 / nx))
    }
    
    // Add Y values with negative weight
    for (y in Y) {
        combined.add(Pair(y, -1.0 / ny))
    }
    
    // Sort by value
    combined.sortBy { it.first }
    
    var distance = 0.0
    var cumulativeWeight = 0.0
    
    // Compute Wasserstein distance
    for (i in 0 until combined.size - 1) {
        cumulativeWeight += combined[i].second
        val width = combined[i + 1].first - combined[i].first
        distance += abs(cumulativeWeight) * width
    }
    
    return distance
}
/**
 * Wasserstein Dist algorithms
 */
@Suppress("unused")
class WassersteinDistAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}