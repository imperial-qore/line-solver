/**
 * @file Kuiper statistical distance test
 * 
 * Implements the Kuiper test statistic, a rotation-invariant variant of the Kolmogorov-Smirnov
 * test that considers both positive and negative deviations between empirical distributions.
 * Particularly useful for circular or cyclic data where traditional tests may fail.
 *
 * @since LINE 3.0
 */
package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.min
import kotlin.math.pow

/**
 * Kuiper distance between two empirical distributions.
 * Computes the sum of the maximum positive and negative deviations between empirical CDFs.
 * This is a rotation-invariant variant of the Kolmogorov-Smirnov test.
 * 
 * @param XX first sample matrix (if 2D, each column is a sample)
 * @param YY second sample matrix (if 2D, each column is a sample)
 * @return Kuiper distance
 */
fun ms_kuiper(XX: Matrix, YY: Matrix): Double {
    // Handle 2D matrices - compute Kuiper for each column and return max
    if (XX.getNumCols() > 1 || YY.getNumCols() > 1) {
        require(XX.getNumCols() == YY.getNumCols()) { "Matrices must have the same number of columns" }
        
        var maxKuiper = 0.0
        for (col in 0 until XX.getNumCols()) {
            val xCol = XX.getColumn(col)
            val yCol = YY.getColumn(col)
            val kuiper = ms_kuiper_1D(xCol, yCol)
            maxKuiper = max(maxKuiper, kuiper)
        }
        return maxKuiper
    }
    
    // Handle 1D case
    return ms_kuiper_1D(XX, YY)
}

private fun ms_kuiper_1D(XX: Matrix, YY: Matrix): Double {
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
    val n = nx + ny
    
    // Combine samples and create indicators
    val combined = mutableListOf<Triple<Double, Double, Double>>() // (value, X_indicator, Y_indicator)
    
    // Add X values
    for (x in X) {
        combined.add(Triple(x, 1.0 / nx, 0.0))
    }
    
    // Add Y values
    for (y in Y) {
        combined.add(Triple(y, 0.0, 1.0 / ny))
    }
    
    // Sort by value
    combined.sortBy { it.first }
    
    var eCDF = 0.0
    var fCDF = 0.0
    var maxUp = 0.0    // Maximum positive deviation (F - E)
    var maxDown = 0.0  // Maximum negative deviation (F - E)
    val power = 1.0
    
    for (i in 0 until n - 1) {
        eCDF += combined[i].second
        fCDF += combined[i].third
        
        // Only update at step changes
        if (combined[i + 1].first != combined[i].first) {
            val height = fCDF - eCDF
            
            if (height > maxUp) {
                maxUp = height
            }
            if (height < maxDown) {
                maxDown = height
            }
        }
    }
    
    // Kuiper statistic is the sum of absolute values of max deviations
    return abs(maxDown).pow(power) + abs(maxUp).pow(power)
}
/**
 * Kuiper Dist algorithms
 */
@Suppress("unused")
class KuiperDistAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}