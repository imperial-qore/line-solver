/**
 * @file Cramér-von Mises statistical distance test
 * 
 * Implements the Cramér-von Mises test statistic for comparing two empirical distributions
 * by computing the integrated squared difference between their cumulative distribution functions.
 * More sensitive to differences in the middle of distributions compared to Kolmogorov-Smirnov.
 *
 * @since LINE 3.0
 */
package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.pow

/**
 * Cramér-von Mises distance between two empirical distributions.
 * Computes the integrated squared difference between empirical CDFs.
 * 
 * @param XX first sample matrix (if 2D, each column is a sample)
 * @param YY second sample matrix (if 2D, each column is a sample)
 * @return Cramér-von Mises distance
 */
fun ms_cramer_von_mises(XX: Matrix, YY: Matrix): Double {
    // Handle 2D matrices - compute CVM for each column and return max
    if (XX.getNumCols() > 1 || YY.getNumCols() > 1) {
        require(XX.getNumCols() == YY.getNumCols()) { "Matrices must have the same number of columns" }
        
        var maxCVM = 0.0
        for (col in 0 until XX.getNumCols()) {
            val xCol = XX.getColumn(col)
            val yCol = YY.getColumn(col)
            val cvm = ms_cramer_von_mises_1D(xCol, yCol)
            maxCVM = max(maxCVM, cvm)
        }
        return maxCVM
    }
    
    // Handle 1D case
    return ms_cramer_von_mises_1D(XX, YY)
}

private fun ms_cramer_von_mises_1D(XX: Matrix, YY: Matrix): Double {
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
    
    var result = 0.0
    var eCDF = 0.0
    var fCDF = 0.0
    val power = 2.0
    
    for (i in 0 until n - 1) {
        eCDF += combined[i].second
        fCDF += combined[i].third
        val height = abs(fCDF - eCDF)
        
        // Accumulate squared difference when values differ
        if (combined[i + 1].first != combined[i].first) {
            result += height.pow(power)
        }
    }
    
    return result
}
/**
 * Cramer Von Mises algorithms
 */
@Suppress("unused")
class CramerVonMisesAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}