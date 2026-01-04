package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.pow
import kotlin.math.sqrt

/**
 * Anderson-Darling distance between two empirical distributions.
 * Computes a weighted version of the Kolmogorov-Smirnov statistic.
 * 
 * @param XX first sample matrix (if 2D, each column is a sample)
 * @param YY second sample matrix (if 2D, each column is a sample)
 * @return Anderson-Darling distance
 */
fun ms_anderson_darling(XX: Matrix, YY: Matrix): Double {
    // Handle 2D matrices - compute AD for each column and return max
    if (XX.getNumCols() > 1 || YY.getNumCols() > 1) {
        require(XX.getNumCols() == YY.getNumCols()) { "Matrices must have the same number of columns" }
        
        var maxAD = 0.0
        for (col in 0 until XX.getNumCols()) {
            val xCol = XX.getColumn(col)
            val yCol = YY.getColumn(col)
            val ad = ms_anderson_darling_1D(xCol, yCol)
            maxAD = max(maxAD, ad)
        }
        return maxAD
    }
    
    // Handle 1D case
    return ms_anderson_darling_1D(XX, YY)
}

private fun ms_anderson_darling_1D(XX: Matrix, YY: Matrix): Double {
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
    var gCDF = 0.0
    val power = 2.0
    
    for (i in 0 until n - 1) {
        eCDF += combined[i].second
        fCDF += combined[i].third
        gCDF += 1.0 / n
        
        val sd = sqrt(n * gCDF * abs(1.0 - gCDF))
        val height = abs(fCDF - eCDF)
        
        // Only accumulate when we're at a step (values differ)
        if (combined[i].first != combined[i + 1].first) {
            if (sd > 0) {
                result += (height / sd).pow(power)
            }
        }
    }
    
    return result
}