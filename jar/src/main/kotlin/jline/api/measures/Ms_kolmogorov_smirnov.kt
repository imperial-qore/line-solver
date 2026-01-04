/**
 * @file Kolmogorov-Smirnov goodness-of-fit test statistic
 * 
 * Implements the Kolmogorov-Smirnov test for determining if a sample follows 
 * a hypothesized distribution. Measures the maximum difference between empirical 
 * and theoretical cumulative distribution functions for distribution comparison.
 * 
 * @since LINE 3.0
 */
package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.abs
import kotlin.math.max

/**
 * Kolmogorov-Smirnov distance between two empirical distributions.
 * Computes the maximum absolute difference between the empirical CDFs.
 * 
 * @param XX first sample matrix (if 2D, each column is a sample)
 * @param YY second sample matrix (if 2D, each column is a sample)
 * @return KS distance
 */
fun ms_kolmogorov_smirnov(XX: Matrix, YY: Matrix): Double {
    // Handle 2D matrices - compute KS for each column and return max
    if (XX.getNumCols() > 1 || YY.getNumCols() > 1) {
        require(XX.getNumCols() == YY.getNumCols()) { "Matrices must have the same number of columns" }
        
        var maxKS = 0.0
        for (col in 0 until XX.getNumCols()) {
            val xCol = XX.getColumn(col)
            val yCol = YY.getColumn(col)
            val ks = ms_kolmogorov_smirnov_1D(xCol, yCol)
            maxKS = max(maxKS, ks)
        }
        return maxKS
    }
    
    // Handle 1D case
    return ms_kolmogorov_smirnov_1D(XX, YY)
}

private fun ms_kolmogorov_smirnov_1D(XX: Matrix, YY: Matrix): Double {
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
    
    // Sort both arrays
    val xSorted = X.sorted()
    val ySorted = Y.sorted()
    
    // Combine and sort all values
    val combined = (X + Y).sorted()
    
    // Compute empirical CDFs
    var eCDF = 0.0
    var fCDF = 0.0
    var maxDist = 0.0
    
    var xi = 0
    var yi = 0
    
    for (value in combined) {
        // Count how many X values are <= current value
        while (xi < nx && xSorted[xi] <= value) {
            eCDF += 1.0 / nx
            xi++
        }
        
        // Count how many Y values are <= current value
        while (yi < ny && ySorted[yi] <= value) {
            fCDF += 1.0 / ny
            yi++
        }
        
        val dist = abs(eCDF - fCDF)
        maxDist = max(maxDist, dist)
    }
    
    return maxDist
}
/**
 * Kolmogorov Smirnov Dist algorithms
 */
@Suppress("unused")
class KolmogorovSmirnovDistAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}