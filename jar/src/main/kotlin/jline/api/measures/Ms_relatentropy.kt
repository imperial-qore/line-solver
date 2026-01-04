/**
 * @file Relative entropy (Kullback-Leibler divergence) for discrete variables
 * 
 * Computes relative entropy KL(P||Q) = ∑ p(x) log₂(p(x)/q(x)) measuring the information
 * lost when using distribution Q to approximate distribution P. Also known as KL divergence,
 * it quantifies the dissimilarity between probability distributions in machine learning.
 *
 * @since LINE 3.0
 */
package jline.api.measures

import jline.GlobalConstants.Inf

import jline.util.matrix.Matrix
import kotlin.math.log2
import kotlin.math.max
import kotlin.math.min

/**
 * Compute relative entropy (a.k.a KL divergence) z=KL(p(x)||p(y)) of two discrete variables x and y.
 * 
 * @param x first matrix
 * @param y second matrix of the same length
 * @return relative entropy (a.k.a KL divergence) z=KL(p(x)||p(y))
 */
fun ms_relatentropy(x: Matrix, y: Matrix): Double {
    require(x.length() == y.length()) { "Input matrices must have the same length" }
    
    if (x.isEmpty) return 0.0
    
    // Convert to integer values and find min/max
    val xVals = mutableListOf<Int>()
    val yVals = mutableListOf<Int>()
    var minVal = Int.MAX_VALUE
    var maxVal = Int.MIN_VALUE
    
    for (i in 0 until x.length()) {
        val xi = x.get(i).toInt()
        val yi = y.get(i).toInt()
        xVals.add(xi)
        yVals.add(yi)
        minVal = min(minVal, min(xi, yi))
        maxVal = max(maxVal, max(xi, yi))
    }
    
    val range = maxVal - minVal + 1
    
    // Count occurrences for each value in the range
    val countsX = IntArray(range)
    val countsY = IntArray(range)
    
    for (value in xVals) {
        countsX[value - minVal]++
    }
    for (value in yVals) {
        countsY[value - minVal]++
    }
    
    val n = x.length().toDouble()
    var klDivergence = 0.0
    
    // Calculate KL divergence
    for (i in 0 until range) {
        val px = countsX[i] / n
        val py = countsY[i] / n
        
        if (px > 0 && py > 0) {
            klDivergence += px * (log2(px) - log2(py))
        } else if (px > 0 && py == 0.0) {
            // KL divergence is infinite when px > 0 but py = 0
            return Inf
        }
    }
    
    return max(0.0, klDivergence)
}
/**
 * Relatentropy algorithms
 */
@Suppress("unused")
class RelatentropyAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}