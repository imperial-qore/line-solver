package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.log2
import kotlin.math.max
import kotlin.math.min

/**
 * Compute joint entropy z=H(x,y) of two discrete variables x and y.
 * 
 * @param x first matrix
 * @param y second matrix of the same length
 * @return joint entropy z=H(x,y)
 */
fun ms_jointentropy(x: Matrix, y: Matrix): Double {
    require(x.length() == y.length()) { "Input matrices must have the same length" }
    
    if (x.isEmpty) return 0.0
    
    // Convert to integer values and find min
    val xVals = mutableListOf<Int>()
    val yVals = mutableListOf<Int>()
    var minVal = Int.MAX_VALUE
    
    for (i in 0 until x.length()) {
        val xi = x.get(i).toInt()
        val yi = y.get(i).toInt()
        xVals.add(xi)
        yVals.add(yi)
        minVal = min(minVal, min(xi, yi))
    }
    
    // Shift values to start from 1 (like in MATLAB)
    val xShifted = xVals.map { it - minVal + 1 }
    val yShifted = yVals.map { it - minVal + 1 }
    
    // Count joint occurrences
    val jointCounts = mutableMapOf<Pair<Int, Int>, Int>()
    for (i in xVals.indices) {
        val pair = Pair(xShifted[i], yShifted[i])
        jointCounts[pair] = jointCounts.getOrDefault(pair, 0) + 1
    }
    
    val n = x.length().toDouble()
    var jointEntropy = 0.0
    
    // Calculate joint entropy
    for (count in jointCounts.values) {
        val p = count / n
        if (p > 0) {
            jointEntropy -= p * log2(p)
        }
    }
    
    return max(0.0, jointEntropy)
}