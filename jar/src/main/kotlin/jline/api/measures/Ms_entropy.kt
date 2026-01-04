package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.log2
import kotlin.math.max

/**
 * Compute entropy z=H(x) of a discrete variable x.
 * 
 * @param x a matrix representing discrete values
 * @return entropy z=H(x)
 */
fun ms_entropy(x: Matrix): Double {
    if (x.isEmpty) return 0.0
    
    // Count occurrences of each unique value
    val counts = mutableMapOf<Int, Int>()
    for (i in 0 until x.length()) {
        val value = x.get(i).toInt()
        counts[value] = counts.getOrDefault(value, 0) + 1
    }
    
    val n = x.length().toDouble()
    var entropy = 0.0
    
    // Calculate entropy
    for (count in counts.values) {
        val p = count / n
        if (p > 0) {
            entropy -= p * log2(p)
        }
    }
    
    return max(0.0, entropy)
}