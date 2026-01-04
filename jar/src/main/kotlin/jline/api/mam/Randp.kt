/**
 * @file Random sampling with relative probabilities
 * 
 * Provides random value selection based on relative probability distributions.
 * Essential utility for stochastic simulation and random sampling in MAM algorithms.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import kotlin.random.Random

/**
 * Pick random values with relative probability.
 * 
 * Returns integers in the range from 1 to PROB.size with a relative probability,
 * so that the value X is present approximately (PROB[X-1]/sum(PROB)) times.
 * 
 * All values of PROB should be equal to or larger than 0.
 *
 * @param P Probability vector (all values should be >= 0)
 * @param rows Number of rows for output matrix
 * @param cols Number of columns for output matrix (defaults to rows if not specified)
 * @return Matrix of random integers with specified probabilities
 */
fun randp(P: DoubleArray, rows: Int, cols: Int = rows): Matrix {
    require(P.all { it >= 0.0 }) { "All probabilities should be 0 or larger" }
    
    val result = Matrix.zeros(rows, cols)
    
    if (P.isEmpty() || P.sum() == 0.0) {
        // All zero probabilities - return matrix of zeros
        return result
    }
    
    // Compute cumulative distribution
    val cumsum = DoubleArray(P.size + 1)
    cumsum[0] = 0.0
    val totalSum = P.sum()
    
    for (i in P.indices) {
        cumsum[i + 1] = cumsum[i] + P[i] / totalSum
    }
    
    // Generate random samples
    for (i in 0 until rows) {
        for (j in 0 until cols) {
            val u = Random.nextDouble()
            
            // Find which bin the random number falls into
            var bin = 1
            for (k in 1 until cumsum.size) {
                if (u <= cumsum[k]) {
                    bin = k
                    break
                }
            }
            
            result[i, j] = bin.toDouble()
        }
    }
    
    return result
}

/**
 * Pick random values with relative probability (Matrix version).
 * 
 * @param P Probability matrix (flattened column-wise)
 * @param rows Number of rows for output matrix
 * @param cols Number of columns for output matrix (defaults to rows if not specified)
 * @return Matrix of random integers with specified probabilities
 */
fun randp(P: Matrix, rows: Int, cols: Int = rows): Matrix {
    // Convert Matrix to DoubleArray (column-wise flattening)
    val probArray = DoubleArray(P.numRows * P.numCols)
    var idx = 0
    
    for (j in 0 until P.numCols) {
        for (i in 0 until P.numRows) {
            probArray[idx++] = P[i, j]
        }
    }
    
    return randp(probArray, rows, cols)
}

/**
 * Generate a single random sample with relative probability.
 * 
 * @param P Probability vector
 * @return Single random integer with specified probability
 */
fun randp(P: DoubleArray): Int {
    val result = randp(P, 1, 1)
    return result[0, 0].toInt()
}

/**
 * Generate a single random sample with relative probability (Matrix version).
 * 
 * @param P Probability matrix
 * @return Single random integer with specified probability  
 */
fun randp(P: Matrix): Int {
    val result = randp(P, 1, 1)
    return result[0, 0].toInt()
}
/**
 * Randp algorithms
 */
@Suppress("unused")
class Randp {
    companion object {
        // Class documentation marker for Dokka
    }
}