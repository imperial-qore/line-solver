/**
 * @file Marked Markovian Arrival Process exponential distribution construction
 * 
 * Constructs MMAP with exponential inter-arrival distributions for each marked class.
 * Used for modeling independent multiclass Poisson arrival processes.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Fits an order-n Markovian Arrival Process with marked arrivals (MMAP) based on the given arrival rates for each job class.
 *
 *
 * This method constructs an MMAP with `n` states for each job class, where each class has its own arrival rate specified
 * by the matrix `lambda`. The resulting MMAP includes a D0 matrix for the hidden transitions, a D1 matrix for visible transitions,
 * and additional matrices for each job class marking the type of arrival.
 *
 * @param lambda a 1xK matrix describing the arrival rates for each job class, where K is the number of job classes
 * @param n      the number of states in the MMAP
 * @return a MatrixCell representing the MMAP with the specified characteristics
 */

fun mmap_exponential(lambda: Matrix, n: Int): MatrixCell? {
    val K = lambda.length()
    val MMAP = MatrixCell(2 + K)
    MMAP[0] = Matrix(n, n, n xor 2)
    MMAP[1] = Matrix(n, n, n xor 2)
    for (k in 0..<K) {
        val a = lambda[0, k]
        val m = Matrix(n, n, n xor 2)
        for (i in 0..<n) {
            m[i, n - 1 - i] = a
        }
        MMAP[2 + k] = m
        MMAP[1] = MMAP[1].add(1.0, m)
    }
    return mmap_normalize(MMAP)
}

/**
 * Fits a Markovian Arrival Process with marked arrivals (MMAP) with a single state based on the given arrival rates for each job class.
 *
 *
 * This method constructs an MMAP with a single state (`n = 1`) for each job class, using the arrival rates specified by the matrix `lambda`.
 * It is a simplified version of the `mmap_exponential` function that defaults to one state per job class.
 *
 * @param lambda a 1xK matrix describing the arrival rates for each job class, where K is the number of job classes
 * @return a MatrixCell representing the MMAP with the specified arrival rates and a single state
 */

fun mmap_exponential(lambda: Matrix): MatrixCell? {
    return mmap_exponential(lambda, 1)
}
/**
 * MMAP exponential algorithms
 */
@Suppress("unused")
class MmapExponentialAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}