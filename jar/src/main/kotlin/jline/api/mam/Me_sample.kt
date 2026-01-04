package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import java.util.Random
import kotlin.math.sqrt

/**
 * ME sampling algorithms using inverse CDF interpolation.
 *
 * Provides methods for generating random samples from Matrix Exponential (ME) distributions
 * using the inverse CDF (Cumulative Distribution Function) method with linear interpolation.
 *
 * This approach is more accurate than phase-based sampling for ME distributions with
 * non-normalized alpha vectors, as it directly uses the CDF rather than simulating phases.
 *
 * @since LINE 3.0
 */

/**
 * Generates random samples from a Matrix Exponential (ME) distribution using inverse CDF interpolation.
 *
 * Algorithm:
 * 1. Precompute a dense grid of CDF values from 0 to mean + 10*sigma
 * 2. For each sample, generate u ~ Uniform(0,1)
 * 3. Find the corresponding time value by:
 *    a. Binary search in the CDF grid to locate u
 *    b. Linear interpolation between adjacent grid points
 *
 * @param alpha The initial vector of the ME distribution
 * @param A The matrix parameter of the ME distribution
 * @param n The number of samples to generate
 * @param random The random number generator to use
 * @return Array of n samples from the ME distribution
 */
fun me_sample(alpha: Matrix, A: Matrix, n: Long, random: Random): DoubleArray {
    val samples = DoubleArray(n.toInt())
    val nPhases = A.numRows

    // Compute mean and standard deviation for grid bounds
    val mean = me_mean(alpha, A)
    val variance = me_var(alpha, A)
    val sigma = sqrt(variance)

    // Create grid from 0 to mean + 10*sigma
    val gridSize = 1000
    val maxT = mean + 10.0 * sigma
    val tGrid = DoubleArray(gridSize)
    val cdfGrid = DoubleArray(gridSize)

    val e = Matrix.ones(nPhases, 1)

    // Precompute CDF at grid points
    // CDF(t) = 1 - alpha * exp(A*t) * e
    for (i in 0 until gridSize) {
        val t = i * maxT / (gridSize - 1)
        tGrid[i] = t

        if (t == 0.0) {
            cdfGrid[i] = 0.0
        } else {
            // Compute exp(A*t) using matrix exponential
            val At = A.scale(t)
            val expAt = At.expm_higham()  // Matrix exponential
            val result = alpha.mult(expAt).mult(e).get(0, 0)
            cdfGrid[i] = 1.0 - result
        }
    }

    // Generate samples using inverse CDF with binary search and interpolation
    for (i in 0 until n.toInt()) {
        val u = random.nextDouble()

        // Handle edge cases
        if (u <= cdfGrid[0]) {
            samples[i] = tGrid[0]
            continue
        }
        if (u >= cdfGrid[gridSize - 1]) {
            samples[i] = tGrid[gridSize - 1]
            continue
        }

        // Binary search to find interval containing u
        var low = 0
        var high = gridSize - 1

        while (low < high - 1) {
            val mid = (low + high) / 2
            if (cdfGrid[mid] < u) {
                low = mid
            } else {
                high = mid
            }
        }

        // Linear interpolation between tGrid[low] and tGrid[high]
        if (cdfGrid[high] > cdfGrid[low]) {
            val frac = (u - cdfGrid[low]) / (cdfGrid[high] - cdfGrid[low])
            samples[i] = tGrid[low] + frac * (tGrid[high] - tGrid[low])
        } else {
            // CDF values are equal (shouldn't happen, but handle gracefully)
            samples[i] = tGrid[low]
        }
    }

    return samples
}

/**
 * Generates random samples from an ME distribution using matrices stored in a MatrixCell.
 *
 * Note: For ME distributions stored in MatrixCell format {D0=A, D1=-A*e*alpha'},
 * we need to reconstruct alpha. However, for sampling purposes, it's better to
 * use the explicit (alpha, A) overload.
 *
 * @param ME The Matrix Exponential distribution stored in a MatrixCell
 * @param n The number of samples to generate
 * @param random The random number generator to use
 * @return Array of n samples from the ME distribution
 */
fun me_sample(ME: MatrixCell, n: Long, random: Random): DoubleArray {
    // For now, delegate to map_sample which works with the process representation
    // This is less accurate for non-normalized alpha, but works with the MatrixCell format
    // TODO: Extract alpha from the process representation for better accuracy
    return map_sample(ME[0], ME[1], n, random)
}

/**
 * ME sampling algorithms documentation marker for Dokka.
 */
@Suppress("unused")
class MeSample {
    companion object {
        // Class documentation marker for Dokka
    }
}
