/**
 * @file Cache Probability Computation via Importance Sampling
 *
 * Computes cache hit probabilities using Monte Carlo importance sampling.
 * For each item, estimates the probability of being cached at each level
 * or missing from the cache.
 *
 * @since LINE 3.0
 */
package jline.api.cache

import jline.io.line_warning
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import java.util.Random

/**
 * Computes cache hit probabilities using Monte Carlo importance sampling.
 *
 * @param gamma Matrix representing the cache access factors (n x h).
 * @param m Matrix representing the cache capacity vector (1 x h).
 * @param samples Number of Monte Carlo samples (default: 100000).
 * @return Matrix (n x h+1) where prob[i,0] = miss probability for item i,
 *         prob[i,j+1] = hit probability for item i at level j.
 */
fun cache_prob_is(gamma: Matrix, m: Matrix, samples: Int = 100000): Matrix {
    val n = gamma.numRows
    val h = gamma.numCols
    val mt = m.elementSum().toInt()

    // Initialize probability matrix
    val prob = Matrix(n, h + 1)

    // Edge cases
    if (n == 0 || mt == 0) {
        // All items miss
        for (i in 0 until n) {
            prob[i, 0] = 1.0
        }
        return prob
    }

    if (n < mt) {
        line_warning("cache_prob_is", "Number of items (%d) less than cache capacity (%d).", n, mt)
        for (i in 0 until n) {
            prob[i, 0] = 1.0
        }
        return prob
    }

    if (n == mt) {
        // All items must be in cache - use exact method
        return cache_prob_erec(gamma, m)
    }

    // Pre-compute log(gamma)
    val logGamma = Matrix(n, h)
    for (i in 0 until n) {
        for (j in 0 until h) {
            logGamma[i, j] = FastMath.log(gamma[i, j] + 1e-300)
        }
    }

    // Pre-compute log(m(j)!)
    var logMFact = 0.0
    for (j in 0 until h) {
        logMFact += factln(m[j].toInt())
    }

    // Log of C(n, mt) and multinomial
    val logCombinations = factln(n) - factln(mt) - factln(n - mt)
    val logMultinomial = factln(mt) - logMFact

    // Accumulators for importance sampling estimates
    val itemLevelWeight = Matrix(n, h)
    var totalWeight = 0.0

    val random = Random()

    for (s in 0 until samples) {
        // Sample mt items uniformly without replacement
        val selected = sampleWithoutReplacement(n, mt, random)

        // Assign items to levels uniformly at random
        val assignment = assignItemsToLevels(m, selected, random)

        // Compute log of unnormalized state probability
        var logStateProb = logMFact
        for (j in 0 until h) {
            for (item in assignment[j]) {
                logStateProb += logGamma[item, j]
            }
        }

        // Compute proposal probability
        val logProposal = -logCombinations - logMultinomial

        // Importance weight for this sample (scaled to avoid overflow)
        val logIsWeight = logStateProb - logProposal
        val isWeight = FastMath.exp(logIsWeight - 50)

        // Update accumulators
        totalWeight += isWeight
        for (j in 0 until h) {
            for (item in assignment[j]) {
                itemLevelWeight[item, j] = itemLevelWeight[item, j] + isWeight
            }
        }
    }

    // Compute probabilities from accumulated weights
    if (totalWeight > 0) {
        for (i in 0 until n) {
            var hitSum = 0.0
            for (j in 0 until h) {
                val hitProb = itemLevelWeight[i, j] / totalWeight
                prob[i, j + 1] = hitProb
                hitSum += hitProb
            }
            prob[i, 0] = FastMath.max(0.0, 1.0 - hitSum)
        }
    } else {
        // Fallback: all items miss
        for (i in 0 until n) {
            prob[i, 0] = 1.0
        }
    }

    return prob
}

/**
 * Sample k items from 0..n-1 uniformly without replacement.
 */
private fun sampleWithoutReplacement(n: Int, k: Int, random: Random): IntArray {
    val selected = IntArray(k)
    val available = (0 until n).toMutableList()

    for (i in 0 until k) {
        val idx = random.nextInt(available.size)
        selected[i] = available[idx]
        available.removeAt(idx)
    }

    return selected
}

/**
 * Assign selected items to cache levels uniformly at random.
 */
private fun assignItemsToLevels(m: Matrix, selected: IntArray, random: Random): Array<IntArray> {
    val h = m.numElements
    val assignment = Array(h) { IntArray(0) }

    // Shuffle selected items
    val shuffled = selected.toMutableList()
    for (i in shuffled.size - 1 downTo 1) {
        val j = random.nextInt(i + 1)
        val temp = shuffled[i]
        shuffled[i] = shuffled[j]
        shuffled[j] = temp
    }

    // Assign to levels
    var idx = 0
    for (j in 0 until h) {
        val levelSize = m[j].toInt()
        assignment[j] = IntArray(levelSize) { shuffled[idx + it] }
        idx += levelSize
    }

    return assignment
}

/**
 * Compute log(n!) using Stirling's approximation for large n.
 */
private fun factln(n: Int): Double {
    if (n < 0) return Double.NaN
    if (n <= 1) return 0.0
    if (n < 20) {
        var result = 0.0
        for (i in 2..n) {
            result += FastMath.log(i.toDouble())
        }
        return result
    }
    val nd = n.toDouble()
    return nd * FastMath.log(nd) - nd + 0.5 * FastMath.log(2 * FastMath.PI * nd)
}

/**
 * Cache probability importance sampling algorithms
 */
@Suppress("unused")
class CacheProbIsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
