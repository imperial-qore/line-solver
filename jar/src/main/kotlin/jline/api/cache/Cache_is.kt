/**
 * @file Cache Analysis via Importance Sampling
 *
 * Implements Monte Carlo importance sampling for cache normalizing constant
 * estimation. Provides an alternative to exact enumeration (EREC) and
 * saddle-point approximation (SPM) that scales better for large item counts.
 *
 * @since LINE 3.0
 */
package jline.api.cache

import jline.io.Ret
import jline.io.line_warning
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import java.util.Random

/**
 * Estimate the normalizing constant of the cache steady state distribution
 * using Monte Carlo importance sampling.
 *
 * @param gamma Matrix representing the cache access factors (n x h).
 * @param m Matrix representing the cache capacity vector (1 x h).
 * @param samples Number of Monte Carlo samples (default: 100000).
 * @return Ret.cacheIs containing the estimated normalizing constant E and its logarithm lE.
 */
fun cache_is(gamma: Matrix, m: Matrix, samples: Int = 100000): Ret.cacheIs {
    // Remove items with zero gamma
    val rowsToKeep = BooleanArray(gamma.numRows)
    val colsToKeep = BooleanArray(gamma.numCols) { true }
    var validRows = 0
    for (i in 0 until gamma.numRows) {
        if (gamma.getRow(i).elementSum() > 0) {
            rowsToKeep[i] = true
            validRows++
        }
    }

    val filteredGamma = if (validRows < gamma.numRows) {
        gamma.getSlice(rowsToKeep, colsToKeep)
    } else {
        gamma
    }

    val n = filteredGamma.numRows
    val h = filteredGamma.numCols
    val mt = m.elementSum().toInt()

    // Edge cases
    if (n == 0 || mt == 0) {
        return Ret.cacheIs(1.0, 0.0)
    }

    if (n < mt) {
        line_warning("cache_is", "Number of items (%d) less than cache capacity (%d).", n, mt)
        return Ret.cacheIs(0.0, Double.NEGATIVE_INFINITY)
    }

    if (n == mt) {
        // All items must be in cache - use exact method
        val E = cache_erec(filteredGamma, m).value()
        return Ret.cacheIs(E, FastMath.log(E))
    }

    // Pre-compute log(gamma)
    val logGamma = Matrix(n, h)
    for (i in 0 until n) {
        for (j in 0 until h) {
            logGamma[i, j] = FastMath.log(filteredGamma[i, j] + 1e-300)
        }
    }

    // Pre-compute log(m(j)!)
    var logMFact = 0.0
    for (j in 0 until h) {
        logMFact += factln(m[j].toInt())
    }

    // Log of C(n, mt)
    val logCombinations = factln(n) - factln(mt) - factln(n - mt)

    // Log of multinomial coefficient mt! / prod(m!)
    val logMultinomial = factln(mt) - logMFact

    val lZSamples = DoubleArray(samples)
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

        // Proposal probability is 1/(C(n,mt) * multinomial(mt; m))
        val logProposal = -logCombinations - logMultinomial

        lZSamples[s] = logStateProb - logProposal
    }

    // Compute mean in log space using log-sum-exp trick
    val lE = logMeanExp(lZSamples)
    val E = FastMath.exp(lE)

    return Ret.cacheIs(E, lE)
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
    // Use Stirling's approximation for large n
    val nd = n.toDouble()
    return nd * FastMath.log(nd) - nd + 0.5 * FastMath.log(2 * FastMath.PI * nd)
}

/**
 * Compute log(mean(exp(x))) in a numerically stable way.
 */
private fun logMeanExp(x: DoubleArray): Double {
    val maxVal = x.maxOrNull() ?: return Double.NEGATIVE_INFINITY
    if (maxVal == Double.NEGATIVE_INFINITY) return Double.NEGATIVE_INFINITY

    var sum = 0.0
    for (xi in x) {
        sum += FastMath.exp(xi - maxVal)
    }

    return maxVal + FastMath.log(sum / x.size)
}

/**
 * Cache importance sampling algorithms
 */
@Suppress("unused")
class CacheIsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
