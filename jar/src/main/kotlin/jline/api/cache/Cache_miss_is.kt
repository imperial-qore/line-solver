/**
 * @file Cache Miss Analysis via Importance Sampling
 *
 * Computes cache miss rates using Monte Carlo importance sampling.
 * Provides global, per-user, and per-item miss rate estimates.
 *
 * @since LINE 3.0
 */
package jline.api.cache

import jline.io.Ret
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes cache miss rates using Monte Carlo importance sampling.
 *
 * @param gamma Matrix representing the cache access factors (n x h).
 * @param m Matrix representing the cache capacity vector (1 x h).
 * @param lambda MatrixCell representing the request rates for different users or items.
 * @param samples Number of Monte Carlo samples (default: 100000).
 * @return Ret.cacheMissSpm containing miss rate metrics (M, MU, MI, pi0, lE).
 */
fun cache_miss_is(gamma: Matrix, m: Matrix, lambda: MatrixCell, samples: Int = 100000): Ret.cacheMissSpm {
    val n = gamma.numRows

    // Compute normalizing constant via importance sampling
    val isResult = cache_is(gamma, m, samples)
    val lE = isResult.lE

    // Compute hit probabilities via importance sampling
    val pij = cache_prob_is(gamma, m, samples)

    // Extract miss probabilities (first column)
    val pi0 = DoubleArray(n) { i -> pij[i, 0] }

    val u = lambda.size()

    // Per-user miss rate
    val MU = DoubleArray(u)
    for (v in 0 until u) {
        for (k in 0 until n) {
            MU[v] += lambda[v][k, 0] * pi0[k]
        }
    }

    // Per-item miss rate
    val MI = DoubleArray(n)
    for (k in 0 until n) {
        MI[k] = lambda.cellsum(k, 0) * pi0[k]
    }

    // Global miss rate (sum of per-item miss rates)
    val M = MI.sum()

    return Ret.cacheMissSpm(M, MU, MI, pi0, lE)
}

/**
 * Cache miss importance sampling algorithms
 */
@Suppress("unused")
class CacheMissIsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
