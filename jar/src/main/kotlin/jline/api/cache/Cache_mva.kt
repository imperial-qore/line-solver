/**
 * @file Cache Analysis via Mean Value Analysis
 * 
 * Implements Mean Value Analysis algorithms for cache system performance
 * evaluation. Provides efficient recursive methods for computing throughput,
 * response times, and occupancy metrics in hierarchical cache architectures.
 * 
 * @since LINE 3.0
 */
package jline.api.cache

import jline.io.Ret
import jline.util.matrix.Matrix

/**
 * Exact recursive solution of the caching model.
 *
 * @param gamma - Matrix representing access factors for item i to access list j.
 * @param m     - Matrix representing the cache capacity vector.
 * @return cacheMVAReturn - An object containing cache state probabilities at equilibrium and other related metrics.
 * Reference: Giuliano Casale, Nicolas Gast. Performance analysis methods for list-based caches
 * with non-uniform access. IEEE/ACM Transactions on Networking, 2020, pp.1-18, for details.
 */

fun cache_mva(gamma: Matrix, m: Matrix): Ret.cacheMVA {
    val n = gamma.numRows
    val h = gamma.numCols
    var SS = Matrix(0, 0)
    for (l in 0..<h) {
        val arg = Matrix(m[l].toInt() + 1, 1)
        var i = 1
        while (i <= m[l] + 1) {
            arg[i - 1] = i.toDouble()
            i++
        }
        SS = Matrix.cartesian(SS, arg)
    }
    for (i in 0..<SS.numRows) {
        for (j in 0..<SS.numCols) {
            SS[i, j] = SS[i, j] - 1
        }
    }
    var pi = Matrix(SS.numRows, n)
    val pij = arrayOfNulls<Matrix>(SS.numRows)
    for (i in 0..<SS.numRows) {
        pij[i] = Matrix(n, h)
    }
    val x = Matrix(1, h)
    val E = 1
    for (s in 0..<SS.numRows) {
        val mcur = Matrix.extractRows(SS, s, s + 1, null)
        for (l in 0..<h) {
            val mcur_l = Matrix.oner(mcur, ArrayList(listOf(l)))
            val s_l = Matrix.matchrow(SS, mcur_l)
            if (s_l >= 0) {
                val one_pi = Matrix.extractRows(pi, s_l, s_l + 1, null)
                for (i in 0..<one_pi.numCols) {
                    one_pi[i] = 1 - one_pi[i]
                }
                x[l] = mcur[l] / Matrix.extractColumn(gamma, l, null).transpose().mult(one_pi.transpose())[0]
                val elMult = Matrix.extractColumn(gamma, l, null).transpose().elementMult(one_pi, null)
                for (i in 0..<n) {
                    pij[s]!![i, l] = elMult[i] * x[l]
                    pi[s, i] = pi[s, i] + pij[s]!![i, l]
                }
            }
        }
    }
    val s = Matrix.matchrow(SS, m)
    pi = Matrix.extractRows(pi, s, s + 1, null).transpose()
    val newPij = Matrix(n, h)
    for (i in 0..<n) {
        for (j in 0..<h) {
            newPij[i, j] = pij[s]!![i, j]
        }
    }
    val pi0 = Matrix(pi.numRows, pi.numCols)
    for (i in 0..<pi0.numRows) {
        for (j in 0..<pi0.numCols) {
            pi0[i, j] = 1 - pi[i, j]
        }
    }
    val u = Matrix(n, h)
    for (l in 0..<h) {
        for (k in 0..<n) {
            u[k, l] = x[l] * gamma[k, l]
        }
    }
    return Ret.cacheMVA(pi, pi0, newPij, x, u, E)
}
/**
 * Cache mva algorithms
 */
@Suppress("unused")
class CacheMvaAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}