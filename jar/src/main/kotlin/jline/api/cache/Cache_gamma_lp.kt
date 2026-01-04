/**
 * @file Cache Access Factor Computation via Linear Programming
 * 
 * Computes cache access factors using linear programming optimization methods.
 * Access factors are essential parameters that characterize request patterns
 * and routing probabilities in multi-level cache systems.
 * 
 * @since LINE 3.0
 */
package jline.api.cache

import jline.io.Ret
import jline.util.matrix.Matrix

/**
 * Computes access factors for the cache.
 *
 * @param lambda - Matrix array representing request arrival rates from users to items of individual lists.
 * @param R      - 2D Matrix array representing the reachability graph of a list for different streams and items.
 * @return cacheGammaLpReturn - An object containing access factors (gamma), the number of users (u),
 * the number of items (n), and the number of lists (h).
 */

fun cache_gamma_lp(lambda: Array<Matrix>, R: Array<Array<Matrix>>): Ret.cacheGamma {
    val u = lambda.size // Number of users
    val n = lambda[0].numRows // Number of items
    val h = lambda[0].numCols - 1 // Number of lists
    val gamma = Matrix(n, h)

    for (i in 0..<n) { // for all items
        for (j in 0..<h) { // for all levels
            // Compute gamma(i,j)
            var Rvi = Matrix(R[0][i].numRows, R[0][i].numCols)
            for (v in 0..<u) {
                Rvi = Rvi.add(1.0, R[v][i])
            }
            val Pij = ArrayList<Int>()
            Pij.add(1 + j)
            var pr_j = cache_par(Rvi, 1 + j)
            while (pr_j.size > 0) {
                Pij.add(0, pr_j[0])
                pr_j = cache_par(Rvi, pr_j[0])
            }
            if (Pij.size == 0) {
                gamma[i, j] = 0
            } else {
                gamma[i, j] = 1
                for (li in 1..<Pij.size) { // For all levels up to the current one
                    var y = 0.0
                    val l_1 = Pij[li - 1]
                    val l = Pij[li]
                    for (v in 0..<u) { // For all streams
                        for (t in 0..l_1) {
                            y += lambda[v][i, t] * R[v][i][t, l]
                        }
                    }
                    gamma[i, j] = gamma[i, j] * y
                }
            }
        }
    }
    return Ret.cacheGamma(gamma, u, n, h)
}

/**
 * Finds the parent of a given list index according to the access probabilities in an access probability matrix.
 *
 * @param R - Matrix representing the access probabilities.
 * @param j - Integer representing the list for which the parent will be found.
 * @return ArrayList<Integer> - A list containing the parent of the specified list index j.
 * @throws RuntimeException if a list has more than one parent.
</Integer> */
fun cache_par(R: Matrix, j: Int): ArrayList<Int> {
    val parent = ArrayList<Int>()
    for (i in 0..j - 1) {
        if (R[i, j] != 0.0) {
            parent.add(i)
        }
    }
    if (parent.size > 1) {
        throw RuntimeException("A cache has a list with more than one parent, but the structure must be a tree.")
    }
    return parent
}
/**
 * Cache gamma lp algorithms
 */
@Suppress("unused")
class CacheGammaLpAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}