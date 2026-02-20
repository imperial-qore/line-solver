/**
 * @file H-LRU Cache Response Time Analysis
 * 
 * Computes response time metrics for Hierarchical Least Recently Used (H-LRU)
 * cache systems. H-LRU extends traditional LRU policies to multi-level cache
 * hierarchies with configurable promotion and demotion strategies.
 * 
 * @since LINE 3.0
 */
package jline.api.cache

import de.xypron.jcobyla.Calcfc
import de.xypron.jcobyla.Cobyla
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import kotlin.math.floor
import kotlin.math.log10
import kotlin.math.pow

/**
 * Computes the characteristic time of each list in the TTL approximation of h-LRU.
 *
 * @param gamma - Matrix representing access factors for items across different levels.
 * @param m     - Matrix representing cache capacity vector.
 * @return Matrix - A matrix containing the computed TTL values.
 */

fun cache_t_hlru(gamma: Matrix, m: Matrix): Matrix {
    val n = gamma.length()
    val h = gamma.numCols


    val x = DoubleArray(h)
    for (i in 0..<h) {
        x[i] = 1.0
    }


    val rhobeg = 0.5
    val rhoend = 1e-6
    val maxFunEvals =
        500 * m.length() * 10.0.pow(floor(log10(n.toDouble()) / 4.0)).toInt() // maximum number of function evaluations

    // Define the objective function
    val objectiveFunction = Calcfc { comn, comm, comx, comcon ->
        val result = hlruTime(comx, gamma, m, n, h)
        val objVal = result.norm()
        objVal
    }

    Cobyla.findMinimum(objectiveFunction,
        h,
        h,
        x,
        rhobeg,
        rhoend,
        0,
        maxFunEvals) //1: print result; 0: no print result


    // Extract the optimized characteristic times
    val t = Matrix(1, h)
    for (i in 0..<h) {
        t[0, i] = x[i]
    }

    return t
}

fun hlruTime(x: DoubleArray, gamma: Matrix, m: Matrix, n: Int, h: Int): Matrix {
    val F = Matrix(1, h)

    for (a in 0..<h) {
        val temp1 = Matrix.ones(n, 1)
        val temp2 = Matrix(n, 1)
        val probh = Matrix(n, 1)

        for (k in 0..<n) {
            for (s in 0..a) {
                temp1[k, 0] = temp1[k, 0] * (1 - FastMath.exp(-gamma[k, s] * x[s]))
            }

            var middtemp = 1.0
            var middtemp2 = 0.0
            for (l in 0..a - 1) {
                for (s in 0..l) {
                    middtemp *= (1 - FastMath.exp(-gamma[k, s] * x[s]))
                }
                middtemp2 += middtemp
            }
            temp2[k, 0] = FastMath.exp(-gamma[k, a] * x[a]) * (1 + middtemp2)

            probh[k, 0] = temp1[k, 0] / (temp1[k, 0] + temp2[k, 0])
        }

        val sumProbh = probh.elementSum()
        F[0, a] = m[a] - sumProbh
    }

    return F
}

/**
 * Auxiliary function for cache_t_hlru inner optimization.
 *
 * @param x     - Array of double values representing certain parameters for the computation.
 * @param gamma - Matrix representing access factors for items across different levels.
 * @param m     - Matrix representing cache capacity vector.
 * @param n     - Integer representing the number of items.
 * @param h     - Integer representing the number of levels.
 * @return Matrix - A matrix containing the computed objective function.
 */

fun cache_t_hlru_aux(x: DoubleArray, gamma: Matrix, m: Matrix, n: Int, h: Int): Matrix {
    val F = Matrix(1, h)

    for (a in 0..<h) {
        val temp1 = Matrix.ones(n, 1)
        val temp2 = Matrix(n, 1)
        val probh = Matrix(n, 1)

        for (k in 0..<n) {
            for (s in 0..a) {
                temp1[k, 0] = temp1[k, 0] * (1 - FastMath.exp(-gamma[k, s] * x[s]))
            }

            var middtemp = 1.0
            var middtemp2 = 0.0
            for (l in 0..a - 1) {
                for (s in 0..l) {
                    middtemp *= (1 - FastMath.exp(-gamma[k, s] * x[s]))
                }
                middtemp2 += middtemp
            }
            temp2[k, 0] = FastMath.exp(-gamma[k, a] * x[a]) * (1 + middtemp2)

            probh[k, 0] = temp1[k, 0] / (temp1[k, 0] + temp2[k, 0])
        }

        val sumProbh = probh.elementSum()
        F[0, a] = m[a] - sumProbh
    }

    return F
}
/**
 * Cache t hlru algorithms
 */
@Suppress("unused")
class CacheTHlruAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}