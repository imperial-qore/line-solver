/**
 * @file Markovian Arrival Process point process probability computation using iteration
 * 
 * Computes exact arrival probabilities using iterative numerical methods based on Neuts and Li.
 * Essential for precise transient analysis and finite-time arrival probability calculations.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import kotlin.math.*

/**
 * Probability of having exactly na arrivals within time interval t using iterative method.
 * Based on Neuts and Li, MAM1.
 * 
 * @param MAP MAP process as Array<Matrix> where MAP[0] = D0 and MAP[1] = D1
 * @param na Number of arrivals
 * @param t Time interval length
 * @param M Optional bisection parameter (auto-calculated if not provided)
 * @return Probability matrix
 */
fun map_pntiter(MAP: Array<Matrix>, na: Int, t: Double, M: Int? = null): Matrix {
    if (MAP.size != 2) {
        throw IllegalArgumentException("MAP must contain exactly 2 matrices [D0, D1]")
    }
    
    val mapMean = map_mean(MAP[0], MAP[1])
    val actualM = M ?: ceil(log2(t * 100 / mapMean)).toInt()
    
    return if (actualM < 0) {
        val result = map_pntbisect(MAP, na, t)
        result.first
    } else {
        val initialResult = map_pntbisect(MAP, na, t / 2.0.pow(actualM))
        var P = initialResult.second.toMutableList()
        
        for (i in 1..actualM) {
            val Pold = P.map { it.copy() }
            for (n in 0..na) {
                P[n] = P[n].copy()
                P[n].fill(0.0)
                
                for (j in 0..n) {
                    P[n] = P[n].add(Pold[j].mult(Pold[n - j]))
                }
            }
        }
        
        P[na]
    }
}

/**
 * Helper function implementing the bisection algorithm.
 * Returns Pair<Matrix, List<Matrix>> where first is Pnt and second is P array
 */
private fun map_pntbisect(MAP: Array<Matrix>, na: Int, t: Double): Pair<Matrix, List<Matrix>> {
    val D0 = MAP[0]
    val D1 = MAP[1]
    
    // tau = max(diag(-D0))
    var tau = 0.0
    for (i in 0 until D0.numRows) {
        val diagVal = -D0[i, i]
        if (diagVal > tau) tau = diagVal
    }
    
    val N = findN(tau, t)
    val V = Array(na + 1) { Array(N + 1) { Matrix(D0.numRows, D0.numCols) } }
    val P = Array(na + 1) { Matrix(D0.numRows, D0.numCols) }
    
    val I = Matrix.eye(D0.numRows)
    val K = D0.scale(1.0 / tau).add(I)
    val K1 = D1.scale(1.0 / tau)
    
    // Initialize
    for (n in 0..na) {
        P[n] = Matrix(D0.numRows, D0.numCols)
        for (k in 0..N) {
            V[n][k] = Matrix(D0.numRows, D0.numCols)
        }
    }
    
    // Algorithm
    V[0][0] = I.copy()
    P[0] = V[0][0].scale(br(tau, t, 0))
    
    for (n in 1..na) {
        V[n][0] = Matrix(D0.numRows, D0.numCols)
        P[n] = Matrix(D0.numRows, D0.numCols)
        
        for (k in 1..N) {
            V[n][k] = V[n][k - 1].mult(K).add(V[n - 1][k - 1].mult(K1))
            P[n] = P[n].add(V[n][k].scale(br(tau, t, k)))
        }
    }
    
    return Pair(P[na], P.toList())
}

/**
 * Find optimal N parameter for the bisection algorithm.
 */
private fun findN(tau: Double, t: Double): Int {
    val epsilon = Double.MIN_VALUE
    val Nmax = 100
    
    for (N in 1..Nmax) {
        var S = 0.0
        for (n in (N + 1)..Nmax) {
            S += br(tau, t, n)
        }
        if (S < epsilon) {
            return N
        }
    }
    
    return Nmax
}

/**
 * Bernoulli coefficient: exp(-tau*t) * (tau*t)^r / r!
 */
private fun br(tau: Double, t: Double, r: Int): Double {
    val tauT = tau * t
    return exp(-tauT) * tauT.pow(r) / factorial(r)
}

/**
 * Factorial helper function.
 */
private fun factorial(n: Int): Double {
    if (n < 0) throw IllegalArgumentException("Factorial of negative number")
    if (n == 0 || n == 1) return 1.0
    
    var result = 1.0
    for (i in 2..n) {
        result *= i
    }
    return result
}
/**
 * MAP pntiter algorithms
 */
@Suppress("unused")
class MapPntiterAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}