/**
 * @file Markovian Arrival Process joint moment analysis
 * 
 * Computes joint moments of MAP inter-arrival times for advanced statistical characterization.
 * Essential for analyzing dependencies between consecutive arrivals in stochastic processes.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the joint moments of a Markovian Arrival Process (MAP).
 *
 * This function calculates the joint moment E[(X_a1)^i1*(X_{a1+a2})^i2*(X_{a1+a2+a3})^i3*...]
 * where X_k represents the k-th inter-arrival time.
 *
 * @param MAP The Markovian Arrival Process stored in a MatrixCell, containing the (D0, D1) matrices
 * @param a A vector (a1, a2, a3, ...) which are the subscripts of each term in the joint moment
 * @param i A vector (i1, i2, i3, ...) specifying the power of each element in the joint moment
 * @return The joint moment value
 */
fun map_joint(MAP: MatrixCell, a: IntArray, i: IntArray): Double {
    return map_joint(MAP[0], MAP[1], a, i)
}

/**
 * Computes the joint moments of a Markovian Arrival Process (MAP).
 *
 * This function calculates the joint moment E[(X_a1)^i1*(X_{a1+a2})^i2*(X_{a1+a2+a3})^i3*...]
 * where X_k represents the k-th inter-arrival time.
 *
 * @param D0 The hidden transition matrix of the MAP
 * @param D1 The visible transition matrix of the MAP
 * @param a A vector (a1, a2, a3, ...) which are the subscripts of each term in the joint moment
 * @param i A vector (i1, i2, i3, ...) specifying the power of each element in the joint moment
 * @return The joint moment value
 */
fun map_joint(D0: Matrix, D1: Matrix, a: IntArray, i: IntArray): Double {
    if (a.size != i.size) {
        throw IllegalArgumentException("Vectors a and i must have the same length")
    }
    
    // Compute the cumulative sum of vector a
    val aCumSum = IntArray(a.size)
    aCumSum[0] = a[0]
    for (k in 1 until a.size) {
        aCumSum[k] = aCumSum[k-1] + a[k]
    }
    
    // Compute P = inv(-D0) * D1
    val minusD0 = D0.scale(-1.0)
    val invD0 = minusD0.inv()
    val P = invD0.mult(D1)
    
    // Initialize joint moment computation as identity matrix (JM = 1 in MATLAB)
    var JM = Matrix.eye(D0.numRows)
    val K = a.size
    
    // Compute the product for k = 1 to K-1
    // MATLAB: JM = JM * factorial(i(k)) * invD0^(i(k)) * (P^(a(k+1)-a(k)))
    for (k in 0 until (K-1)) {
        val factorial_i_k = factorial(i[k]).toDouble()
        val invD0_power = matrixPower(invD0, i[k])
        val P_power = matrixPower(P, (aCumSum[k+1] - aCumSum[k]))
        
        // Chain multiplication: JM = JM * factorial(i[k]) * invD0^(i[k]) * P^(a[k+1] - a[k])
        JM = JM.mult(invD0_power.scale(factorial_i_k)).mult(P_power)
    }
    
    // Final computation: map_pie(MAP) * JM * factorial(i[K]) * invD0^(i[K]) * ones
    val pie = map_pie(D0, D1)
    val factorial_i_K = factorial(i[K-1]).toDouble()
    val invD0_power_final = matrixPower(invD0, i[K-1])
    val ones = Matrix.ones(P.numRows, 1)
    
    val finalResult = pie.mult(JM)
        .mult(invD0_power_final.scale(factorial_i_K))
        .mult(ones)
    
    return finalResult.get(0, 0)
}

/**
 * Computes the factorial of a non-negative integer.
 */
private fun factorial(n: Int): Long {
    if (n < 0) throw IllegalArgumentException("Factorial is not defined for negative numbers")
    if (n == 0 || n == 1) return 1
    
    var result = 1L
    for (i in 2..n) {
        result *= i
    }
    return result
}

/**
 * Computes matrix raised to an integer power.
 */
private fun matrixPower(matrix: Matrix, power: Int): Matrix {
    if (power < 0) throw IllegalArgumentException("Negative powers not supported")
    if (power == 0) return Matrix.eye(matrix.numRows)
    if (power == 1) return matrix.copy()
    
    var result = Matrix.eye(matrix.numRows)
    var base = matrix.copy()
    var exp = power
    
    while (exp > 0) {
        if (exp % 2 == 1) {
            result = result.mult(base)
        }
        base = base.mult(base)
        exp /= 2
    }
    
    return result
}
/**
 * MAP joint algorithms
 */
@Suppress("unused")
class MapJointAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}