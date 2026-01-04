/**
 * @file Absorbing Phase-type distribution Bernstein polynomial approximation
 * 
 * Constructs APH distributions using Bernstein exponential approximation methods.
 * Advanced technique for approximating arbitrary density functions with phase-type distributions.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import kotlin.math.ln

/**
 * Fits an Acyclic Phase-type distribution using Bernstein's approximation
 *
 * Based on: András Horváth, Enrico Vicario: Construction of Phase Type Distributions 
 * by Bernstein Exponentials. EPEW 2023: 201-215.
 *
 * @param f Function handle for the density function f(x) to approximate
 * @param order Approximation order (e.g., 100)
 * @return Pair containing (D0, D1) matrices representing the APH distribution
 */
fun aph_bernstein(f: (Double) -> Double, order: Int): Pair<Matrix, Matrix> {
    val n = order
    
    // Compute normalization constant c
    var c = 0.0
    for (i in 1..n) {
        c += f(-ln(i.toDouble() / n)) / i
    }
    
    // Create the generator matrix T (diagonal with -[1:n] and super-diagonal with [1:(n-1)])
    val T = Matrix.zeros(n, n)
    for (i in 0 until n) {
        T[i, i] = -(i + 1).toDouble()
        if (i < n - 1) {
            T[i, i + 1] = (i + 1).toDouble()
        }
    }
    
    // Compute initial probability vector alpha
    val alpha = Matrix.zeros(1, n)
    for (i in 1..n) {
        alpha[0, i - 1] = f(-ln(i.toDouble() / n)) / (i * c)
    }
    
    // Create P matrix (replicated alpha vector)
    val P = Matrix.zeros(n, n)
    for (i in 0 until n) {
        for (j in 0 until n) {
            P[i, j] = alpha[0, j]
        }
    }
    
    // Compute D0 and D1 matrices
    val D0 = T
    val D1 = T.scale(-1.0).mult(P)
    
    return Pair(D0, D1)
}
/**
 * APH  bernstein algorithms
 */
@Suppress("unused")
class AphBernsteinAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}