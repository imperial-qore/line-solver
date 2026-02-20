/**
 * @file Markovian Arrival Process point process probability computation using quadrature
 * 
 * Computes MAP point process probabilities using ODE quadrature methods with Runge-Kutta integration.
 * Provides alternative numerical approach for arrival probability calculations with controlled accuracy.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import kotlin.math.abs

/**
 * Compute MAP point process probabilities using ODE quadrature method.
 * This is a simplified implementation using Runge-Kutta 4th order method.
 * 
 * @param MAP MAP process as Array<Matrix> where MAP[0] = D0 and MAP[1] = D1
 * @param na Maximum number of arrivals to consider
 * @param t Time interval
 * @return Pair<Matrix, Matrix> containing (Pnt, S) where:
 *         Pnt: Probability matrix for exactly na arrivals
 *         S: Verification sum matrix (should sum to identity)
 */
fun map_pntquad(MAP: Array<Matrix>, na: Int, t: Double): Pair<Matrix, Matrix> {
    if (MAP.size != 2) {
        throw IllegalArgumentException("MAP must contain exactly 2 matrices [D0, D1]")
    }
    
    val D0 = MAP[0]
    val D1 = MAP[1]
    val Ki = D0.numRows
    
    if (D0.numCols != Ki || D1.numRows != Ki || D1.numCols != Ki) {
        throw IllegalArgumentException("D0 and D1 must be square matrices of the same size")
    }
    
    // Initialize state vector P0
    val stateSize = Ki * Ki * (na + 1)
    val P0 = DoubleArray(stateSize)
    
    // Set initial condition: P0(1:Ki^2) = reshape(eye(Ki), 1, Ki^2)
    val eye = Matrix.eye(Ki)
    for (i in 0..<Ki) {
        for (j in 0..<Ki) {
            P0[i * Ki + j] = eye[i, j]
        }
    }
    
    // All other components start at zero (already initialized)
    
    // Solve ODE using 4th-order Runge-Kutta method
    val result = rungeKutta4(P0, t, Ki, na, D0, D1)
    
    // Extract results
    val Pnt = Matrix(Ki, Ki)
    val S = Matrix(Ki, Ki)
    
    // Extract all probability matrices and compute sum
    for (n in 0..na) {
        val startIdx = n * Ki * Ki
        val Pn = Matrix(Ki, Ki)
        
        for (i in 0..<Ki) {
            for (j in 0..<Ki) {
                val idx = startIdx + i * Ki + j
                Pn[i, j] = result[idx]
            }
        }
        
        S.addEq(Pn)
        
        // Keep the last one (na arrivals) as Pnt
        if (n == na) {
            for (i in 0..<Ki) {
                for (j in 0..<Ki) {
                    Pnt[i, j] = Pn[i, j]
                }
            }
        }
    }
    
    return Pair(Pnt, S)
}

/**
 * 4th-order Runge-Kutta method for solving the ODE system.
 */
private fun rungeKutta4(
    P0: DoubleArray, 
    t: Double, 
    Ki: Int, 
    na: Int, 
    D0: Matrix, 
    D1: Matrix
): DoubleArray {
    val h = t / 1000.0  // Step size
    val steps = (t / h).toInt()
    
    var P = P0.copyOf()
    var currentT = 0.0
    
    for (step in 0..<steps) {
        val k1 = pntOde(currentT, P, Ki, na, D0, D1)
        
        val P_k1 = DoubleArray(P.size)
        for (i in P.indices) P_k1[i] = P[i] + h * k1[i] / 2.0
        val k2 = pntOde(currentT + h / 2.0, P_k1, Ki, na, D0, D1)
        
        val P_k2 = DoubleArray(P.size)
        for (i in P.indices) P_k2[i] = P[i] + h * k2[i] / 2.0
        val k3 = pntOde(currentT + h / 2.0, P_k2, Ki, na, D0, D1)
        
        val P_k3 = DoubleArray(P.size)
        for (i in P.indices) P_k3[i] = P[i] + h * k3[i]
        val k4 = pntOde(currentT + h, P_k3, Ki, na, D0, D1)
        
        // Update P
        for (i in P.indices) {
            P[i] += h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
            // Ensure non-negativity
            if (P[i] < 0.0) P[i] = 0.0
        }
        
        currentT += h
    }
    
    return P
}

/**
 * ODE right-hand side function.
 * Implements the differential equation system for the probability evolution.
 */
private fun pntOde(t: Double, P: DoubleArray, Ki: Int, na: Int, D0: Matrix, D1: Matrix): DoubleArray {
    val dP = DoubleArray(P.size)
    
    // For n=0: dP(1:Ki^2) = reshape(reshape(P(1:Ki^2), Ki, Ki) * D0, Ki^2, 1)
    val P0_matrix = Matrix(Ki, Ki)
    for (i in 0..<Ki) {
        for (j in 0..<Ki) {
            P0_matrix[i, j] = P[i * Ki + j]
        }
    }
    
    val dP0_matrix = P0_matrix.mult(D0)
    for (i in 0..<Ki) {
        for (j in 0..<Ki) {
            dP[i * Ki + j] = dP0_matrix[i, j]
        }
    }
    
    // For n=1 to na
    for (n in 1..na) {
        val startIdx = n * Ki * Ki
        val prevStartIdx = (n - 1) * Ki * Ki
        
        // Extract Pn-1(t) and Pn(t)
        val Pn_1t = Matrix(Ki, Ki)
        val Pnt = Matrix(Ki, Ki)
        
        for (i in 0..<Ki) {
            for (j in 0..<Ki) {
                Pn_1t[i, j] = P[prevStartIdx + i * Ki + j]
                Pnt[i, j] = P[startIdx + i * Ki + j]
            }
        }
        
        // Compute derivative: dPn/dt = Pn * D0 + Pn-1 * D1
        val dPn_matrix = Pnt.mult(D0).add(1.0, Pn_1t.mult(D1))
        
        // Store back in dP array
        for (i in 0..<Ki) {
            for (j in 0..<Ki) {
                dP[startIdx + i * Ki + j] = dPn_matrix[i, j]
            }
        }
    }
    
    return dP
}
/**
 * MAP pntquad algorithms
 */
@Suppress("unused")
class MapPntquad {
    companion object {
        // Class documentation marker for Dokka
    }
}