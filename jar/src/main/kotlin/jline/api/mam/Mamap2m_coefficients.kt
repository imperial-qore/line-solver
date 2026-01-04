/**
 * @file Markovian Arrival MAP with Marked arrivals coefficient computation
 * 
 * Computes coefficients for MAMAP(2,m) fitting formulas in canonical forms.
 * Essential mathematical support for advanced multiclass arrival process parameter estimation.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix

/**
 * Returns the coefficients used in the direct and inverse formulas for
 * fitting a MAMAP(2,m) in first canonical form (gamma > 0).
 * 
 * @param h1 First holding time parameter of the underlying AMAP(2) with gamma > 0
 * @param h2 Second holding time parameter of the underlying AMAP(2) with gamma > 0  
 * @param r1 First transition probability parameter of the underlying AMAP(2)
 * @param r2 Second transition probability parameter of the underlying AMAP(2)
 * @return Triple of (G coefficients of p1,p11,F11,B11, U coefficients to fit, Y denominators)
 */
fun mamap2m_can1_coefficients(h1: Double, h2: Double, r1: Double, r2: Double): Triple<Matrix, Matrix, Matrix> {
    val G = Matrix.zeros(15, 1)
    val U = Matrix.zeros(12, 1)  
    val Y = Matrix.zeros(3, 1)
    
    // G coefficients
    G[0, 0] = 1.0 - r1 / (r2 * (r1 - 1.0) + 1.0)
    G[1, 0] = -(r1 * (r2 - 1.0)) / (r1 * r2 - r2 + 1.0)
    G[2, 0] = (r1 * r2) / (r1 * r2 - r2 + 1.0)
    G[3, 0] = (r1 * (r1 - 1.0)) / (r2 * (r1 - 1.0) + 1.0) - r1 + 1.0
    G[4, 0] = -(r1 * (r1 - 1.0) * (r2 - 1.0) * (r2 - 2.0)) / (r1 * r2 - r2 + 1.0)
    G[5, 0] = (r1 * r2 * (r1 - 1.0) * (r2 - 1.0)) / (r1 * r2 - r2 + 1.0)
    G[6, 0] = (r1 * r1 * (r2 - 1.0) * (r2 - 1.0)) / (r2 * (r1 - 1.0) + 1.0)
    G[7, 0] = -(r1 * r2 * (r1 + 1.0) * (r2 - 1.0)) / (r1 * r2 - r2 + 1.0)
    G[8, 0] = (r1 * r2 * r2) / (r1 * r2 - r2 + 1.0)
    G[9, 0] = h1 - (h1 * r1) / (r2 * (r1 - 1.0) + 1.0)
    G[10, 0] = -(r1 * (r2 - 1.0) * (h1 + h2 - h1 * r2)) / (r1 * r2 - r2 + 1.0)
    G[11, 0] = (r1 * r2 * (h1 + h2 - h1 * r2)) / (r1 * r2 - r2 + 1.0)
    G[12, 0] = ((h1 + h2 * r1) * (r1 - 1.0) * (r2 - 1.0)) / (r1 * r2 - r2 + 1.0)
    G[13, 0] = -(r1 * (h1 + h2 * r1) * (r2 - 1.0)) / (r1 * r2 - r2 + 1.0)
    G[14, 0] = (h2 * r1 * r2) / (r1 * r2 - r2 + 1.0)
    
    // U coefficients  
    U[0, 0] = (r1 * r2 - r2 + 1.0).let { it * it }
    U[1, 0] = -(r1 * r2 - r2 + 1.0) * (2.0 * h1 - h1 * r1 - 2.0 * h1 * r2 + 3.0 * h2 * r1 - h2 * r1 * r1 + h2 * r1 * r1 * r2 + h1 * r1 * r2 - h2 * r1 * r2)
    U[2, 0] = r1 * (r2 - 1.0) * (h1 - h2 + h2 * r1).let { it * it }
    U[3, 0] = (r1 * r2 - r2 + 1.0) * (h2 * h2 * r1 - h1 * h1 * r2 + h1 * h1 + h1 * h2 * r1 - h1 * h2 * r1 * r2)
    U[4, 0] = -r1 * (r2 - 1.0) * (r1 * r2 - r2 + 1.0) * (h1 - h2 + h2 * r1)
    U[5, 0] = r1 * (r2 - 1.0) * (h1 - h1 * r2 + h2 * r1) * (h1 - h2 + h2 * r1)
    U[6, 0] = (r1 * r2 - r2 + 1.0).let { it * it }
    U[7, 0] = -(r1 * r2 - r2 + 1.0) * (2.0 * h1 - 2.0 * h1 * r2 + h2 * r1 - h1 * r1 * r2 * r2 + h1 * r1 * r2 + h2 * r1 * r2)
    U[8, 0] = r1 * (h2 - h1 * r2).let { it * it } * (r2 - 1.0)
    U[9, 0] = (r1 * r2 - r2 + 1.0) * (h2 * h2 * r1 - h1 * h1 * r2 + h1 * h1 + h1 * h2 * r1 - h1 * h2 * r1 * r2)
    U[10, 0] = -r1 * (h2 - h1 * r2) * (r2 - 1.0) * (r1 * r2 - r2 + 1.0)
    U[11, 0] = r1 * (h2 - h1 * r2) * (r2 - 1.0) * (h1 - h1 * r2 + h2 * r1)
    
    // Y denominators
    Y[0, 0] = G[0, 0] * G[10, 0] * G[14, 0] - G[0, 0] * G[11, 0] * G[13, 0] - G[1, 0] * G[9, 0] * G[14, 0] + G[1, 0] * G[11, 0] * G[12, 0] + G[2, 0] * G[9, 0] * G[13, 0] - G[2, 0] * G[10, 0] * G[12, 0]
    Y[1, 0] = G[2, 0] * G[12, 0] - G[0, 0] * G[14, 0]
    Y[2, 0] = G[9, 0] * G[2, 0] - G[11, 0] * G[0, 0]
    
    return Triple(G, U, Y)
}

/**
 * Returns the coefficients used in the direct and inverse formulas for
 * fitting a MAMAP(2,m) in second canonical form (gamma < 0).
 * 
 * @param h1 First holding time parameter of the underlying AMAP(2) with gamma < 0
 * @param h2 Second holding time parameter of the underlying AMAP(2) with gamma < 0
 * @param r1 First transition probability parameter of the underlying AMAP(2)  
 * @param r2 Second transition probability parameter of the underlying AMAP(2)
 * @return Triple of (E coefficients of p1,p11,F11,B11, V coefficients to fit, Z denominators)
 */
fun mamap2m_can2_coefficients(h1: Double, h2: Double, r1: Double, r2: Double): Triple<Matrix, Matrix, Matrix> {
    val E = Matrix.zeros(14, 1)  // Note: dev has 15 but only uses 14
    val V = Matrix.zeros(12, 1)
    val Z = Matrix.zeros(3, 1)
    
    // E coefficients
    E[0, 0] = 1.0 - 1.0 / (r2 * (r1 - 1.0) - r1 + 2.0)
    E[1, 0] = -(r2 - 1.0) / (r1 * (r2 - 1.0) - r2 + 2.0)
    E[2, 0] = r2 / (r1 * (r2 - 1.0) - r2 + 2.0)
    E[3, 0] = (r2 - 2.0) / (r1 * (r2 - 1.0) - r2 + 2.0) - r2 + 2.0
    E[4, 0] = r2 - r2 / (r1 * (r2 - 1.0) - r2 + 2.0)
    E[5, 0] = -(r1 * (r2 - 1.0) * (r2 - 1.0)) / (r1 + r2 - r1 * r2 - 2.0)
    E[6, 0] = -r2 - (r2 * (2.0 * r2 - 3.0)) / (r1 * (r2 - 1.0) - r2 + 2.0)
    E[7, 0] = r2 * r2 / (r2 * (r1 - 1.0) - r1 + 2.0)
    E[8, 0] = h1 - h1 / (r2 * (r1 - 1.0) - r1 + 2.0)
    E[9, 0] = h1 * (r2 - 1.0) - ((r2 - 1.0) * (2.0 * h1 + h2 - h1 * r2)) / (r1 * (r2 - 1.0) - r2 + 2.0)
    E[10, 0] = (r2 * (2.0 * h1 + h2 - h1 * r2)) / (r1 * (r2 - 1.0) - r2 + 2.0) - h1 * r2
    E[11, 0] = h2 - h2 / (r2 * (r1 - 1.0) - r1 + 2.0)
    E[12, 0] = ((h1 + h2 * r1) * (r2 - 1.0)) / (r1 + r2 - r1 * r2 - 2.0)
    E[13, 0] = (h2 * r2) / (r1 * (r2 - 1.0) - r2 + 2.0)
    
    // V coefficients
    V[0, 0] = -(r1 + r2 - r1 * r2 - 2.0).let { it * it }
    V[1, 0] = -(r1 + r2 - r1 * r2 - 2.0) * (2.0 * h1 + 2.0 * h2 - h1 * r2 - h2 * r2 + h2 * r1 * r2)
    V[2, 0] = h2 * (2.0 * h1 - h1 * r2 + h2 * r1) * (r1 + r2 - r1 * r2 - 2.0)
    V[3, 0] = (r2 - 1.0) * (h1 - h2 + h2 * r1).let { it * it }
    V[4, 0] = (h1 - h2 + h2 * r1) * (2.0 * r2 - r1 * r2 + r1 * r2 * r2 - r2 * r2)
    V[5, 0] = -(h1 * r2 + h2 * r2 - h1 * r2 * r2) * (h1 - h2 + h2 * r1)
    V[6, 0] = -(r1 + r2 - r1 * r2 - 2.0).let { it * it }
    V[7, 0] = -(r1 + r2 - r1 * r2 - 2.0) * (2.0 * h1 + 2.0 * h2 - h1 * r2 - h2 * r2 + h1 * r1 * r2 * r2 - h1 * r1 * r2)
    V[8, 0] = h1 * (r1 + r2 - r1 * r2 - 2.0) * (2.0 * h2 + h1 * r1 - h2 * r2 + h1 * r1 * r2 * r2 - 2.0 * h1 * r1 * r2)
    V[9, 0] = (r2 - 1.0) * (h1 - h2 - h1 * r1 + h1 * r1 * r2).let { it * it }
    V[10, 0] = -r2 * (h1 - h2 - h1 * r1 + h1 * r1 * r2) * (r1 + r2 - r1 * r2 - 2.0)
    V[11, 0] = -r2 * (h1 + h2 - h1 * r2) * (h1 - h2 - h1 * r1 + h1 * r1 * r2)
    
    // Z denominators  
    Z[0, 0] = E[9, 0] * E[11, 0] * E[2, 0] - E[9, 0] * E[13, 0] * E[0, 0] - E[10, 0] * E[11, 0] * E[1, 0] + E[10, 0] * E[12, 0] * E[0, 0] - E[12, 0] * E[2, 0] * E[8, 0] + E[13, 0] * E[1, 0] * E[8, 0]
    Z[1, 0] = E[11, 0] * E[1, 0] - E[12, 0] * E[0, 0]
    Z[2, 0] = E[9, 0] * E[0, 0] - E[1, 0] * E[8, 0]
    
    return Triple(E, V, Z)
}
/**
 * MAMAP 2m coefficients algorithms
 */
@Suppress("unused")
class Mamap2mCoefficients {
    companion object {
        // Class documentation marker for Dokka
    }
}