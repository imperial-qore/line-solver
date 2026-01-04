package jline.api.mc

import jline.util.matrix.Matrix

/**
 * Result class for DTMC uniformization analysis containing the probability vector and maximum iterations used.
 */
data class DtmcUniformizationResult(
    val pi: Matrix,
    val kmax: Int
)

/**
 * Compute the transient probability distribution of a DTMC using uniformization.
 * 
 * This method converts the DTMC transition matrix to a CTMC infinitesimal generator
 * and then applies the CTMC uniformization algorithm.
 *
 * @param pi0 Initial state probability vector of the DTMC
 * @param P Transition matrix of the DTMC
 * @param t Time point for transient analysis (default: 10000.0)
 * @param tol Tolerance for convergence (default: 1e-12)
 * @param maxiter Maximum number of iterations (default: 100)
 * @return DtmcUniformizationResult containing the probability vector at time t and maximum iterations used
 */
fun dtmc_uniformization(
    pi0: Matrix,
    P: Matrix,
    t: Double = 1e4,
    tol: Double = 1e-12,
    maxiter: Int = 100
): DtmcUniformizationResult {
    // Convert DTMC transition matrix to CTMC infinitesimal generator
    val Q = ctmc_makeinfgen(P)
    
    // The existing ctmc_uniformization only returns pi, not kmax
    // So we need to implement the MATLAB logic: [pi,kmax] = ctmc_uniformization(pi0,Q,t,tol,maxiter)
    val (pi, kmax) = ctmc_uniformization_with_kmax(pi0, Q, t, tol, maxiter)
    
    return DtmcUniformizationResult(pi, kmax)
}

/**
 * Extended CTMC uniformization that returns both the probability vector and kmax.
 * This is a modified version of the existing ctmc_uniformization to also return kmax.
 */
private fun ctmc_uniformization_with_kmax(
    pi0: Matrix,
    Q: Matrix,
    t: Double,
    tol: Double,
    maxiter: Int
): Pair<Matrix, Int> {
    val n = Q.numCols
    var q = 0.0
    
    // Find the maximum diagonal element (in absolute value) and scale by 1.1
    for (i in 0..<n) {
        q = kotlin.math.max(q, 1.1 * kotlin.math.abs(Q[i, i]))
    }
    
    // Create the uniformized matrix Qs = I + Q/q
    val Qs = Matrix.eye(n)
    Qs.addEq(1.0 / q, Q)
    
    // Find the truncation point kmax
    var k = 0
    var s = 1.0
    var r = 1.0
    var iter = 0
    var kmax = 1
    
    while (iter < maxiter) {
        iter++
        k++
        r = r * (q * t) / k
        s = s + r
        if (1 - kotlin.math.exp(-q * t) * s <= tol) {
            kmax = k
            break
        }
    }
    
    // Compute the transient probability using uniformization
    // MATLAB: pi=pi0*(exp(-q*t));
    var pi = Matrix(1, n)
    pi0.scaleEq(kotlin.math.exp(-q * t), pi)
    
    // MATLAB: P=pi0;
    val P = Matrix(pi0)
    
    // MATLAB: ri=exp(-q*t);
    var ri = kotlin.math.exp(-q * t)
    
    // MATLAB: for j=1:kmax ... end
    for (j in 1..kmax) {
        // MATLAB: P=P*Qs;
        P.multEq(Qs)
        // MATLAB: ri=ri*(q*t/j);
        ri = ri * (q * t / j)
        // MATLAB: pi=pi+ri*P;
        pi = pi.add(ri, P)
    }
    
    return Pair(pi, kmax)
}
/**
 * DTMC uniformization algorithms
 */
@Suppress("unused")
class DtmcUniformizationAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}