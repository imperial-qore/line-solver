/**
 * @file Markovian Arrival MAP with Marked arrivals gamma forward-backward trace fitting
 *
 * Performs approximate fitting of a marked trace, yielding a second-order
 * acyclic MMAP that matches the class probabilities, the forward and
 * backward moments.
 *
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import jline.api.trace.*

/**
 * Performs approximate fitting of a marked trace, yielding a second-order
 * acyclic MMAP that matches the class probabilities, the forward and
 * backward moments.
 *
 * @param T Inter-arrival times as array
 * @param A Class marks of each job as array
 * @return Fitted MAMAP[m]
 */
fun mamap2m_fit_gamma_fb_trace(T: DoubleArray, A: IntArray): MatrixCell {
    // Compute moments
    val M1 = T.average()
    val M2 = T.map { it * it }.average()
    val M3 = T.map { it * it * it }.average()
    val GAMMA = trace_gamma(T)

    // Compute class characteristics
    val P = mtrace_pc(T, A).toArray1D()
    val F = mtrace_forward_moment(T, A, intArrayOf(1), 1).getColumn(0).toArray1D()
    val B = mtrace_backward_moment(T, A, 1) as DoubleArray

    return mamap2m_fit_gamma_fb(M1, M2, M3, GAMMA, P, F, B)
}

/**
 * Performs approximate fitting of a marked trace from Matrix inputs.
 *
 * @param T Inter-arrival times as Matrix
 * @param A Class marks of each job as Matrix
 * @return Fitted MAMAP[m]
 */
fun mamap2m_fit_gamma_fb_trace(T: Matrix, A: Matrix): MatrixCell {
    val tArray = T.toArray1D()
    val aArray = IntArray(A.numRows * A.numCols) { i -> A.toArray1D()[i].toInt() }
    return mamap2m_fit_gamma_fb_trace(tArray, aArray)
}

/**
 * Computes autocorrelation decay rate (gamma) from trace.
 *
 * @param T Inter-arrival times
 * @return Autocorrelation decay rate gamma
 */
fun trace_gamma(T: DoubleArray): Double {
    if (T.size < 3) return 0.0

    val mean = T.average()
    val variance = T.map { (it - mean) * (it - mean) }.average()

    if (variance < 1e-10) return 0.0

    // Compute lag-1 autocorrelation
    var cov1 = 0.0
    for (i in 0 until T.size - 1) {
        cov1 += (T[i] - mean) * (T[i + 1] - mean)
    }
    cov1 /= (T.size - 1)

    val rho1 = cov1 / variance

    // Compute lag-2 autocorrelation
    var cov2 = 0.0
    for (i in 0 until T.size - 2) {
        cov2 += (T[i] - mean) * (T[i + 2] - mean)
    }
    cov2 /= (T.size - 2)

    val rho2 = cov2 / variance

    // Gamma is approximately rho2/rho1 (autocorrelation decay rate)
    return if (Math.abs(rho1) > 1e-10) {
        Math.max(-0.99, Math.min(0.99, rho2 / rho1))
    } else {
        0.0
    }
}

/**
 * MAMAP 2m fit gamma fb trace algorithms
 */
@Suppress("unused")
class Mamap2mFitGammaFbTraceAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
