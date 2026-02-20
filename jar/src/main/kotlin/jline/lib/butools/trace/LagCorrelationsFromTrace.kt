/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.trace

import jline.util.matrix.Matrix

/**
 * Returns the lag-k autocorrelation of a trace.
 *
 * @param trace The trace data
 * @param K The number of lags to compute (default 3)
 * @return The lag-k autocorrelation function of the trace up to lag K
 */
fun lagCorrelationsFromTrace(trace: DoubleArray, K: Int = 3): DoubleArray {
    val n = trace.size
    val mean = trace.sum() / n

    // Compute variance
    var sumSq = 0.0
    for (sample in trace) {
        sumSq += (sample - mean) * (sample - mean)
    }
    val variance = sumSq / (n - 1)

    val acf = DoubleArray(K)
    for (k in 1..K) {
        var sum = 0.0
        for (i in 0 until n - k) {
            sum += trace[i] * trace[i + k]
        }
        val covariance = sum / (n - k) - mean * mean
        acf[k - 1] = covariance / variance
    }

    return acf
}

/**
 * Returns the lag-L joint moments of a trace.
 *
 * It is computed as Nm_{i,j} = (1/(N-L)) * sum_{k=0}^{N-L-1} x_k^i * x_{k+L}^j
 *
 * @param trace The trace data
 * @param K The joint moments are computed up to order K
 * @param L The lag at which the joint moments are computed (default 1)
 * @return Matrix of shape (K+1, K+1) containing the lag-L joint moments, starting from moment 0
 */
@JvmOverloads
fun lagkJointMomentsFromTrace(trace: DoubleArray, K: Int, L: Int = 1): Matrix {
    val n = trace.size
    val size = K + 1
    val Nm = Matrix(size, size)

    // Precompute powers of trace[0..n-L-1] and trace[L..n-1]
    // Nm[i,j] = (1/(n-L)) * sum_{k=0}^{n-L-1} trace[k]^i * trace[k+L]^j
    val count = n - L
    for (i in 0 until size) {
        for (j in 0 until size) {
            var sum = 0.0
            for (k in 0 until count) {
                sum += Math.pow(trace[k], i.toDouble()) * Math.pow(trace[k + L], j.toDouble())
            }
            Nm[i, j] = sum / count
        }
    }

    return Nm
}
