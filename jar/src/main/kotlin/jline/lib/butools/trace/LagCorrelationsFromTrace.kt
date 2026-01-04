/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.trace

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
 * Returns the lag-k joint moments of a trace.
 *
 * @param trace The trace data
 * @param K The number of lags to compute
 * @return The lag-k joint moments of the trace
 */
fun lagkJointMomentsFromTrace(trace: DoubleArray, K: Int): DoubleArray {
    val n = trace.size
    val jointMoms = DoubleArray(K)

    for (k in 1..K) {
        var sum = 0.0
        for (i in 0 until n - k) {
            sum += trace[i] * trace[i + k]
        }
        jointMoms[k - 1] = sum / (n - k)
    }

    return jointMoms
}
