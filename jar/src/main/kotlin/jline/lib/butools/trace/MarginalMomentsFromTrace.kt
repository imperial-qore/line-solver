/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.trace

import kotlin.math.pow

/**
 * Returns the marginal moments of a trace.
 *
 * @param trace The trace data
 * @param K The number of moments to compute (default 5)
 * @return The (raw) moments of the trace
 */
fun marginalMomentsFromTrace(trace: DoubleArray, K: Int = 5): DoubleArray {
    val moms = DoubleArray(K)
    val n = trace.size.toDouble()

    for (i in 0 until K) {
        var sum = 0.0
        for (sample in trace) {
            sum += sample.pow(i + 1)
        }
        moms[i] = sum / n
    }

    return moms
}

/**
 * Returns the marginal moments of a weighted trace.
 *
 * @param trace The trace data
 * @param weights The weights for each sample
 * @param K The number of moments to compute (default 5)
 * @return The (raw) moments of the weighted trace
 */
fun marginalMomentsFromWeightedTrace(trace: DoubleArray, weights: DoubleArray, K: Int = 5): DoubleArray {
    val moms = DoubleArray(K)
    val totalWeight = weights.sum()

    for (i in 0 until K) {
        var sum = 0.0
        for (j in trace.indices) {
            sum += weights[j] * trace[j].pow(i + 1)
        }
        moms[i] = sum / totalWeight
    }

    return moms
}
