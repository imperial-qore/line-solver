/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.trace

/**
 * Result class for CdfFromTrace containing x and y arrays.
 */
data class CdfResult(val x: DoubleArray, val y: DoubleArray)

/**
 * Returns the empirical distribution function of the trace.
 *
 * @param trace The trace data
 * @return CdfResult containing x (points) and y (values) of the empirical CDF
 */
fun cdfFromTrace(trace: DoubleArray): CdfResult {
    val sorted = trace.sorted()
    val n = trace.size

    val x = sorted.toDoubleArray()
    val y = DoubleArray(n) { i -> i.toDouble() / (n - 1) }

    return CdfResult(x, y)
}

/**
 * Returns the empirical distribution function of a weighted trace.
 *
 * @param trace The trace data
 * @param weights The weights for each sample
 * @return CdfResult containing x (points) and y (values) of the empirical CDF
 */
fun cdfFromWeightedTrace(trace: DoubleArray, weights: DoubleArray): CdfResult {
    // Sort trace and weights together by trace values
    val indices = trace.indices.sortedBy { trace[it] }
    val x = DoubleArray(trace.size) { i -> trace[indices[i]] }
    val sortedWeights = DoubleArray(trace.size) { i -> weights[indices[i]] }

    val totalWeight = sortedWeights.sum()
    var cumWeight = 0.0
    val y = DoubleArray(trace.size)
    for (i in trace.indices) {
        cumWeight += sortedWeights[i]
        y[i] = cumWeight / totalWeight
    }

    return CdfResult(x, y)
}
