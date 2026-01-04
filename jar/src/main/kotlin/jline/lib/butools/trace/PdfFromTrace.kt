/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.trace

/**
 * Result class for PdfFromTrace containing x and y arrays.
 */
data class PdfResult(val x: DoubleArray, val y: DoubleArray)

/**
 * Returns the empirical density function of a trace.
 *
 * @param trace The trace data
 * @param intBounds The array of interval boundaries. The pdf is the
 *        number of samples falling into an interval divided by the interval length.
 * @return PdfResult containing x (center of intervals) and y (values) of the empirical PDF
 */
fun pdfFromTrace(trace: DoubleArray, intBounds: DoubleArray): PdfResult {
    val numIntervals = intBounds.size - 1
    val hist = IntArray(numIntervals)

    // Count samples in each interval
    for (sample in trace) {
        for (i in 0 until numIntervals) {
            if (sample >= intBounds[i] && sample < intBounds[i + 1]) {
                hist[i]++
                break
            }
        }
        // Handle upper boundary
        if (trace.isNotEmpty() && trace.max() == intBounds.last()) {
            if (sample == intBounds.last()) {
                hist[numIntervals - 1]++
            }
        }
    }

    val n = trace.size.toDouble()
    val x = DoubleArray(numIntervals)
    val y = DoubleArray(numIntervals)

    for (i in 0 until numIntervals) {
        val intLen = intBounds[i + 1] - intBounds[i]
        x[i] = (intBounds[i] + intBounds[i + 1]) / 2.0
        y[i] = hist[i] / intLen / n
    }

    return PdfResult(x, y)
}

/**
 * Returns the empirical density function of a weighted trace.
 *
 * @param trace The trace data
 * @param weights The weights for each sample
 * @param intBounds The array of interval boundaries.
 * @return PdfResult containing x (center of intervals) and y (values) of the empirical PDF
 */
fun pdfFromWeightedTrace(trace: DoubleArray, weights: DoubleArray, intBounds: DoubleArray): PdfResult {
    val numIntervals = intBounds.size - 1
    val hist = DoubleArray(numIntervals)

    // Sum weights in each interval
    for (j in trace.indices) {
        val sample = trace[j]
        for (i in 0 until numIntervals) {
            if (sample >= intBounds[i] && sample < intBounds[i + 1]) {
                hist[i] += weights[j]
                break
            }
        }
    }

    val totalWeight = weights.sum()
    val x = DoubleArray(numIntervals)
    val y = DoubleArray(numIntervals)

    for (i in 0 until numIntervals) {
        val intLen = intBounds[i + 1] - intBounds[i]
        x[i] = (intBounds[i] + intBounds[i + 1]) / 2.0
        y[i] = hist[i] / intLen / totalWeight
    }

    return PdfResult(x, y)
}
