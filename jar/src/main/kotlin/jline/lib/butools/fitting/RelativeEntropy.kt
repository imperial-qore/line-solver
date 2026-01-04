/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.fitting

import kotlin.math.abs
import kotlin.math.ln

/**
 * Returns the relative entropy (aka Kullback-Leibler divergence) of two vectors.
 *
 * @param p1 The first vector
 * @param p2 The second vector
 * @return The relative entropy calculated as sum(p1_i * |log(p1_i/p2_i)|)
 */
fun relativeEntropy(p1: DoubleArray, p2: DoubleArray): Double {
    require(p1.size == p2.size) { "Vectors must have the same length" }

    var re = 0.0
    for (i in p1.indices) {
        if (p1[i] > 0.0 && p2[i] > 0.0) {
            re += p1[i] * abs(ln(p1[i] / p2[i]))
        }
    }
    return re
}

/**
 * Returns the empirical relative entropy using trace data.
 *
 * @param p1 The first vector (from empirical data)
 * @param p2 The second vector (from model)
 * @return The relative entropy
 */
fun empiricalRelativeEntropy(p1: DoubleArray, p2: DoubleArray): Double {
    return relativeEntropy(p1, p2)
}
