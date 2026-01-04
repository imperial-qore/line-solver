/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.fitting

/**
 * Returns the squared difference between two vectors.
 *
 * @param p1 The first vector
 * @param p2 The second vector
 * @return The squared difference calculated as sum((p1_i - p2_i)^2)
 */
fun squaredDifference(p1: DoubleArray, p2: DoubleArray): Double {
    require(p1.size == p2.size) { "Vectors must have the same length" }

    var sd = 0.0
    for (i in p1.indices) {
        val diff = p1[i] - p2[i]
        sd += diff * diff
    }
    return sd
}

/**
 * Returns the empirical squared difference using trace data.
 *
 * @param p1 The first vector (from empirical data)
 * @param p2 The second vector (from model)
 * @return The squared difference
 */
fun empiricalSquaredDifference(p1: DoubleArray, p2: DoubleArray): Double {
    return squaredDifference(p1, p2)
}
