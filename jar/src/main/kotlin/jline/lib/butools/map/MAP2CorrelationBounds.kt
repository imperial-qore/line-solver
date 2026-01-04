/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map

import kotlin.math.sqrt

/**
 * Returns the upper and lower correlation bounds for a MAP(2)
 * given the three marginal moments.
 *
 * @param moms First three marginal moments of the inter-arrival times
 * @return Pair of (lower bound, upper bound) for correlation
 */
fun map2CorrelationBounds(moms: DoubleArray): Pair<Double, Double> {
    val m1 = moms[0]
    val m2 = moms[1]
    val m3 = moms[2]

    val h2 = m2 / (2.0 * m1 * m1) - 1
    val h3 = m3 / (6.0 * m1 * m1 * m1) - m2 * m2 / (4.0 * m1 * m1 * m1 * m1)
    val cv2 = m2 / m1 / m1 - 1.0

    val gub = if (h2 >= 0) {
        h2
    } else {
        -(h2 + sqrt(-h3)) * (h2 + sqrt(-h3))
    }

    val glb = if (h2 <= 0 || h3 / h2 + h2 < 1) {
        -h3 - h2 * h2
    } else {
        val temp = h3 + h2 * h2 - h2
        h2 * (temp - sqrt(temp * temp + 4.0 * h2 * h2 * h2)) / (temp + sqrt(temp * temp + 4.0 * h2 * h2 * h2))
    }

    return if (h2 >= 0) {
        Pair(glb / cv2, gub / cv2)
    } else {
        Pair(gub / cv2, glb / cv2)
    }
}
