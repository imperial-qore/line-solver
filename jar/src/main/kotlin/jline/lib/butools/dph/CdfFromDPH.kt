/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph

import jline.util.matrix.Matrix

/**
 * Returns the cumulative distribution function of a discrete phase-type distribution.
 *
 * @param alpha The initial probability vector of the discrete phase-type distribution.
 * @param A The transition probability matrix of the discrete phase-type distribution.
 * @param x Vector of non-negative integers at which to compute the CDF.
 * @return The probabilities that the discrete phase type distributed random variable
 *         is less or equal to the corresponding "x" values.
 */
fun cdfFromDPH(alpha: Matrix, A: Matrix, x: IntArray): DoubleArray {
    return cdfFromMG(alpha, A, x)
}

/**
 * Overload for DoubleArray alpha.
 */
fun cdfFromDPH(alpha: DoubleArray, A: Matrix, x: IntArray): DoubleArray {
    return cdfFromDPH(Matrix(alpha), A, x)
}

/**
 * Overload for single integer x.
 */
fun cdfFromDPH(alpha: Matrix, A: Matrix, x: Int): Double {
    return cdfFromDPH(alpha, A, intArrayOf(x))[0]
}

/**
 * Overload for DoubleArray alpha and single integer x.
 */
fun cdfFromDPH(alpha: DoubleArray, A: Matrix, x: Int): Double {
    return cdfFromDPH(Matrix(alpha), A, intArrayOf(x))[0]
}
