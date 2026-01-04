/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph

import jline.util.matrix.Matrix

/**
 * Returns the probability mass function of a discrete phase-type distribution.
 *
 * @param alpha The initial probability vector of the discrete phase-type distribution.
 *        The sum of the entries of alpha is less or equal to 1.
 * @param A The transient generator matrix of the discrete phase-type distribution.
 * @param x Vector of non-negative integers at which to compute the PMF.
 * @return The probabilities that the discrete phase type distributed random variable
 *         takes the corresponding "x" values.
 */
fun pmfFromDPH(alpha: Matrix, A: Matrix, x: IntArray): DoubleArray {
    return pmfFromMG(alpha, A, x)
}

/**
 * Overload for DoubleArray alpha.
 */
fun pmfFromDPH(alpha: DoubleArray, A: Matrix, x: IntArray): DoubleArray {
    return pmfFromDPH(Matrix(alpha), A, x)
}

/**
 * Overload for single integer x.
 */
fun pmfFromDPH(alpha: Matrix, A: Matrix, x: Int): Double {
    return pmfFromDPH(alpha, A, intArrayOf(x))[0]
}

/**
 * Overload for DoubleArray alpha and single integer x.
 */
fun pmfFromDPH(alpha: DoubleArray, A: Matrix, x: Int): Double {
    return pmfFromDPH(Matrix(alpha), A, intArrayOf(x))[0]
}
