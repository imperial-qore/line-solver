/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.fitting

import jline.lib.butools.mc.dtmcSolve
import jline.util.matrix.Matrix
import kotlin.math.*

/**
 * Evaluates the log-likelihood of a trace with the given PH distribution or MAP.
 * The result is divided by the length of the trace.
 *
 * @param trace The samples of the trace
 * @param alpha The initial probability vector of the PH distribution
 * @param A The transient generator of the PH distribution
 * @param prec Numerical precision used by the randomization (default 1e-14)
 * @return The log likelihood divided by the size of the trace
 */
fun likelihoodFromTracePH(
    trace: DoubleArray,
    alpha: Matrix,
    A: Matrix,
    prec: Double = 1e-14
): Double {
    val sortedTrace = trace.sorted().toDoubleArray()
    val N = A.numRows

    // Find lambda as max absolute value of diagonal
    var lambda = 0.0
    for (i in 0 until N) {
        val absVal = abs(A[i, i])
        if (absVal > lambda) lambda = absVal
    }

    // P = A/lambda + I
    val P = A.scale(1.0 / lambda).add(Matrix.eye(N))

    // a = -A * ones (closing vector)
    val a = Matrix(N, 1)
    for (i in 0 until N) {
        var rowSum = 0.0
        for (j in 0 until N) {
            rowSum += A[i, j]
        }
        a[i, 0] = -rowSum
    }

    val eps = max(prec, 10.0.pow(log10(prec) + log10(lambda)))
    val L = sortedTrace.size

    // Initialize Poisson probabilities
    val lpoi = DoubleArray(L) { -lambda * sortedTrace[it] }
    val logTrace = DoubleArray(L) { ln(sortedTrace[it]) }
    val poi = DoubleArray(L) { exp(lpoi[it]) }
    val spoi = poi.copyOf()

    // fx = poi * (alpha * a)
    val alphaA = alpha.mult(a)[0, 0]
    val fx = DoubleArray(L) { poi[it] * alphaA }

    var k = 1
    var first = 0
    var coeffv = alpha.copy()
    val maxIter = 10000

    while (first < L && k < maxIter) {
        coeffv = coeffv.mult(P)
        val coeffvA = coeffv.mult(a)[0, 0]

        for (i in first until L) {
            lpoi[i] += ln(lambda) + logTrace[i] - ln(k.toDouble())
            poi[i] = exp(lpoi[i])
            spoi[i] += poi[i]
            fx[i] += poi[i] * coeffvA
        }

        k++

        // Find new first index where spoi < 1 - eps
        var newFirst = L
        for (i in first until L) {
            if (spoi[i] < 1 - eps) {
                newFirst = i
                break
            }
        }
        first = newFirst
    }

    // Compute log-likelihood
    var logLikelihood = 0.0
    for (i in 0 until L) {
        if (fx[i] > 0) {
            logLikelihood += ln(fx[i])
        }
    }

    return logLikelihood / L
}

/**
 * Evaluates the log-likelihood of a trace with the given MAP.
 * The result is divided by the length of the trace.
 *
 * @param trace The samples of the trace
 * @param D0 The D0 matrix of the MAP
 * @param D1 The D1 matrix of the MAP
 * @param prec Numerical precision used by the randomization (default 1e-14)
 * @return The log likelihood divided by the size of the trace
 */
fun likelihoodFromTraceMAP(
    trace: DoubleArray,
    D0: Matrix,
    D1: Matrix,
    prec: Double = 1e-14
): Double {
    val N = D0.numRows
    val L = trace.size

    // Sort trace and keep original indices
    val indexed = trace.withIndex().sortedBy { it.value }
    val sortedTrace = indexed.map { it.value }.toDoubleArray()
    val ix = indexed.map { it.index }.toIntArray()

    // Find lambda
    var lambda = 0.0
    for (i in 0 until N) {
        val absVal = abs(D0[i, i])
        if (absVal > lambda) lambda = absVal
    }

    // P = D0/lambda + I
    val P = D0.scale(1.0 / lambda).add(Matrix.eye(N))

    val eps = max(prec, 10.0.pow(log10(prec) + log10(lambda)))

    // Initialize
    val lpoi = DoubleArray(L) { -lambda * sortedTrace[it] }
    val logTrace = DoubleArray(L) { ln(sortedTrace[it]) }
    val poi = DoubleArray(L) { exp(lpoi[it]) }
    val spoi = poi.copyOf()

    // fx is L x N x N (stored as list of matrices)
    val fx = Array(L) { i -> D1.scale(poi[i]) }

    var k = 1
    var first = 0
    var coeffv = D1.copy()
    val maxIter = 10000

    while (first < L && k < maxIter) {
        coeffv = P.mult(coeffv)

        for (i in first until L) {
            lpoi[i] += ln(lambda) + logTrace[i] - ln(k.toDouble())
            poi[i] = exp(lpoi[i])
            spoi[i] += poi[i]
            fx[i] = fx[i].add(coeffv.scale(poi[i]))
        }

        k++

        var newFirst = L
        for (i in first until L) {
            if (spoi[i] < 1 - eps) {
                newFirst = i
                break
            }
        }
        first = newFirst
    }

    // Compute stationary distribution
    val Pi = D0.neg().inv().mult(D1)
    val alpha = dtmcSolve(Pi)

    // Compute log-likelihood by multiplying matrices in original order
    var l = alpha.copy()
    var sc = 0

    // Create reverse mapping
    val ixrev = IntArray(L)
    for (i in 0 until L) {
        ixrev[ix[i]] = i
    }

    for (i in 0 until L) {
        l = l.mult(fx[ixrev[i]])

        // Rescale periodically to avoid numerical issues
        if ((i + 1) % 10 == 0) {
            val sumL = l.elementSum()
            if (sumL > 0) {
                val scale = ceil(log2(sumL)).toInt()
                if (scale > 1) {
                    l = l.scale(1.0 / (1 shl scale).toDouble())
                    sc += scale
                }
                if (scale < -10) {
                    val adjScale = scale + 10
                    l = l.scale(1.0 / 2.0.pow(adjScale))
                    sc += adjScale
                }
            }
        }
    }

    return (ln(l.elementSum()) + sc * ln(2.0)) / L
}
