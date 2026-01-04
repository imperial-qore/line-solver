/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 *
 * Reference:
 * G Horváth, M Telek, "A minimal representation of Markov arrival processes
 * and a moments matching method," Performance Evaluation 64:(9-12) pp. 1153-1168. (2007)
 */
package jline.lib.butools.dph

import jline.lib.butools.reptrans.findMarkovianRepresentation
import jline.util.matrix.Matrix
import kotlin.math.min

/**
 * Obtains a Markovian representation of a matrix-geometric distribution
 * of the same size, if possible.
 *
 * @param alpha The initial vector of the matrix-geometric distribution.
 * @param A The matrix parameter of the matrix-geometric distribution.
 * @param precision A representation is considered to be a Markovian one
 *        if it is closer than the precision (default 1e-14)
 * @return The MGRepresentation containing beta (initial probability vector) and B (transition probability matrix)
 */
fun dphFromMG(alpha: Matrix, A: Matrix, precision: Double = 1e-14): MGRepresentation {
    if (!checkMGRepresentation(alpha, A)) {
        throw IllegalArgumentException("DPHFromMG: Input isn't a valid MG distribution!")
    }

    // Transform function: apply similarity transformation B to representation
    val transfun: (List<Matrix>, Matrix) -> List<Matrix> = { rep, B ->
        val newAlpha = rep[0].mult(B)
        val newA = B.inv().mult(rep[1]).mult(B)
        listOf(newAlpha, newA)
    }

    // Evaluation function: measure distance from Markovian representation
    val evalfun: (List<Matrix>, Int) -> Double = { rep, k ->
        val ao = rep[0]
        val Ao = rep[1]
        val N = Ao.numRows

        // av = 1 - sum(Ao, 2) (closing vector)
        val av = DoubleArray(N)
        for (i in 0 until N) {
            var rowSum = 0.0
            for (j in 0 until N) {
                rowSum += Ao[i, j]
            }
            av[i] = 1.0 - rowSum
        }

        // Ad = Ao - diag(diag(Ao)) (off-diagonal elements)
        val Ad = Matrix.zeros(N, N)
        for (i in 0 until N) {
            for (j in 0 until N) {
                if (i != j) {
                    Ad[i, j] = Ao[i, j]
                }
            }
        }

        if (k % 2 == 0) {
            // Return max negative value
            var minVal = Double.MAX_VALUE
            for (i in 0 until ao.length()) {
                if (ao[i] < minVal) minVal = ao[i]
            }
            for (i in 0 until N) {
                if (av[i] < minVal) minVal = av[i]
            }
            for (i in 0 until N) {
                for (j in 0 until N) {
                    if (Ad[i, j] < minVal) minVal = Ad[i, j]
                }
            }
            -min(0.0, minVal)
        } else {
            // Return sum of negative values
            var negSum = 0.0
            for (i in 0 until ao.length()) {
                if (ao[i] < 0) negSum -= ao[i]
            }
            for (i in 0 until N) {
                if (av[i] < 0) negSum -= av[i]
            }
            for (i in 0 until N) {
                for (j in 0 until N) {
                    if (Ad[i, j] < 0) negSum -= Ad[i, j]
                }
            }
            negSum
        }
    }

    val rep = listOf(alpha, A)
    val nrep = findMarkovianRepresentation(rep, transfun, evalfun, precision)

    return MGRepresentation(nrep[0], nrep[1])
}

/**
 * Overload for DoubleArray alpha.
 */
fun dphFromMG(alpha: DoubleArray, A: Matrix, precision: Double = 1e-14): MGRepresentation {
    return dphFromMG(Matrix(alpha), A, precision)
}
