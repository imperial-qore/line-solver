/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.lib.butools.map.canonicalFromMAP2
import jline.lib.butools.mc.drpSolve
import jline.util.matrix.Matrix
import org.apache.commons.math3.complex.Complex
import kotlin.math.abs
import kotlin.math.sqrt

/**
 * Returns the canonical form of an order-2 discrete Markovian arrival process.
 *
 * @param D0 The D0 matrix of the DMAP(2)
 * @param D1 The D1 matrix of the DMAP(2)
 * @param prec Numerical precision to check the input
 * @return Pair of (G0, G1) matrices of the canonical DMAP(2)
 */
fun canonicalFromDMAP2(D0: Matrix, D1: Matrix, prec: Double = 1e-14): Pair<Matrix, Matrix> {
    if (!checkDMAPRepresentation(D0, D1, prec)) {
        throw IllegalArgumentException("CanonicalFromDMAP2: Input isn't a valid DMAP representation!")
    }

    if (D0.numRows != 2) {
        throw IllegalArgumentException("CanonicalFromDMAP2: Size is not 2!")
    }

    // Sort eigenvalues by descending real part
    val eigenvalues: List<Complex> = D0.eig()
    val sortedEig = eigenvalues.sortedByDescending { ev: Complex -> ev.real }
    val s1 = sortedEig[0].real
    val s2 = sortedEig[1].real

    if (s2 >= 0) {
        // Use continuous-time canonical form
        val I = Matrix.eye(2)
        val (G0, G1) = canonicalFromMAP2(D0.sub(I), D1, prec)
        return Pair(G0.add(I), G1)
    }

    // s2 is negative
    val I = Matrix.eye(2)
    val av = drpSolve(I.sub(D0).inv().mult(D1))

    // Get eigenvalue of transition matrix
    val P = I.sub(D0).inv().mult(D1)
    val gammaEig: List<Complex> = P.eig()
    val sortedGamma = gammaEig.sortedByDescending { ev: Complex -> abs(ev.real) }
    val gamma = sortedGamma[1].real

    // w1 = 1/(s1-s2)*(sum(D0,2)-s2*[1;1])
    val rowSum = Matrix(2, 1)
    rowSum[0, 0] = D0[0, 0] + D0[0, 1]
    rowSum[1, 0] = D0[1, 0] + D0[1, 1]

    val w1 = Matrix(2, 1)
    w1[0, 0] = (rowSum[0, 0] - s2) / (s1 - s2)
    w1[1, 0] = (rowSum[1, 0] - s2) / (s1 - s2)

    val w2 = Matrix(2, 1)
    w2[0, 0] = 1.0 - w1[0, 0]
    w2[1, 0] = 1.0 - w1[1, 0]

    // W = [w1 w2]
    val W = Matrix(2, 2)
    W[0, 0] = w1[0, 0]
    W[0, 1] = w2[0, 0]
    W[1, 0] = w1[1, 0]
    W[1, 1] = w2[1, 0]

    // A = (1-s1)*(av*W)
    val avW = av.mult(W)
    val A = avW.scale(1.0 - s1)
    val a1 = A[0, 0]

    val G0: Matrix
    val G1: Matrix

    if (gamma >= 0) {
        val temp1 = 1.0 + s1 * gamma - a1 * s1 + 5.0 * s1 * s1 - a1 * s1 * s1 - 2.0 * s1 * s1 * s1 - 2.0 * s2 - a1 * s2 + 5.0 * s1 * s2 - 3.0 * s1 * s1 * s2 + s2 * s2 + a1 * s2 * s2 - s1 * s2 * s2 - gamma + 3.0 * s1 * gamma - a1 * s1 * gamma - 3.0 * s1 * s1 * gamma + a1 * s1 * s1 * gamma + s1 * s1 * s1 * gamma + s2 * gamma + a1 * s2 * gamma - 2.0 * s1 * s2 * gamma + s1 * s1 * s2 * gamma - a1 * s2 * s2 * gamma
        val innerExpr1 = -1.0 + s1 * s1 * (-2.0 + gamma) + gamma + s2 * (1.0 + a1 - a1 * gamma) + s1 * (3.0 - a1 - s2 - 2.0 * gamma + a1 * gamma)
        val innerExpr2 = -s1 * s1 * s1 * (-1.0 + gamma) + a1 * (-1.0 + s2) * s2 * (-1.0 + gamma) + s1 * s1 * (-2.0 + a1 + s2 + 2.0 * gamma - a1 * gamma) + s1 * (1.0 - a1 - s2 - gamma + a1 * gamma)
        val sqrtExpr = (-1.0 + s1 + s2) * (-1.0 + s1 + s2) * (innerExpr1 * innerExpr1 - 4.0 * (-1.0 + s1) * innerExpr2)

        val a = -(1.0 / (2.0 * (-1.0 + s1) * (-1.0 + s1 + s2) * (-1.0 + s1 + s2))) * (1.0 - 4.0 * s1 + a1 * s1 + 5.0 * s1 * s1 - a1 * s1 * s1 - 2.0 * s1 * s1 * s1 - 2.0 * s2 - a1 * s2 + 5.0 * s1 * s2 - 3.0 * s1 * s1 * s2 + s2 * s2 + a1 * s2 * s2 - s1 * s2 * s2 - gamma + 3.0 * s1 * gamma - a1 * s1 * gamma - 3.0 * s1 * s1 * gamma + a1 * s1 * s1 * gamma + s1 * s1 * s1 * gamma + s2 * gamma + a1 * s2 * gamma - 2.0 * s1 * s2 * gamma + s1 * s1 * s2 * gamma - a1 * s2 * s2 * gamma + sqrt(sqrtExpr))
        val b = 1.0 + (a * (-1.0 + s1 + s2 - s1 * s2) * gamma) / ((a - 1.0) * (-s1 * s2 + a * (-1.0 + s1 + s2)))

        G0 = Matrix(2, 2)
        G0[0, 0] = s1 + s2
        G0[0, 1] = a * (1.0 - s1 - s2)
        G0[1, 0] = s1 * s2 / (a * (s1 + s2 - 1.0))
        G0[1, 1] = 0.0

        G1 = Matrix(2, 2)
        G1[0, 0] = (1.0 - a) * (1.0 - s1 - s2)
        G1[0, 1] = 0.0
        G1[1, 0] = b * (1.0 + s1 * s2 / (a * (1.0 - s1 - s2)))
        G1[1, 1] = (1.0 - b) * (1.0 + s1 * s2 / (a * (1.0 - s1 - s2)))
    } else {
        // gamma < 0
        val innerNum = a1 * s1 - a1 * s1 * s1 + s2 - a1 * s2 - 3.0 * s1 * s2 + 2.0 * s1 * s1 * s2 - s2 * s2 + a1 * s2 * s2 + s1 * s2 * s2 + s1 * gamma - a1 * s1 * gamma - 2.0 * s1 * s1 * gamma + a1 * s1 * s1 * gamma + s1 * s1 * s1 * gamma + a1 * s2 * gamma - a1 * s2 * s2 * gamma
        val sqrtArg = -4.0 * (-1.0 + s1) * s1 * s2 * (-1.0 + s1 + s2) * (a1 * (s1 - s2) * (-1.0 + gamma) + (-1.0 + s1) * (s2 + (-1.0 + s1) * gamma))
        val innerExpr = a1 * (-s1 + s1 * s1 + s2 - s2 * s2) * (-1.0 + gamma) + (-1.0 + s1) * ((-1.0 + 2.0 * s1) * s2 + s2 * s2 + (-1.0 + s1) * s1 * gamma)
        val sqrtFull = sqrtArg + innerExpr * innerExpr
        val denom = 2.0 * (-1.0 + s1 + s2) * (a1 * (s1 - s2) * (-1.0 + gamma) + (-1.0 + s1) * (s2 + (-1.0 + s1) * gamma))

        val a = (innerNum + sqrt(sqrtFull)) / denom
        val b = -((a * (1.0 - s1) * (1.0 - s2) * gamma) / ((a - 1.0) * (-a + a * s1 + a * s2 - s1 * s2)))

        G0 = Matrix(2, 2)
        G0[0, 0] = s1 + s2
        G0[0, 1] = a * (1.0 - s1 - s2)
        G0[1, 0] = s1 * s2 / (a * (s1 + s2 - 1.0))
        G0[1, 1] = 0.0

        G1 = Matrix(2, 2)
        G1[0, 0] = 0.0
        G1[0, 1] = (1.0 - a) * (1.0 - s1 - s2)
        G1[1, 0] = b * (1.0 - s1 * s2 / (a * (s1 + s2 - 1.0)))
        G1[1, 1] = (1.0 - b) * (1.0 - s1 * s2 / (a * (s1 + s2 - 1.0)))
    }

    return Pair(G0, G1)
}
