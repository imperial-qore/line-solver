/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Obtains a Markovian representation of a discrete rational
 * arrival process of the same size, if possible.
 *
 * @param H0 The H0 matrix of the discrete rational arrival process
 * @param H1 The H1 matrix of the discrete rational arrival process
 * @param prec A representation is considered to be Markovian if it is closer than this precision
 * @return Pair of (D0, D1) matrices of the discrete Markovian arrival process
 */
fun dmapFromDRAP(H0: Matrix, H1: Matrix, prec: Double = 1e-14): Pair<Matrix, Matrix> {
    if (!checkDRAPRepresentation(H0, H1, prec)) {
        throw IllegalArgumentException("DMAPFromDRAP: Input isn't a valid DRAP representation!")
    }

    val H = MatrixCell(2)
    H[0] = H0
    H[1] = H1

    val Y = dmmapFromDMRAP(H, prec)

    return Pair(Y[0], Y[1])
}
