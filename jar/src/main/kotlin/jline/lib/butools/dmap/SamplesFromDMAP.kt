/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import java.util.Random

/**
 * Generates random samples from a discrete Markovian arrival process.
 *
 * @param D0 The D0 matrix of the discrete MAP
 * @param D1 The D1 matrix of the discrete MAP
 * @param K The number of samples to generate
 * @param initial Optional initial state (1-indexed)
 * @param prec Numerical precision for validation
 * @param random Random number generator
 * @return Vector of random samples (inter-arrival times)
 */
fun samplesFromDMAP(
    D0: Matrix,
    D1: Matrix,
    K: Int,
    initial: Int? = null,
    prec: Double = 1e-14,
    random: Random = Random()
): IntArray {
    if (!checkDMAPRepresentation(D0, D1, prec)) {
        throw IllegalArgumentException("SamplesFromDMAP: Input isn't a valid DMAP representation!")
    }

    val D = MatrixCell(2)
    D[0] = D0
    D[1] = D1

    return samplesFromDMMAP(D, K, initial, prec, random) as IntArray
}
