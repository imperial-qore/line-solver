/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.util.matrix.Matrix
import java.util.Random

/**
 * Returns a random discrete Markovian arrival process.
 *
 * @param order The size of the DMAP
 * @param mean The mean inter-arrival time of the DMAP (default: 10.0)
 * @param zeroEntries The number of zero entries in the D0 and D1 matrices (default: 0)
 * @param maxTrials Maximum number of trials to find a proper DMAP (default: 1000)
 * @param prec Numerical precision for checking irreducibility (default: 1e-7)
 * @param random Random number generator
 * @return Pair of (D0, D1) matrices of the DMAP
 */
fun randomDMAP(
    order: Int,
    mean: Double = 10.0,
    zeroEntries: Int = 0,
    maxTrials: Int = 1000,
    prec: Double = 1e-7,
    random: Random = Random()
): Pair<Matrix, Matrix> {
    val D = randomDMMAP(order, 1, mean, zeroEntries, maxTrials, prec, random)
    return Pair(D[0], D[1])
}
