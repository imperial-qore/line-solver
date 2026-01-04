/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Checks if the input matrices define a discrete time MMAP.
 *
 * All matrices D0...DK must have the same size, D0 must be a transient
 * probability matrix, D1...DK have only non-negative elements, and the
 * rowsum of D0+D1+...+DK is 1 (up to the numerical precision).
 *
 * @param D The D0...DK matrices of the DMMAP to check (as MatrixCell)
 * @param prec Numerical precision, default value is 1e-14
 * @return True if the matrices define a valid DMMAP
 */
fun checkDMMAPRepresentation(D: MatrixCell, prec: Double = 1e-14): Boolean {
    if (D.size() < 2) {
        return false
    }

    // Check if any D1...DK has negative elements
    for (i in 1 until D.size()) {
        if (D[i].elementMin() < -prec) {
            return false
        }
    }

    // Sum D1...DK
    var sumD = D[1].copy()
    for (i in 2 until D.size()) {
        sumD = sumD.add(D[i])
    }

    return checkDMAPRepresentation(D[0], sumD, prec)
}

/**
 * Overload for Array<Matrix>.
 */
fun checkDMMAPRepresentation(D: Array<Matrix>, prec: Double = 1e-14): Boolean {
    if (D.size < 2) {
        return false
    }

    // Check if any D1...DK has negative elements
    for (i in 1 until D.size) {
        if (D[i].elementMin() < -prec) {
            return false
        }
    }

    // Sum D1...DK
    var sumD = D[1].copy()
    for (i in 2 until D.size) {
        sumD = sumD.add(D[i])
    }

    return checkDMAPRepresentation(D[0], sumD, prec)
}
