/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Checks if the input matrices define a discrete time MRAP.
 *
 * All matrices H0...HK must have the same size, the dominant eigenvalue
 * of H0 is real and less than 1, and the rowsum of H0+H1+...+HK is 1
 * (up to the numerical precision).
 *
 * @param H The H0...HK matrices of the DMRAP to check (as MatrixCell)
 * @param prec Numerical precision, default value is 1e-14
 * @return True if the matrices define a valid DMRAP
 */
fun checkDMRAPRepresentation(H: MatrixCell, prec: Double = 1e-14): Boolean {
    if (H.size() < 2) {
        return false
    }

    // Sum H1...HK
    var sumH = H[1].copy()
    for (i in 2 until H.size()) {
        sumH = sumH.add(H[i])
    }

    return checkDRAPRepresentation(H[0], sumH, prec)
}

/**
 * Overload for Array<Matrix>.
 */
fun checkDMRAPRepresentation(H: Array<Matrix>, prec: Double = 1e-14): Boolean {
    if (H.size < 2) {
        return false
    }

    // Sum H1...HK
    var sumH = H[1].copy()
    for (i in 2 until H.size) {
        sumH = sumH.add(H[i])
    }

    return checkDRAPRepresentation(H[0], sumH, prec)
}
