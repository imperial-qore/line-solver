/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 *
 * Reference:
 * Mocanu, S., Commault, C.: "Sparse representations of phase-type distributions,"
 * Stoch. Models 15, 759-778 (1999)
 */
package jline.lib.butools.reptrans

import jline.util.matrix.Matrix
import kotlin.math.abs

/**
 * Transforms an arbitrary matrix to a Markovian bi-diagonal matrix.
 *
 * @param A Matrix parameter of the initial representation (shape N,N)
 * @param maxSize The maximal order of the resulting Markovian representation (default 100)
 * @param precision Matrix entries smaller than the precision are considered to be zeros (default 1e-14)
 * @return Transient (bi-diagonal) generator matrix of the Markovian acyclic representation.
 *
 * Note: Calls the 'transformToMonocyclic' procedure if all the eigenvalues are real,
 * otherwise it raises an error if no Markovian acyclic generator has been found.
 *
 * @throws IllegalArgumentException if complex eigenvalues are found (no acyclic representation exists).
 */
fun transformToAcyclic(A: Matrix, maxSize: Int = 100, precision: Double = 1e-14): Matrix {
    // Check if any eigenvalue has non-zero imaginary part
    val eigenvalues = A.eig()
    for (ev in eigenvalues) {
        if (abs(ev.imaginary) >= precision) {
            throw IllegalArgumentException("TransformToAcyclic: Complex eigenvalue found, no acyclic representation exists.")
        }
    }

    return transformToMonocyclic(A, maxSize, precision)
}
