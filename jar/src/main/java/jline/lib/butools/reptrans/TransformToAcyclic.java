/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 *
 * Reference:
 * Mocanu, S., Commault, C.: "Sparse representations of phase-type distributions,"
 * Stoch. Models 15, 759-778 (1999)
 */
package jline.lib.butools.reptrans;

import java.util.List;

import org.apache.commons.math3.complex.Complex;

import jline.util.matrix.Matrix;

public final class TransformToAcyclic {
    private TransformToAcyclic() {}

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
    public static Matrix transformToAcyclic(Matrix A, int maxSize, double precision) {
        // Check if any eigenvalue has non-zero imaginary part
        List<Complex> eigenvalues = A.eig();
        for (Complex ev : eigenvalues) {
            if (Math.abs(ev.getImaginary()) >= precision) {
                throw new IllegalArgumentException("TransformToAcyclic: Complex eigenvalue found, no acyclic representation exists.");
            }
        }

        return TransformToMonocyclic.transformToMonocyclic(A, maxSize, precision);
    }

    public static Matrix transformToAcyclic(Matrix A, int maxSize) {
        return transformToAcyclic(A, maxSize, 1e-14);
    }

    public static Matrix transformToAcyclic(Matrix A) {
        return transformToAcyclic(A, 100, 1e-14);
    }
}
