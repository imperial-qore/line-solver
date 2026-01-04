/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.io.Serializable;

/**
 * A Marked Markovian Arrival Process (MMAP)
 *
 * MMAP is a point process where arrivals are categorized into K different types (marks).
 * Uses the standard MMAP representation:
 * - D0: infinitesimal generator for transitions without arrivals
 * - D1: rate matrix for transitions generating type-1 arrivals
 * - D2: rate matrix for transitions generating type-2 arrivals
 * - ...
 * - Dk: rate matrix for transitions generating type-k arrivals
 *
 * This class is an alias for MarkedMAP, providing a consistent naming convention
 * with other distribution classes in LINE.
 */
public class MMAP extends MarkedMAP implements Serializable {

    /**
     * Default constructor creating an empty MMAP
     */
    public MMAP() {
        super();
    }

    /**
     * Construct an MMAP from a MatrixCell containing {D0, D1, D2, ..., Dk}
     *
     * @param mmap MatrixCell where:
     *             - mmap.get(0) = D0 (transitions without arrivals)
     *             - mmap.get(k) = Dk (transitions generating type-k arrivals) for k >= 1
     */
    public MMAP(MatrixCell mmap) {
        super(mmap);
    }

    /**
     * Construct an MMAP from D0 and variable number of Dk matrices
     *
     * @param D0 the infinitesimal generator for transitions without arrivals
     * @param Dk variable number of matrices where Dk[i-1] represents transitions generating type-i arrivals
     */
    public MMAP(Matrix D0, Matrix... Dk) {
        super(createMatrixCell(D0, Dk));
    }

    /**
     * Returns the name of this distribution type.
     *
     * @return "MMAP" to identify this as a Marked Markovian Arrival Process
     */
    @Override
    public String getName() {
        return "MMAP";
    }

    /**
     * Helper method to create MatrixCell from D0 and Dk matrices
     */
    private static MatrixCell createMatrixCell(Matrix D0, Matrix... Dk) {
        MatrixCell mmap = new MatrixCell();
        mmap.set(0, D0);
        for (int k = 0; k < Dk.length; k++) {
            mmap.set(k + 1, Dk[k]);
        }
        return mmap;
    }

    /**
     * Get the number of arrival types (marks) in this MMAP
     *
     * @return the number of arrival types K
     */
    public int getNumTypes() {
        return getProcess().size() - 1;
    }

    /**
     * Get the D_k matrix for arrival type k
     *
     * @param type the arrival type (1, 2, 3, ...)
     * @return D_k matrix for the specified type, or null if out of range
     */
    public Matrix getTypeMatrix(int type) {
        if (type < 1 || type > getNumTypes()) {
            return null;
        }
        return getProcess().get(type);
    }

    /**
     * Factory method to create a simple 2-type MMAP from two MAPs
     *
     * @param D0 the shared D0 matrix (transitions without arrivals)
     * @param D1 rate matrix for type-1 arrivals
     * @param D2 rate matrix for type-2 arrivals
     * @return an MMAP with two arrival types
     */
    public static MMAP fromMAPs(Matrix D0, Matrix D1, Matrix D2) {
        return new MMAP(D0, D1, D2);
    }

    @Override
    public String toString() {
        return "MMAP(phases=" + nPhases + ", types=" + getNumTypes() + ")";
    }
}
