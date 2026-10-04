/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.sections;

import jline.util.matrix.Matrix;

import java.io.Serializable;

/**
 * Input class for ClassSwitcher function.
 *
 * This class is the input parameter for the class switching function (csFun).
 * It contains the row index (r), column index (s), and optional state information.
 *
 * Note: This was previously an inner class of ClassSwitcher but has been moved
 * to a top-level class to avoid classloader issues in MATLAB's Java integration.
 */
public class CSFunInput implements Serializable {
    private static final long serialVersionUID = 1L;

    /** Row index (source class) */
    public int r;

    /** Column index (destination class) */
    public int s;

    /** Optional state matrix */
    public Matrix state;

    /** Optional state-dependent matrix */
    public Matrix statedep;

    /**
     * Creates a new CSFunInput.
     *
     * @param r Row index (source class)
     * @param s Column index (destination class)
     * @param state Optional state matrix
     * @param statedep Optional state-dependent matrix
     */
    public CSFunInput(int r, int s, Matrix state, Matrix statedep) {
        this.r = r;
        this.s = s;
        this.state = state;
        this.statedep = statedep;
    }
}
