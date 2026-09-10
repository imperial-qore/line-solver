/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.mc;

import jline.util.matrix.Matrix;

/** Per-source and initial-law-weighted passage-time moments. */
public class PassageMomentsResult {
    /**
     * (nstates x nmax): row i holds the moments for a passage STARTED IN STATE
     * i, zero on target states and infinite where the target cannot be reached.
     */
    public final Matrix mall;
    /** The pi0-weighted moment vector, length nmax. */
    public final double[] m;

    public PassageMomentsResult(Matrix mall, double[] m) {
        this.mall = mall;
        this.m = m;
    }
}
