/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.pfqn.lcfs;

import jline.util.matrix.Matrix;

/**
 * Result class for LCFS MVA algorithm.
 */
public final class LcfsqnMvaResult {
    private final Matrix T;
    private final Matrix Q;
    private final Matrix U;
    private final Matrix B;

    public LcfsqnMvaResult(Matrix T, Matrix Q, Matrix U, Matrix B) {
        this.T = T;
        this.Q = Q;
        this.U = U;
        this.B = B;
    }

    public Matrix getT() { return T; }
    public Matrix getQ() { return Q; }
    public Matrix getU() { return U; }
    public Matrix getB() { return B; }
}
