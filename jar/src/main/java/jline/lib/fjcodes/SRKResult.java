/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.lib.fjcodes;

import jline.util.matrix.Matrix;

/**
 * Result of constructSRK.
 */
public final class SRKResult {
    public final Matrix Se;
    public final Matrix Sestar;
    public final Matrix R0;
    public final Matrix Ke;
    public final Matrix Kc;
    public final Matrix S_long;
    public final Matrix S_last;

    public SRKResult(Matrix Se, Matrix Sestar, Matrix R0, Matrix Ke, Matrix Kc,
                     Matrix S_long, Matrix S_last) {
        this.Se = Se;
        this.Sestar = Sestar;
        this.R0 = R0;
        this.Ke = Ke;
        this.Kc = Kc;
        this.S_long = S_long;
        this.S_last = S_last;
    }

    public Matrix getSe() { return Se; }
    public Matrix getSestar() { return Sestar; }
    public Matrix getR0() { return R0; }
    public Matrix getKe() { return Ke; }
    public Matrix getKc() { return Kc; }
    public Matrix getS_long() { return S_long; }
    public Matrix getS_last() { return S_last; }
}
