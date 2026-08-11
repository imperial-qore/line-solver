/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.lib.fjcodes;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Result of build_SA.
 */
public final class SAResult {
    private final Matrix S;
    private final Matrix A_jump;

    public SAResult(Matrix S, Matrix A_jump) {
        this.S = S;
        this.A_jump = A_jump;
    }

    public Matrix getS() { return S; }
    public Matrix getA_jump() { return A_jump; }

    public Matrix component1() { return S; }
    public Matrix component2() { return A_jump; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof SAResult)) return false;
        SAResult that = (SAResult) o;
        return Objects.equals(S, that.S) && Objects.equals(A_jump, that.A_jump);
    }

    @Override
    public int hashCode() {
        return Objects.hash(S, A_jump);
    }

    @Override
    public String toString() {
        return "SAResult(S=" + S + ", A_jump=" + A_jump + ")";
    }
}
