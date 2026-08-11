/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Result class for PH3From5Moments.
 */
public final class PH3Representation {
    public final Matrix alpha;
    public final Matrix A;

    public PH3Representation(Matrix alpha, Matrix A) {
        this.alpha = alpha;
        this.A = A;
    }

    public Matrix getAlpha() { return alpha; }
    public Matrix getA() { return A; }

    public Matrix component1() { return alpha; }
    public Matrix component2() { return A; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof PH3Representation)) return false;
        PH3Representation that = (PH3Representation) o;
        return Objects.equals(alpha, that.alpha) && Objects.equals(A, that.A);
    }

    @Override
    public int hashCode() {
        return Objects.hash(alpha, A);
    }

    @Override
    public String toString() {
        return "PH3Representation(alpha=" + alpha + ", A=" + A + ")";
    }
}
