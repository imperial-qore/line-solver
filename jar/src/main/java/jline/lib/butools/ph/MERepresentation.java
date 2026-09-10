/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.butools.ph;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Result class for MEFromMoments containing both alpha and A.
 */
public final class MERepresentation {
    public final Matrix alpha;
    public final Matrix A;

    public MERepresentation(Matrix alpha, Matrix A) {
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
        if (!(o instanceof MERepresentation)) return false;
        MERepresentation that = (MERepresentation) o;
        return Objects.equals(alpha, that.alpha) && Objects.equals(A, that.A);
    }

    @Override
    public int hashCode() {
        return Objects.hash(alpha, A);
    }

    @Override
    public String toString() {
        return "MERepresentation(alpha=" + alpha + ", A=" + A + ")";
    }
}
