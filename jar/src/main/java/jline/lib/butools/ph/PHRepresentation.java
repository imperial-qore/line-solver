/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Result class for PH representations containing alpha (initial vector) and A (generator matrix).
 */
public final class PHRepresentation {
    public final Matrix alpha;
    public final Matrix A;

    public PHRepresentation(Matrix alpha, Matrix A) {
        this.alpha = alpha;
        this.A = A;
    }

    public Matrix getAlpha() {
        return alpha;
    }

    public Matrix getA() {
        return A;
    }

    public Matrix component1() {
        return alpha;
    }

    public Matrix component2() {
        return A;
    }

    public PHRepresentation copy(Matrix alpha, Matrix A) {
        return new PHRepresentation(alpha, A);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof PHRepresentation)) return false;
        PHRepresentation that = (PHRepresentation) o;
        return Objects.equals(alpha, that.alpha) && Objects.equals(A, that.A);
    }

    @Override
    public int hashCode() {
        return Objects.hash(alpha, A);
    }

    @Override
    public String toString() {
        return "PHRepresentation(alpha=" + alpha + ", A=" + A + ")";
    }
}
