/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.butools.reptrans;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Result class for ExtendToMarkovian.
 */
public final class MarkovianRepresentation {
    private final Matrix beta;
    private final Matrix B;

    public MarkovianRepresentation(Matrix beta, Matrix B) {
        this.beta = beta;
        this.B = B;
    }

    public Matrix getBeta() { return beta; }
    public Matrix getB() { return B; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof MarkovianRepresentation)) return false;
        MarkovianRepresentation that = (MarkovianRepresentation) o;
        return Objects.equals(beta, that.beta) && Objects.equals(B, that.B);
    }

    @Override
    public int hashCode() {
        return Objects.hash(beta, B);
    }

    @Override
    public String toString() {
        return "MarkovianRepresentation(beta=" + beta + ", B=" + B + ")";
    }
}
