/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Result class for RandomMAP.
 */
public final class MAPRepresentation {
    private final Matrix D0;
    private final Matrix D1;

    public MAPRepresentation(Matrix D0, Matrix D1) {
        this.D0 = D0;
        this.D1 = D1;
    }

    public Matrix getD0() { return D0; }
    public Matrix getD1() { return D1; }

    public Matrix component1() { return D0; }
    public Matrix component2() { return D1; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof MAPRepresentation)) return false;
        MAPRepresentation that = (MAPRepresentation) o;
        return Objects.equals(D0, that.D0) && Objects.equals(D1, that.D1);
    }

    @Override
    public int hashCode() {
        return Objects.hash(D0, D1);
    }

    @Override
    public String toString() {
        return "MAPRepresentation(D0=" + D0 + ", D1=" + D1 + ")";
    }
}
