/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mam.handlers;

import java.util.Arrays;
import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Container for transient QBD analysis results.
 */
public final class TransientResult {
    private final Matrix[][] Qt;
    private final Matrix[][] Ut;
    private final Matrix[][] Tt;

    public TransientResult(Matrix[][] Qt, Matrix[][] Ut, Matrix[][] Tt) {
        this.Qt = Qt;
        this.Ut = Ut;
        this.Tt = Tt;
    }

    public Matrix[][] getQt() { return Qt; }
    public Matrix[][] getUt() { return Ut; }
    public Matrix[][] getTt() { return Tt; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof TransientResult)) return false;
        TransientResult that = (TransientResult) o;
        return Arrays.deepEquals(Qt, that.Qt)
                && Arrays.deepEquals(Ut, that.Ut)
                && Arrays.deepEquals(Tt, that.Tt);
    }

    @Override
    public int hashCode() {
        return Objects.hash(Arrays.deepHashCode(Qt), Arrays.deepHashCode(Ut), Arrays.deepHashCode(Tt));
    }
}
