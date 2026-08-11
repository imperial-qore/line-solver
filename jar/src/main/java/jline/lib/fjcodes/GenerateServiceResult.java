/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.lib.fjcodes;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Result of generateService.
 */
public final class GenerateServiceResult {
    private final Matrix T;
    private final int newdim;
    private final int dim_notbusy;

    public GenerateServiceResult(Matrix T, int newdim, int dim_notbusy) {
        this.T = T;
        this.newdim = newdim;
        this.dim_notbusy = dim_notbusy;
    }

    public Matrix getT() { return T; }
    public int getNewdim() { return newdim; }
    public int getDim_notbusy() { return dim_notbusy; }

    public Matrix component1() { return T; }
    public int component2() { return newdim; }
    public int component3() { return dim_notbusy; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof GenerateServiceResult)) return false;
        GenerateServiceResult that = (GenerateServiceResult) o;
        return newdim == that.newdim && dim_notbusy == that.dim_notbusy && Objects.equals(T, that.T);
    }

    @Override
    public int hashCode() {
        return Objects.hash(T, newdim, dim_notbusy);
    }

    @Override
    public String toString() {
        return "GenerateServiceResult(T=" + T + ", newdim=" + newdim + ", dim_notbusy=" + dim_notbusy + ")";
    }
}
