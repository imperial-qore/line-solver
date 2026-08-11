/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.lib.fjcodes;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Service representation for 2-node FJ job.
 */
public final class FJServiceH {
    private final Matrix service_phases;  // Possible phase combinations (rows are states)
    private final Matrix beta;            // Initial probability for 2-node job (row)
    private final Matrix S;               // PH sub-generator for 2-node job

    public FJServiceH(Matrix service_phases, Matrix beta, Matrix S) {
        this.service_phases = service_phases;
        this.beta = beta;
        this.S = S;
    }

    public Matrix getService_phases() { return service_phases; }
    public Matrix getBeta() { return beta; }
    public Matrix getS() { return S; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof FJServiceH)) return false;
        FJServiceH that = (FJServiceH) o;
        return Objects.equals(service_phases, that.service_phases)
                && Objects.equals(beta, that.beta)
                && Objects.equals(S, that.S);
    }

    @Override
    public int hashCode() {
        return Objects.hash(service_phases, beta, S);
    }

    @Override
    public String toString() {
        return "FJServiceH(service_phases=" + service_phases + ", beta=" + beta + ", S=" + S + ")";
    }
}
