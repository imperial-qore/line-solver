/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.butools.mam;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Result class for GeneralFluidSolve containing the parameters of the
 * matrix-exponentially distributed stationary distribution.
 */
public final class GeneralFluidSolution {
    private final Matrix mass0;
    private final Matrix ini;
    private final Matrix K;
    private final Matrix clo;

    public GeneralFluidSolution(Matrix mass0, Matrix ini, Matrix K, Matrix clo) {
        this.mass0 = mass0;
        this.ini = ini;
        this.K = K;
        this.clo = clo;
    }

    public Matrix getMass0() { return mass0; }
    public Matrix getIni() { return ini; }
    public Matrix getK() { return K; }
    public Matrix getClo() { return clo; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof GeneralFluidSolution)) return false;
        GeneralFluidSolution that = (GeneralFluidSolution) o;
        return Objects.equals(mass0, that.mass0)
                && Objects.equals(ini, that.ini)
                && Objects.equals(K, that.K)
                && Objects.equals(clo, that.clo);
    }

    @Override
    public int hashCode() {
        return Objects.hash(mass0, ini, K, clo);
    }

    @Override
    public String toString() {
        return "GeneralFluidSolution(mass0=" + mass0 + ", ini=" + ini
                + ", K=" + K + ", clo=" + clo + ")";
    }
}
