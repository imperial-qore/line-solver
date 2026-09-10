/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.mam;

import java.util.Map;

import jline.lib.butools.FluidFundamentalMatrices;
import jline.lib.butools.mc.CTMCSolve;
import jline.util.matrix.Matrix;

public final class FluidSolve {
    private FluidSolve() {}

    /**
     * Result class for FluidSolve.
     */
    public static final class FluidSolution {
        public final Matrix mass0;
        public final Matrix ini;
        public final Matrix K;
        public final Matrix clo;

        public FluidSolution(Matrix mass0, Matrix ini, Matrix K, Matrix clo) {
            this.mass0 = mass0;
            this.ini = ini;
            this.K = K;
            this.clo = clo;
        }
    }

    /**
     * Returns the parameters of the matrix-exponentially distributed stationary
     * distribution of a canonical Markovian fluid model.
     */
    public static FluidSolution fluidSolve(Matrix Fpp, Matrix Fpm, Matrix Fmp, Matrix Fmm, double prec) {
        int Np = Fpp.getNumRows();
        int Nm = Fmm.getNumRows();

        Map<String, Matrix> result = FluidFundamentalMatrices.FluidFundamentalMatrices(Fpp, Fpm, Fmp, Fmm, prec, null, null);
        Matrix Psi = result.get("P");
        Matrix K = result.get("K");
        Matrix U = result.get("U");

        Matrix mass0Minus = CTMCSolve.ctmcSolve(U, 1e-14);

        Matrix Ki = K.neg().inv();
        double nr = mass0Minus.elementSum() + 2 * mass0Minus.mult(Fmp).mult(Ki).elementSum();
        mass0Minus = mass0Minus.scale(1.0 / nr);

        Matrix ini = mass0Minus.mult(Fmp);

        Matrix clo = Matrix.zeros(Np, Np + Nm);
        for (int i = 0; i < Np; i++) {
            clo.set(i, i, 1.0);
            for (int j = 0; j < Nm; j++) {
                clo.set(i, Np + j, Psi.get(i, j));
            }
        }

        Matrix mass0 = Matrix.zeros(1, Np + Nm);
        for (int i = 0; i < Nm; i++) {
            mass0.set(0, Np + i, mass0Minus.get(0, i));
        }

        return new FluidSolution(mass0, ini, K, clo);
    }

    public static FluidSolution fluidSolve(Matrix Fpp, Matrix Fpm, Matrix Fmp, Matrix Fmm) {
        return fluidSolve(Fpp, Fpm, Fmp, Fmm, 1e-14);
    }
}
