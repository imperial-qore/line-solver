/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.pfqn.ld;

import jline.io.Ret;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Load-dependent normalizing constant dispatcher.
 *
 * The single-queueing-station-with-delay case is the canonical LD repairman
 * model and had no test, which is how a wrong rate matrix (the delay-augmented
 * muz where the unaugmented mu is required) reached that branch and returned a
 * silently wrong lG. The fixtures below are MATLAB values.
 */
public class Pfqn_ncldTest {

    private static final double TOL = 1e-9;

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) m.set(0, i, v[i]);
        return m;
    }

    /**
     * M==1 with a delay: dispatches to exact/comomld. MATLAB pfqn_ncld returns
     * 0.17227122094 here, which pfqn_comomrm_ld(L,N,Z,mu) reproduces exactly
     * and pfqn_comomrm_ld(L,N,Z,muz) does not (it gives 0.486123011126).
     */
    @Test
    public void ldRepairmanWithDelay() {
        SolverOptions options = new SolverOptions();
        Ret.pfqnNc ret = Pfqn_ncld.pfqn_ncld(row(0.6, 0.4), row(2, 1), row(1, 0.5),
                row(1, 2, 2), options);
        assertEquals(0.17227122094, ret.lG, TOL);
    }

    /**
     * Same demands and population with a unit rate vector, where the model
     * degenerates to the constant-rate repairman: lG must fall back to the
     * value the whole comom family agrees on.
     */
    @Test
    public void unitRatesMatchConstantRateRepairman() {
        SolverOptions options = new SolverOptions();
        Ret.pfqnNc ret = Pfqn_ncld.pfqn_ncld(row(0.6, 0.4), row(2, 1), row(1, 0.5),
                row(1, 1, 1), options);
        assertEquals(0.610851937833, ret.lG, TOL);
    }
}
