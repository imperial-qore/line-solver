/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.pfqn.mva;

import jline.api.pfqn.nc.Pfqn_ncoi;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.function.ToDoubleFunction;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * Cross-validation of the OI mean-value routine against its normalizing-constant
 * sibling, plus the input guard and the population identity.
 *
 * pfqn_mvaoi and pfqn_ncoi compute the same exact per-class throughput by two
 * disjoint routes: the former from mean quantities only, the latter as the
 * normalizing-constant ratio G(N - e_r)/G(N). Agreement between them is a strong
 * check, since no intermediate is shared.
 */
public class Pfqn_mvaoiTest {

    private static final double TOL = 1e-12;

    /** A single-server OI station: unit total rate whenever it is nonempty. */
    private static ToDoubleFunction<int[]> unitRate() {
        return new ToDoubleFunction<int[]>() {
            @Override
            public double applyAsDouble(int[] n) {
                for (int v : n) if (v > 0) return 1.0;
                return 0.0;
            }
        };
    }

    /**
     * Reference values obtained from MATLAB pfqn_mvaoi for Z = [1 0.5], N = [2 1]
     * with one unit-rate OI station.
     */
    @Test
    public void twoClassSingleOiMatchesMatlabReference() {
        double[] Z = {1.0, 0.5};
        int[] N = {2, 1};
        Pfqn_mvaoi.Result r = Pfqn_mvaoi.pfqn_mvaoi(Z, N, unitRate());

        assertEquals(0.5925925926, r.X[0], 1e-9);
        assertEquals(0.3703703704, r.X[1], 1e-9);
        assertEquals(1.4074074074, r.Qoi[0][0], 1e-9);
        assertEquals(0.8148148148, r.Qoi[0][1], 1e-9);
        assertEquals(0.5925925926, r.Qdelay[0], 1e-9);
        assertEquals(0.1851851852, r.Qdelay[1], 1e-9);
        assertEquals(0.5925925926, r.getSoi()[0][0], 1e-9);
        assertEquals(0.3703703704, r.getSoi()[0][1], 1e-9);
    }

    /** Population conservation: with one OI station and a delay, Qoi + Qdelay = N. */
    @Test
    public void populationIsConserved() {
        double[] Z = {1.0, 0.5};
        int[] N = {2, 1};
        Pfqn_mvaoi.Result r = Pfqn_mvaoi.pfqn_mvaoi(Z, N, unitRate());
        for (int c = 0; c < N.length; c++) {
            assertEquals(N[c], r.Qoi[0][c] + r.Qdelay[c], TOL);
        }
    }

    /** X_r from mean-value analysis equals the normalizing-constant ratio G(N-e_r)/G(N). */
    @Test
    public void throughputMatchesNormalizingConstantRatio() {
        double[] Z = {1.0, 0.5};
        int[] N = {2, 1};
        List<ToDoubleFunction<int[]>> mu = new ArrayList<ToDoubleFunction<int[]>>();
        mu.add(unitRate());

        double lgN = Pfqn_ncoi.pfqn_ncoi(Z, N, mu).lG;
        for (int c = 0; c < N.length; c++) {
            int[] Nr = N.clone();
            Nr[c]--;
            double ratio = Math.exp(Pfqn_ncoi.pfqn_ncoi(Z, Nr, mu).lG - lgN);
            assertEquals(ratio, Pfqn_mvaoi.pfqn_mvaoi(Z, N, mu, null).X[c], 1e-10);
        }
    }

    /**
     * An empty mu list must be rejected, not silently degenerate to plain LI+delay
     * MVA: with K = 0 every OI loop is skipped and the caller gets a plausible
     * number for a network they did not describe.
     */
    @Test
    public void emptyRateListIsRejected() {
        assertThrows(IllegalArgumentException.class, new org.junit.jupiter.api.function.Executable() {
            @Override
            public void execute() {
                Pfqn_mvaoi.pfqn_mvaoi(new double[]{1.0, 0.5}, new int[]{2, 1},
                        new ArrayList<ToDoubleFunction<int[]>>(), null);
            }
        });
    }
}
