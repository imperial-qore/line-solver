/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import jline.api.npfqn.Npfqn_bnd_bgt;
import jline.api.npfqn.Npfqn_bnd_bpt;
import jline.util.matrix.Matrix;

/**
 * Cross-language regression for the two Bertsimas open-network bounds.
 *
 * <p>{@code npfqn_bnd_bpt} is the first-order LP relaxation of the achievable
 * region (Bertsimas-Paschalidis-Tsitsiklis, Ann. Appl. Prob. 4(1), 1994) and
 * {@code npfqn_bnd_bgt} is the piecewise-linear Lyapunov bound
 * (Bertsimas-Gamarnik-Tsitsiklis, Ann. Appl. Prob. 11(4), 2001). The reference
 * values are MATLAB's (matlab/src/api/npfqn/) and are reproduced by the python
 * and C++ ports to the digits asserted here; see _kb/06-solver-catalog.md under
 * {@code bpt.lower} and {@code bgt.upper}.</p>
 */
public class NpfqnBndAchievableTest {

    private static double bpt(double[] l0, double[] mu, Matrix P, int[] st, double[] c) {
        return Npfqn_bnd_bpt.npfqn_bnd_bpt(l0, mu, P, st, c).zlb;
    }

    /** The achievable-region relaxation is tight on M/M/1 at every load. */
    @Test
    public void bptIsExactOnMM1() {
        double[] rhos = {0.3, 0.5, 0.7, 0.9};
        for (int i = 0; i < rhos.length; i++) {
            double z = bpt(new double[]{rhos[i]}, new double[]{1.0}, new Matrix(1, 1),
                    new int[]{0}, new double[]{1.0});
            assertEquals(1.0 / (1.0 - rhos[i]), z, 1e-9);
        }
    }

    /** The bound sits below the cmu-optimal achievable point, 3.6875. */
    @Test
    public void bptTwoClassesAtOneStation() {
        Matrix P = new Matrix(2, 2);
        double[] l0 = {0.4, 0.4}, mu = {2.0, 1.0};
        int[] st = {0, 0};
        assertEquals(3.4375, bpt(l0, mu, P, st, new double[]{1, 1}), 1e-9);
        assertEquals(0.625, bpt(l0, mu, P, st, new double[]{1, 0}), 1e-9);
        assertEquals(1.0 / 0.6, bpt(l0, mu, P, st, new double[]{0, 1}), 1e-9);
        assertTrue(bpt(l0, mu, P, st, new double[]{1, 1}) <= 3.6875 + 1e-12);
    }

    /** Station 1 is an M/M/1 in isolation; station 2 falls back to 1/mu. */
    @Test
    public void bptTandemFirstStationIsExact() {
        Matrix P = new Matrix(2, 2);
        P.set(0, 1, 1.0);
        double[] l0 = {0.5, 0.0}, mu = {1.0, 1.0};
        int[] st = {0, 1};
        assertEquals(2.0, bpt(l0, mu, P, st, new double[]{1, 0}), 1e-9);
        assertEquals(1.0, bpt(l0, mu, P, st, new double[]{0, 1}), 1e-9);
    }

    @Test
    public void bptRefusesASaturatedStation() {
        assertThrows(RuntimeException.class, () -> bpt(new double[]{1.2}, new double[]{1.0},
                new Matrix(1, 1), new int[]{0}, new double[]{1.0}));
    }

    /** With J = 1 the uniformized drift of L = 1 gives (mu-lambda)/(lambda+mu). */
    @Test
    public void bgtMM1GammaHasAClosedForm() {
        double[] rhos = {0.3, 0.5, 0.7, 0.9};
        double[] qref = {23.9340659340659, 32.6666666666667, 53.6862745098039, 160.105263157895};
        for (int i = 0; i < rhos.length; i++) {
            Npfqn_bnd_bgt.Result r = Npfqn_bnd_bgt.npfqn_bnd_bgt(new double[]{rhos[i]},
                    new double[][]{{1.0}}, new int[][]{{0}}, 1);
            assertEquals((1 - rhos[i]) / (1 + rhos[i]), r.gamma, 1e-9);
            assertEquals(qref[i], r.Qub[0][0], 1e-7);
            // it IS an upper bound on the exact mean queue length
            assertTrue(r.Qub[0][0] >= rhos[i] / (1 - rhos[i]));
        }
    }

    @Test
    public void bgtTandem() {
        Npfqn_bnd_bgt.Result r = Npfqn_bnd_bgt.npfqn_bnd_bgt(new double[]{0.5},
                new double[][]{{1.0, 1.0}}, new int[][]{{0, 1}}, 2);
        assertEquals(0.2, r.gamma, 1e-9);
        assertEquals(1.0, r.Lmax, 1e-9);
        assertEquals(5529.6, r.B, 1e-6);
        assertEquals(5578.0, r.U, 1e-6);
        assertEquals(1.1 / 1.15, r.tailRatio, 1e-9);
    }

    @Test
    public void bgtStableLuKumar() {
        Npfqn_bnd_bgt.Result r = Npfqn_bnd_bgt.npfqn_bnd_bgt(new double[]{1.0},
                new double[][]{{1 / 0.3, 1 / 0.6, 1 / 0.3, 1 / 0.1}},
                new int[][]{{0, 1, 1, 0}}, 2);
        assertEquals(0.00574712643678161, r.gamma, 1e-12);
        assertEquals(0.4, r.rhoStation[0], 1e-12);
        assertEquals(0.9, r.rhoStation[1], 1e-12);
    }

    /**
     * Every station is at 0.7, yet rho_2 + rho_4 = 1.2 &gt; 1 (Rybko-Stolyar), so
     * the network is globally unstable. A per-station load test would accept
     * it; GLP[dm] must not.
     */
    @Test
    public void bgtRefusesAGloballyUnstableLuKumar() {
        assertThrows(RuntimeException.class, () -> Npfqn_bnd_bgt.npfqn_bnd_bgt(new double[]{1.0},
                new double[][]{{1 / 0.1, 1 / 0.6, 1 / 0.1, 1 / 0.6}},
                new int[][]{{0, 1, 1, 0}}, 2));
    }
}
