/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.pfqn.nc;

import jline.io.Ret;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.ToDoubleFunction;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Macrostate OI normalizing constant (Pfqn_ncoi) against the microstate
 * pass-and-swap walk (Pfqn_pas_nc) and against MATLAB fixtures.
 *
 * The two routines must agree to machine precision on every OI model: the
 * macrostate convolution is legitimate exactly because an OI rank rate is
 * permutation-invariant. The PAS fixture checks the placement-order pruning,
 * whose reference is the product-form sum over the communicating class
 * enumerated by SolverCTMC in MATLAB.
 */
public class Pfqn_ncoiTest {

    private static final double TOL = 1e-12;

    /** OI rank rate: sum of the capacities of the classes present in n. */
    private static ToDoubleFunction<int[]> rank(final double[] cap) {
        return new ToDoubleFunction<int[]>() {
            @Override
            public double applyAsDouble(int[] n) {
                double s = 0.0;
                for (int r = 0; r < n.length; r++) {
                    if (n[r] > 0) {
                        s += cap[r];
                    }
                }
                return s;
            }
        };
    }

    private static ToDoubleFunction<int[]> constant(final double c) {
        return new ToDoubleFunction<int[]>() {
            @Override
            public double applyAsDouble(int[] n) {
                return c;
            }
        };
    }

    private static void assertMacroEqualsMicro(double[] Z, int[] N, List<ToDoubleFunction<int[]>> mu) {
        Ret.pfqnOiNc macro = Pfqn_ncoi.pfqn_ncoi(Z, N, mu);
        Ret.pfqnOiNc micro = Pfqn_pas_nc.pfqn_pas_nc(Z, N, mu);
        assertEquals(micro.G, macro.G, TOL * Math.max(1.0, Math.abs(micro.G)));
        assertEquals(micro.lG, macro.lG, TOL * Math.max(1.0, Math.abs(micro.lG)));
    }

    @Test
    public void macrostateMatchesMicrostateOnOiModels() {
        List<ToDoubleFunction<int[]>> one = new ArrayList<ToDoubleFunction<int[]>>();
        one.add(rank(new double[]{0.8, 1.2}));
        assertMacroEqualsMicro(new double[]{0.7, 1.3}, new int[]{3, 2}, one);

        List<ToDoubleFunction<int[]>> two = new ArrayList<ToDoubleFunction<int[]>>();
        two.add(rank(new double[]{0.8, 1.2}));
        two.add(rank(new double[]{1.5, 0.6}));
        assertMacroEqualsMicro(new double[]{1.0, 0.5}, new int[]{2, 3}, two);

        List<ToDoubleFunction<int[]>> three = new ArrayList<ToDoubleFunction<int[]>>();
        three.add(rank(new double[]{0.9, 1.1, 0.4}));
        three.add(rank(new double[]{0.5, 0.7, 1.3}));
        assertMacroEqualsMicro(new double[]{0.4, 0.9, 1.1}, new int[]{2, 1, 2}, three);

        // No delay demand: only the empty delay state is feasible.
        List<ToDoubleFunction<int[]>> nodelay = new ArrayList<ToDoubleFunction<int[]>>();
        nodelay.add(rank(new double[]{1.0, 2.0}));
        assertMacroEqualsMicro(new double[]{0.0, 0.0}, new int[]{4, 4}, nodelay);
    }

    @Test
    public void constantRateReducesToTheLoadIndependentConstant() {
        // Two OI stations with unit rank rate and a delay: MATLAB pfqn_ca fixture.
        List<ToDoubleFunction<int[]>> mu = new ArrayList<ToDoubleFunction<int[]>>();
        mu.add(constant(1.0));
        mu.add(constant(1.0));
        Ret.pfqnOiNc r = Pfqn_ncoi.pfqn_ncoi(new double[]{0.3, 0.7}, new int[]{2, 3}, mu);
        assertEquals(93.4059225, r.G, 1e-7);
    }

    @Test
    public void placementOrderSelectsTheCommunicatingClass() {
        // Four single-job classes on a two-station P&S cycle with a cyclic swap
        // graph. The placement order is total, so the communicating class holds
        // one ordering per split and G_C = sum_k (1/mu1)^k (1/mu2)^(4-k).
        // MATLAB reference (SolverCTMC state space): 3.166240678.
        int R = 4;
        boolean[][] P = new boolean[R][R];
        for (int i = 0; i < R; i++) {
            for (int j = i + 1; j < R; j++) {
                P[i][j] = true;
            }
        }
        boolean[][] Pt = new boolean[R][R];
        for (int i = 0; i < R; i++) {
            for (int j = 0; j < R; j++) {
                Pt[i][j] = P[j][i];
            }
        }
        List<ToDoubleFunction<int[]>> mu = new ArrayList<ToDoubleFunction<int[]>>();
        mu.add(constant(1.0));
        mu.add(constant(1.3));
        Ret.pfqnOiNc r = Pfqn_pas_nc.pfqn_pas_nc(null, new int[]{1, 1, 1, 1}, mu,
                Arrays.asList(P, Pt));
        assertEquals(3.166240678, r.G, 1e-8);

        // Without the placement order the same model is plain OI, which counts
        // every ordering and is therefore strictly larger.
        Ret.pfqnOiNc oi = Pfqn_ncoi.pfqn_ncoi(null, new int[]{1, 1, 1, 1}, mu);
        assertEquals(75.98977627, oi.G, 1e-6);
    }
}
