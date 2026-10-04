/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.examples.java.advanced;

import jline.api.mc.Ctmc_hitting_time;
import jline.api.mc.Ctmc_makeinfgen;
import jline.api.mc.Ctmc_passage_moments;
import jline.api.mc.Ctmc_passage_time;
import jline.api.mc.PassageCurve;
import jline.api.mc.PassageMomentsResult;
import jline.api.mc.Smp_passage_moments;
import jline.api.pfqn.nc.Pfqn_cyclet_ofree;
import jline.api.pfqn.nc.PfqnCycletResult;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * First passage times in a Markov chain, and the exact cycle time along an
 * overtake-free path of a closed tree-like product-form network.
 *
 * <p>Reference: P. G. Harrison and W. J. Knottenbelt, "Passage Time
 * Distributions in Large Markov Chains", 2002.</p>
 *
 * <p>State and node indices are 0-based here and 1-based in the MATLAB
 * reference.</p>
 */
public class PassageTimeExample {

    /**
     * First passage times, semi-Markov passage times and an overtake-free cycle
     * time (passage_firstpassage.m).
     */
    public static void passage_firstpassage() {
        // 1. The time for an M/M/1/K queue to fill from empty.
        // This is a first passage into a STATE SET, which no response-time getter
        // can express: it is the chain reaching a marking, not a job finishing
        // service.
        final int K = 6;
        final double lambda = 1, mu = 1.5;
        final int n = K + 1;
        Matrix Q = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            if (i < n - 1) {
                Q.set(i, i + 1, lambda);
            }
            if (i > 0) {
                Q.set(i, i - 1, mu);
            }
        }
        Q = Ctmc_makeinfgen.ctmc_makeinfgen(Q);

        Matrix pi0 = new Matrix(1, n);
        pi0.set(0, 0, 1.0);                  // start empty
        int[] target = {n - 1};              // the full buffer

        PassageMomentsResult pm = Ctmc_passage_moments.ctmc_passage_moments(Q, pi0, target, 3);
        double[] m = pm.m;
        System.out.printf("Time to fill an M/M/1/%d from empty (lambda=%g, mu=%g)%n", K, lambda, mu);
        System.out.printf("  mean            = %.6f%n", m[0]);
        System.out.printf("  variance        = %.6f%n", m[1] - m[0] * m[0]);
        System.out.printf("  coeff. of var.  = %.6f%n",
                Math.sqrt(m[1] - m[0] * m[0]) / m[0]);
        System.out.printf("  from state K-1  = %.6f (one arrival away)%n", pm.mall.get(n - 2, 0));

        double[] tset = linspace(0, 4 * m[0], 400);
        PassageCurve curve = Ctmc_passage_time.ctmc_passage_time(Q, pi0, target, tset);
        System.out.printf("  P(fill <= mean) = %.6f%n", interp1(tset, curve.F, m[0]));

        // The mean hitting time from every state at once, the CTMC twin of
        // dtmc_hitting_time.
        Matrix h = Ctmc_hitting_time.ctmc_hitting_time(Q, target);
        double[] hr = new double[n];
        for (int i = 0; i < n; i++) {
            hr[i] = Math.round(h.get(i, 0) * 1000.0) / 1000.0;
        }
        System.out.printf("  hitting times   = %s%n", Arrays.toString(hr));

        // 2. The same passage with a non-exponential sojourn (semi-Markov).
        // The embedded chain is unchanged; only the holding-time law moves.
        // Nothing in a generator can express this, which is why the semi-Markov
        // route exists.
        double[] rate = new double[n];
        for (int i = 0; i < n; i++) {
            rate[i] = -Q.get(i, i);
        }
        Matrix P = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            if (rate[i] > 0) {
                for (int j = 0; j < n; j++) {
                    P.set(i, j, Q.get(i, j) / rate[i]);
                }
                P.set(i, i, 0.0);
            } else {
                P.set(i, i, 1.0);
            }
        }
        // Deterministic sojourns of the same mean: same embedded chain, tighter law.
        Matrix hmom = new Matrix(n, 3);
        for (int i = 0; i < n; i++) {
            double d = 1.0 / rate[i];
            hmom.set(i, 0, d);
            hmom.set(i, 1, d * d);
            hmom.set(i, 2, d * d * d);
        }
        double[] mD = Smp_passage_moments.smp_passage_moments(P, hmom, pi0, target, 3).m;
        System.out.printf("%nSame embedded chain, DETERMINISTIC sojourns of equal mean:%n");
        System.out.printf("  mean            = %.6f (unchanged, as it must be)%n", mD[0]);
        System.out.printf("  coeff. of var.  = %.6f (was %.6f)%n",
                Math.sqrt(mD[1] - mD[0] * mD[0]) / mD[0],
                Math.sqrt(m[1] - m[0] * m[0]) / m[0]);

        // 3. The cycle time of the tree-like network of Fig. 6 of the paper.
        double[] mu6 = {3, 5, 4, 6, 2, 1};
        double p12 = 0.2, p13 = 0.5, p14 = 0.3;
        double[] v = {1, p12, p13, p14, p12, p14};
        int N = 18;
        List<int[]> paths = new ArrayList<int[]>();
        paths.add(new int[]{0, 2});
        paths.add(new int[]{0, 1, 4});
        paths.add(new int[]{0, 3, 5});
        double[] pathprob = {p13, p12, p14};
        double[] tt = linspace(0, 40, 161);
        PfqnCycletResult cyc = Pfqn_cyclet_ofree.pfqn_cyclet_ofree(v, mu6, N, paths, tt,
                "auto", 3, pathprob, "euler", Pfqn_cyclet_ofree.DEFAULT_TOL);
        System.out.printf("%nTree network of Fig. 6, N = %d customers%n", N);
        System.out.printf("  moments  = %.5f  %.4f  %.3f%n",
                cyc.mom[0], cyc.mom[1], cyc.mom[2]);
        System.out.printf("  paper    = 6.12717  53.3067  612.887%n");
        System.out.printf("  routes   = %s%n", String.join(", ", cyc.method));
        System.out.printf("  P(cycle <= mean) = %.6f%n", interp1(tt, cyc.F, cyc.mom[0]));
    }

    /** The MATLAB linspace, as an array of npts points from a to b inclusive. */
    private static double[] linspace(double a, double b, int npts) {
        double[] t = new double[npts];
        for (int i = 0; i < npts; i++) {
            t[i] = a + (b - a) * i / (npts - 1.0);
        }
        return t;
    }

    /** Linear interpolation of y on the increasing grid x, the MATLAB interp1. */
    private static double interp1(double[] x, double[] y, double xq) {
        if (xq <= x[0]) {
            return y[0];
        }
        if (xq >= x[x.length - 1]) {
            return y[y.length - 1];
        }
        int k = Arrays.binarySearch(x, xq);
        if (k >= 0) {
            return y[k];
        }
        k = -k - 2;
        double w = (xq - x[k]) / (x[k + 1] - x[k]);
        return y[k] + w * (y[k + 1] - y[k]);
    }

    public static void main(String[] args) {
        passage_firstpassage();
    }
}
