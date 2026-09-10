/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.fes;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.jupiter.api.Test;

import jline.api.mam.Map2_fit_idc;
import jline.api.mam.Map_exponential;
import jline.api.mam.Map_idc;
import jline.api.mam.Map_moment;
import jline.api.mam.Map_scv;
import jline.api.pfqn.mva.Pfqn_mva;
import jline.io.Ret;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Contract of the MAP flow-equivalent server of Casale, Mi, Cherkasova and Smirni, IEEE
 * Trans. Soft. Eng. 37(5), 2011, Section 5.2.
 *
 * The mean inter-departure time of the aggregated subnetwork is the reciprocal of its
 * throughput, so for exponential service the recursion must reproduce exact MVA at every
 * population: the MAP flow-equivalent server degrades to the classic Norton one, which is
 * exact for product-form. Burstiness must survive the aggregation when the subnetwork
 * carries it. Expected values come from the MATLAB reference
 * (matlab/src/api/fes/fes_map_aggregate.m).
 */
public class FesMapAggregateTest {

    private static final double TOL = 1e-9;

    private static MatrixCell exp(double mean) {
        return Map_exponential.map_exponential(mean);
    }

    @Test
    public void testInterdepartureMeanIsInverseThroughput() {
        double mu1 = 1.5;
        double mu2 = 1.0;
        List<MatrixCell> station = Fes_map_levels.fes_map_levels(exp(1 / mu1), 10, 1);
        List<MatrixCell> fes = Fes_map_levels.fes_map_levels(exp(1 / mu2), 10, 1);

        Matrix L = new Matrix(2, 1);
        L.set(0, 0, 1 / mu1);
        L.set(1, 0, 1 / mu2);
        Matrix Z = new Matrix(1, 1);
        Matrix mi = new Matrix(1, 2);
        mi.set(0, 0, 1);
        mi.set(0, 1, 1);

        for (int n : new int[]{1, 2, 5, 10}) {
            MatrixCell T = Fes_map_interdeparture.fes_map_interdeparture(station, fes, n);
            FesMapMomentsResult mom = Fes_map_moments.fes_map_moments(T);
            Matrix N = new Matrix(1, 1);
            N.set(0, 0, n);
            double X = Pfqn_mva.pfqn_mva(L, N, Z, mi).X.get(0);
            assertEquals(1 / X, mom.e1, TOL * (1 / X), "mean inter-departure at n=" + n);
        }
    }

    @Test
    public void testEulerQuadratureMatchesLinearSolve() {
        List<MatrixCell> station = Fes_map_levels.fes_map_levels(exp(0.8), 6, 1);
        List<MatrixCell> fes = Fes_map_levels.fes_map_levels(exp(1.3), 6, 1);
        MatrixCell T = Fes_map_interdeparture.fes_map_interdeparture(station, fes, 6);

        // the quadrature propagates with exp(T0 dt) ~ I + T0 dt, so it is first order in
        // the step: the test asserts that rate rather than an arbitrary tolerance
        FesMapMomentsResult a = Fes_map_moments.fes_map_moments(T.get(0), T.get(1), "ssolve");
        double eCoarse = Math.abs(Fes_map_moments.fes_map_moments(T.get(0), T.get(1), "euler", 0.1).e1 - a.e1);
        double eFine = Math.abs(Fes_map_moments.fes_map_moments(T.get(0), T.get(1), "euler", 0.01).e1 - a.e1);
        assertTrue(eCoarse / a.e1 < 5e-2, "euler mean at the default step");
        assertTrue(eFine < 0.2 * eCoarse, "euler error must fall with the step, "
                + eCoarse + " -> " + eFine);
    }

    @Test
    public void testExponentialSubnetworkReproducesExactMva() {
        double[] rates = {2, 1.5, 1};
        List<MatrixCell> maps = new ArrayList<MatrixCell>();
        Matrix L = new Matrix(3, 1);
        Matrix mi = new Matrix(1, 3);
        for (int i = 0; i < 3; i++) {
            maps.add(exp(1 / rates[i]));
            L.set(i, 0, 1 / rates[i]);
            mi.set(0, i, 1);
        }
        Matrix Z = new Matrix(1, 1);

        int nmax = 8;
        FesMapAggregateResult res = Fes_map_aggregate.fes_map_aggregate(maps, new double[]{1, 1, 1}, nmax);
        for (int k = 1; k <= nmax; k++) {
            Matrix N = new Matrix(1, 1);
            N.set(0, 0, k);
            double X = Pfqn_mva.pfqn_mva(L, N, Z, mi).X.get(0);
            assertEquals(X, res.throughput[k - 1], TOL * X, "aggregate throughput at n=" + k);
        }
    }

    @Test
    public void testDelayAndMultiserverSubnetwork() {
        // delay of rate 2 and a two-server queue of rate 1, exact birth-death reference
        List<MatrixCell> maps = Arrays.asList(exp(0.5), exp(1.0));
        int nmax = 6;
        FesMapAggregateResult res = Fes_map_aggregate.fes_map_aggregate(
                maps, new double[]{Double.POSITIVE_INFINITY, 2}, nmax);

        for (int k = 1; k <= nmax; k++) {
            double[] p = birthDeath(k);
            double X = 0;
            for (int j = 0; j <= k; j++) {
                X += p[j] * Math.min(j, 2);
            }
            assertEquals(X, res.throughput[k - 1], 1e-8 * X, "aggregate throughput at n=" + k);
        }
    }

    @Test
    public void testBurstinessSurvivesAggregation() {
        // a hyper-exponential station with a large index of dispersion, folded with an
        // exponential one, must leave the flow-equivalent server bursty
        MatrixCell bursty = new MatrixCell(
                new Matrix(new double[][]{{-1.9, 0}, {0, -0.1}}),
                new Matrix(new double[][]{{1.71, 0.19}, {0.01, 0.09}}));
        List<MatrixCell> maps = Arrays.asList(bursty, exp(1.2));
        FesMapAggregateResult res = Fes_map_aggregate.fes_map_aggregate(maps, new double[]{1, 1}, 5);

        for (int k = 1; k <= 5; k++) {
            assertEquals(Map2_fit_idc.STATUS_EXACT, res.status[k - 1], "fit status at n=" + k);
            MatrixCell f = res.fes.get(k - 1);
            assertEquals(res.moments[3][k - 1], Map_idc.map_idc(f), 1e-6 * res.moments[3][k - 1],
                    "index of dispersion carried at n=" + k);
            assertEquals(res.moments[0][k - 1], Map_moment.map_moment(f, 1), 1e-9,
                    "mean carried at n=" + k);
            assertTrue(Map_scv.map_scv(f) > 1, "the aggregate must stay overdispersed at n=" + k);
        }
    }

    @Test
    public void testFitFallsBackToExponentialWhenNotBursty() {
        Ret.mamMAPFitIdcReturn fit = Map2_fit_idc.map2_fit_idc(1.0, 1.5, 4.0, 0.7);
        assertEquals(Map2_fit_idc.STATUS_EXPONENTIAL, fit.status);
        assertEquals(1.0, Map_moment.map_moment(fit.MAP, 1), TOL);
    }

    @Test
    public void testGridSkipsLevelsOnlyAboveTwenty() {
        assertEquals(20, Fes_map_grid.fes_map_grid(20).length);
        assertTrue(Fes_map_grid.fes_map_grid(100).length < 100);
        int[] grid = Fes_map_grid.fes_map_grid(100);
        assertEquals(1, grid[0]);
        assertEquals(100, grid[grid.length - 1]);
    }

    @Test
    public void testReducedModelSolveMatchesExactMva() {
        // the closed model of a delay and an aggregated exponential subnetwork is
        // product-form, so the reduced solve must return exact MVA
        double[] rates = {2, 1.5, 1};
        double Z = 0.8;
        List<MatrixCell> maps = new ArrayList<MatrixCell>();
        Matrix L = new Matrix(3, 1);
        Matrix mi = new Matrix(1, 3);
        for (int i = 0; i < 3; i++) {
            maps.add(exp(1 / rates[i]));
            L.set(i, 0, 1 / rates[i]);
            mi.set(0, i, 1);
        }
        Matrix Zm = new Matrix(1, 1);
        Zm.set(0, 0, Z);

        for (int n : new int[]{1, 3, 6, 10}) {
            FesMapAggregateResult agg = Fes_map_aggregate.fes_map_aggregate(maps, new double[]{1, 1, 1}, n);
            FesMapSolveResult sol = Fes_map_solve.fes_map_solve(agg.fes, exp(Z), n);
            Matrix N = new Matrix(1, 1);
            N.set(0, 0, n);
            Ret.pfqnMVA ref = Pfqn_mva.pfqn_mva(L, N, Zm, mi);
            assertEquals(ref.X.get(0), sol.X, TOL * ref.X.get(0), "throughput at N=" + n);
            double R = 0;
            for (int i = 0; i < 3; i++) {
                R += ref.R.get(i, 0);
            }
            assertEquals(R, sol.R, 1e-8 * R, "response time at N=" + n);
            double mass = 0;
            for (int k = 0; k <= n; k++) {
                mass += sol.pk[k];
            }
            assertEquals(1.0, mass, 1e-12, "level law at N=" + n);
        }
    }

    @Test
    public void testDeaggregationMatchesExactMvaPerStation() {
        double[] rates = {2, 1.5, 1};
        double Z = 0.8;
        int n = 8;
        List<MatrixCell> maps = new ArrayList<MatrixCell>();
        Matrix L = new Matrix(3, 1);
        Matrix mi = new Matrix(1, 3);
        double[] demands = new double[3];
        for (int i = 0; i < 3; i++) {
            maps.add(exp(1 / rates[i]));
            demands[i] = 1 / rates[i];
            L.set(i, 0, demands[i]);
            mi.set(0, i, 1);
        }
        Matrix Zm = new Matrix(1, 1);
        Zm.set(0, 0, Z);

        FesMapAggregateResult agg = Fes_map_aggregate.fes_map_aggregate(maps, new double[]{1, 1, 1}, n);
        FesMapSolveResult sol = Fes_map_solve.fes_map_solve(agg.fes, exp(Z), n);
        FesMapDeaggregateResult de = Fes_map_deaggregate.fes_map_deaggregate(
                sol.pk, demands, new double[]{1, 1, 1}, new boolean[]{false, false, false});

        Matrix N = new Matrix(1, 1);
        N.set(0, 0, n);
        Ret.pfqnMVA ref = Pfqn_mva.pfqn_mva(L, N, Zm, mi);
        for (int i = 0; i < 3; i++) {
            assertEquals(ref.Q.get(i, 0), de.Q[i], 1e-8 * ref.Q.get(i, 0), "queue length at station " + i);
            assertEquals(ref.U.get(i, 0), de.U[i], 1e-8 * ref.U.get(i, 0), "utilization at station " + i);
            assertEquals(ref.X.get(0), de.X[i], 1e-8 * ref.X.get(0), "throughput at station " + i);
        }
    }

    /** Stationary distribution of the delay plus two-server queue reference model. */
    private static double[] birthDeath(int k) {
        double[] p = new double[k + 1];
        p[0] = 1;
        double sum = 1;
        for (int j = 1; j <= k; j++) {
            p[j] = p[j - 1] * ((k - j + 1) * 2.0) / (Math.min(j, 2) * 1.0);
            sum += p[j];
        }
        for (int j = 0; j <= k; j++) {
            p[j] /= sum;
        }
        return p;
    }
}
