/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.pfqn.nc;

import jline.api.pfqn.ld.Pfqn_ncld;
import jline.io.Ret;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Chen-O'Cinneide regularization (Pfqn_mcmc).
 *
 * <p>W. Chen, C. A. O'Cinneide, "Towards a Polynomial-Time Randomized Algorithm for
 * Closed Product-Form Networks", ACM TOMACS 8(3):227-253, 1998.</p>
 *
 * <p>Two kinds of assertion are made here, and only the first kind is statistical.</p>
 *
 * <p>The throughput checks compare the estimator against the EXACT ratio G(N-e_r)/G(N)
 * from convolution -- Pfqn_ca for the single-server models, Pfqn_ncld with
 * mu_i(k)=min(k,c_i) for the multiserver one -- at a tolerance the measured error clears
 * with room to spare. They are seeded, so they are deterministic runs of a random
 * algorithm rather than flaky tests.</p>
 *
 * <p>The queue-length check is NOT statistical and holds to machine precision: every
 * state of the regularized chain satisfies sum_i Y(i,r) = N(r), and the reported Q is a
 * weighted average of those states, so with no delay station the columns of Q sum to N
 * exactly whatever the sample path was. That invariant is what catches an indexing or
 * weighting error, which a tolerance on a noisy mean cannot.</p>
 */
public class Pfqn_mcmcTest {

    /** Monte Carlo tolerance on the throughput ratios, well above the measured error. */
    private static final double RTOL = 0.02;

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) m.set(0, i, v[i]);
        return m;
    }

    private static Matrix mat(double[][] v) {
        Matrix m = new Matrix(v.length, v[0].length);
        for (int i = 0; i < v.length; i++) {
            for (int j = 0; j < v[0].length; j++) m.set(i, j, v[i][j]);
        }
        return m;
    }

    private static SolverOptions opts(int samples) {
        SolverOptions o = new SolverOptions();
        o.samples = samples;
        o.seed = 23000;
        return o;
    }

    /** Exact load-dependent convolution, the multiserver reference. */
    private static SolverOptions ncldOpts() {
        SolverOptions o = new SolverOptions();
        o.method = "exact";
        o.tol = 1e-12;
        return o;
    }

    /** X(r) = G(N-e_r)/G(N) by convolution, the quantity Pfqn_mcmc estimates. */
    private static double[] exactRatios(Matrix L, Matrix N, Matrix Z) {
        int R = N.getNumCols();
        double lG = Pfqn_ca.pfqn_ca(L, N, Z).lG;
        double[] X = new double[R];
        for (int r = 0; r < R; r++) {
            Matrix Nr = N.copy();
            Nr.set(0, r, Nr.get(0, r) - 1);
            X[r] = Math.exp(Pfqn_ca.pfqn_ca(L, Nr, Z).lG - lG);
        }
        return X;
    }

    /** The 3 single-server stations plus IS station of Example 5.1 of the paper. */
    private static Matrix example51L() {
        double[] mu = {0.2, 0.5, 0.8};
        int[][] sets = {{1, 1, 1}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1}};
        Matrix L = new Matrix(3, 4);
        for (int c = 0; c < 4; c++) {
            for (int i = 0; i < 3; i++) {
                if (sets[c][i] == 1) L.set(i, c, 1.0 / mu[i]);
            }
        }
        return L;
    }

    private static Matrix example51Z() {
        return row(0.0, 2.0, 2.0, 2.0);
    }

    @Test
    public void example51MatchesTheExactRatios() {
        Matrix L = example51L();
        Matrix N = row(3, 3, 3, 3);
        Matrix Z = example51Z();
        Ret.pfqnMcmc res = Pfqn_mcmc.pfqn_mcmc(L, N, Z, opts(200000));
        double[] exact = exactRatios(L, N, Z);
        for (int r = 0; r < 4; r++) {
            assertEquals(exact[r], res.X.get(0, r), RTOL * exact[r]);
            // Every interval is a real interval around the estimate, and every standard
            // error is a positive number rather than a NaN from a zero denominator.
            assertTrue(res.Xse.get(0, r) > 0);
            assertTrue(res.Xlo.get(0, r) < res.X.get(0, r));
            assertTrue(res.X.get(0, r) < res.Xhi.get(0, r));
        }
        assertEquals(30, res.batches);
        assertEquals(res.samples / 10, res.burnin);
    }

    /**
     * Tenfold more work must buy roughly sqrt(10) less error, not a fixed floor.
     *
     * <p>A biased estimator -- a mis-weighted holding time, a warm-up that never ends --
     * passes the tolerance check above at one sample size and then stops improving.
     * Comparing two sizes separates noise from bias without pinning either run's value.</p>
     */
    @Test
    public void theEstimatorIsConsistent() {
        Matrix L = example51L();
        Matrix N = row(3, 3, 3, 3);
        Matrix Z = example51Z();
        double[] exact = exactRatios(L, N, Z);
        double[] err = new double[2];
        int[] sizes = {50000, 500000};
        for (int k = 0; k < 2; k++) {
            Ret.pfqnMcmc res = Pfqn_mcmc.pfqn_mcmc(L, N, Z, opts(sizes[k]));
            double worst = 0.0;
            for (int r = 0; r < 4; r++) {
                worst = Math.max(worst, Math.abs(res.X.get(0, r) - exact[r]) / exact[r]);
            }
            err[k] = worst;
        }
        assertTrue(err[1] < 0.5 * err[0], "error " + err[1] + " did not shrink below half of " + err[0]);
    }

    /** sum_i Y(i,r) = N(r) in every state, so the weighted average inherits it. */
    @Test
    public void queueLengthsConserveThePopulationExactly() {
        Matrix L = mat(new double[][]{{0.6, 0.2}, {0.3, 0.5}, {0.1, 0.4}});
        Matrix N = row(4, 3);
        Ret.pfqnMcmc res = Pfqn_mcmc.pfqn_mcmc(L, N, null, opts(50000));
        for (int r = 0; r < 2; r++) {
            double sum = 0.0;
            for (int i = 0; i < 3; i++) sum += res.Q.get(i, r);
            assertEquals(N.get(0, r), sum, 1e-12);
        }
    }

    /**
     * The paper's own selling point (its Tables IV and V): exact multiservers.
     *
     * <p>The reference is the load-dependent convolution with mu_i(k) = min(k, c_i), the
     * same product form the regularized chain samples, so any disagreement beyond the
     * Monte Carlo error is an error in the Psi_i(Y_i) = min(s_i, Y_i) rate.</p>
     */
    @Test
    public void multiserverMatchesTheExactLoadDependentConstant() {
        Matrix L = mat(new double[][]{{0.6, 0.2}, {0.3, 0.5}, {0.1, 0.4}});
        Matrix N = row(4, 3);
        Matrix s = new Matrix(3, 1);
        s.set(0, 0, 2.0);
        s.set(1, 0, 1.0);
        s.set(2, 0, 3.0);
        int Ntot = 7;
        Matrix mu = new Matrix(3, Ntot);
        for (int i = 0; i < 3; i++) {
            for (int k = 1; k <= Ntot; k++) mu.set(i, k - 1, Math.min(k, s.get(i, 0)));
        }
        Matrix Z = row(0.0, 0.0);
        double lG = Pfqn_ncld.pfqn_ncld(L, N, Z, mu, ncldOpts()).lG;
        double[] exact = new double[2];
        for (int r = 0; r < 2; r++) {
            Matrix Nr = N.copy();
            Nr.set(0, r, Nr.get(0, r) - 1);
            exact[r] = Math.exp(Pfqn_ncld.pfqn_ncld(L, Nr, Z, mu, ncldOpts()).lG - lG);
        }
        Ret.pfqnMcmc res = Pfqn_mcmc.pfqn_mcmc(L, N, null, s, opts(400000));
        for (int r = 0; r < 2; r++) {
            assertEquals(exact[r], res.X.get(0, r), RTOL * exact[r]);
            double sum = 0.0;
            for (int i = 0; i < 3; i++) sum += res.Q.get(i, r);
            assertEquals(N.get(0, r), sum, 1e-12);
        }
    }

    /**
     * The chain lives on the integer lattice, so a fractional N has no state space.
     *
     * <p>This is what makes mcmc unusable as a SolverLN layer method, which is correct
     * rather than unfortunate: rounding would answer a question nobody asked.</p>
     */
    @Test
    public void aFractionalPopulationIsRefusedNotRounded() {
        Matrix L = mat(new double[][]{{0.6, 0.2}, {0.3, 0.5}});
        RuntimeException e = assertThrows(RuntimeException.class,
                () -> Pfqn_mcmc.pfqn_mcmc(L, row(2.5, 1.0), null, opts(1000)));
        assertTrue(e.getMessage().contains("fractional"));
    }

    @Test
    public void anInfinitePopulationIsRefused() {
        Matrix L = mat(new double[][]{{0.6, 0.2}, {0.3, 0.5}});
        RuntimeException e = assertThrows(RuntimeException.class,
                () -> Pfqn_mcmc.pfqn_mcmc(L, row(Double.POSITIVE_INFINITY, 1.0), null, opts(1000)));
        assertTrue(e.getMessage().contains("closed model"));
    }

    @Test
    public void aPopulatedClassWithNoDemandAnywhereIsRefused() {
        Matrix L = mat(new double[][]{{0.6, 0.0}, {0.3, 0.0}});
        RuntimeException e = assertThrows(RuntimeException.class,
                () -> Pfqn_mcmc.pfqn_mcmc(L, row(2, 1), null, opts(1000)));
        assertTrue(e.getMessage().contains("no demand"));
    }

    @Test
    public void anEmptyNetworkReturnsZerosWithoutSimulating() {
        Matrix L = mat(new double[][]{{0.6, 0.2}, {0.3, 0.5}});
        Ret.pfqnMcmc res = Pfqn_mcmc.pfqn_mcmc(L, row(0, 0), null, opts(1000));
        assertEquals(0.0, res.X.elementSum(), 0.0);
        assertEquals(0.0, res.Q.elementSum(), 0.0);
        assertEquals(0, res.batches);
        assertEquals(0L, res.samples);
        assertEquals(0L, res.burnin);
    }

    /**
     * The method estimates ratios and never forms G, so Pfqn_nc supplies lG from BLE.
     *
     * <p>Pinning it here records that the number is deliberate: it cancels out of every
     * mean value the analyzer reports, and only getProbNormConstAggr reads it.</p>
     */
    @Test
    public void pfqnNcMcmcArmReturnsTheBleConstantNotAPaperResult() {
        Matrix L = example51L();
        Matrix N = row(3, 3, 3, 3);
        Matrix Z = example51Z();
        SolverOptions o = opts(20000);
        o.method = "mcmc";
        Ret.pfqnNcXQ ret = Pfqn_nc.pfqn_nc(new Matrix(1, 4), L, N, Z, o);
        assertEquals(Pfqn_ble.pfqn_ble(L, N, Z).lG, ret.lG, 1e-9);
        // and the X/Q channel carried the simulation's mean values out
        assertTrue(ret.X.elementSum() > 0);
        assertTrue(ret.Q.elementSum() > 0);
    }
}
