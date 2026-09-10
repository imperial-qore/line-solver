/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.dpfqn;

import jline.api.dqsys.Bernoulli1Result;
import jline.api.dqsys.Dqsys_bernoulli1;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * Discrete-time normalizing constants and Bernoulli server, checked against
 * brute-force enumeration of the underlying discrete-time Markov chains.
 *
 * <p>The reference values were produced by solving the chain directly: the
 * closed cycle by its transition matrix over the compositions of N, and the
 * single node by its birth-death chain. They are the same numbers the MATLAB
 * and native-Python ports return.</p>
 */
public class DpfqnNcTest {

    private static final double TOL = 1e-9;

    /** Cycle of three Bernoulli servers, N = 5. */
    @Test
    public void cycleReproducesTheChainStationaryLaw() {
        double[] p = {0.5, 0.25, 0.7};
        int N = 5;
        DpfqnNcResult nc = Dpfqn_nc.dpfqn_nc(p, N);

        // Throughput is the same at every node and equals G_1(N,J)/G(N,J).
        assertEquals(0.24828195, nc.throughput(), 1e-8);
        assertEquals(6.9088216466, nc.lG, 1e-9);

        double[] expectedQ = {0.721909, 3.866843, 0.411247};
        double[] expectedU = {0.496564, 0.993128, 0.354689};
        double total = 0;
        for (int j = 0; j < p.length; j++) {
            double q = 1 - p[j];
            double tail = 0;
            for (int n = 1; n <= N; n++) {
                tail += Math.pow(q / p[j], n) / q * nc.G1[N - n + 1] / nc.G;
            }
            assertEquals(expectedQ[j], tail, 1e-6);
            assertEquals(expectedU[j], nc.throughput() / p[j], 1e-6);
            total += tail;
        }
        assertEquals(N, total, 1e-9, "the marginals must account for every job");
    }

    /** The load-dependent route must reproduce the state independent one. */
    @Test
    public void loadDependentAgreesWithStateIndependent() {
        double[] p = {0.5, 0.25, 0.7};
        int N = 5;
        DpfqnNcResult nc = Dpfqn_nc.dpfqn_nc(p, N);
        double[][] P = new double[p.length][N];
        for (int j = 0; j < p.length; j++) {
            for (int n = 0; n < N; n++) {
                P[j][n] = p[j];
            }
        }
        DpfqnNcLdResult ncld = Dpfqn_ncld.dpfqn_ncld(P, N);
        assertEquals(nc.lG, ncld.lG, TOL);
        for (int j = 0; j < p.length; j++) {
            double[] marg = ncld.marginal(j);
            double sum = 0;
            double q = 0;
            double t = 0;
            for (int n = 0; n <= N; n++) {
                sum += marg[n];
                q += marg[n] * n;
                if (n >= 1) {
                    t += marg[n] * p[j];
                }
            }
            assertEquals(1.0, sum, TOL, "the marginal must be a probability law");
            assertEquals(nc.throughput(), t, 1e-9, "throughput is equal at every node");
            assertEquals(1.0 - marg[0], nc.throughput() / p[j], 1e-9);
            assertEquals(q, q, TOL);
        }
    }

    /** Station 2 at p(n) = p min(n,2), against the brute-force chain. */
    @Test
    public void loadDependentCycleMatchesTheChain() {
        double[] p = {0.5, 0.25, 0.7};
        int N = 5;
        double[][] P = new double[3][N];
        for (int n = 1; n <= N; n++) {
            P[0][n - 1] = p[0];
            P[1][n - 1] = p[1] * Math.min(n, 2);
            P[2][n - 1] = p[2];
        }
        DpfqnNcLdResult nc = Dpfqn_ncld.dpfqn_ncld(P, N);
        double[] expectedQ = {1.793601, 2.391174, 0.815225};
        double[] expectedU = {0.845151, 0.947567, 0.603679};
        for (int j = 0; j < 3; j++) {
            double[] marg = nc.marginal(j);
            double q = 0;
            double t = 0;
            for (int n = 0; n <= N; n++) {
                q += marg[n] * n;
                if (n >= 1) {
                    t += marg[n] * P[j][n - 1];
                }
            }
            assertEquals(expectedQ[j], q, 1e-6);
            assertEquals(expectedU[j], 1 - marg[0], 1e-6);
            assertEquals(0.422576, t, 1e-6);
        }
    }

    /** Geo/Geo/1, the closed form of corollary 2.7. */
    @Test
    public void bernoulliServerUnboundedBuffer() {
        Bernoulli1Result r = Dqsys_bernoulli1.dqsys_bernoulli1(0.2, 0.5);
        assertEquals(0.533333333333, r.meanQueueLength, 1e-9);
        assertEquals(2.666666666667, r.meanSojournTime, 1e-9);
        assertEquals(0.4, r.utilization, 1e-12);
        assertEquals(0.2, r.throughput, 1e-12);
        assertEquals(0.0, r.lossProb, 0.0);
    }

    /** Geo/Geo/1/4, the loss system of corollary 2.8. */
    @Test
    public void bernoulliServerLossSystem() {
        Bernoulli1Result r = Dqsys_bernoulli1.dqsys_bernoulli1(
                new double[]{0.2}, new double[]{0.5}, 4);
        double[] expected = {0.601504, 0.300752, 0.075188, 0.018797, 0.003759};
        for (int n = 0; n <= 4; n++) {
            assertEquals(expected[n], r.pmf[n], 1e-6);
        }
        assertEquals(0.0037594, r.lossProb, 1e-7);
        assertEquals(0.199248, r.throughput, 1e-6);
    }

    /**
     * The arrival law of theorem 2.11 is not the time-stationary law. With a
     * state independent Bernoulli stream it is geometric with ratio
     * r = b(1-p)/(p(1-b)), the EAS-convention queue length law.
     */
    @Test
    public void arrivalLawIsNotTheTimeStationaryLaw() {
        double b = 0.2;
        double p = 0.5;
        int L = 40;
        Bernoulli1Result r = Dqsys_bernoulli1.dqsys_bernoulli1(
                new double[]{b}, new double[]{p}, L);
        double ratio = b * (1 - p) / (p * (1 - b));
        for (int n = 0; n < 6; n++) {
            assertEquals((1 - ratio) * Math.pow(ratio, n), r.arrivalPmf[n], 1e-9);
        }
        // and it genuinely differs from pi
        assertEquals(0.75, r.arrivalPmf[0], 1e-9);
        assertEquals(0.6, r.pmf[0], 1e-6);
    }

    /** The load-dependent Bernoulli server of example 2.10. */
    @Test
    public void multiserverApproximationMatchesTheChain() {
        int L = 20;
        int c = 3;
        double[] p = new double[L];
        for (int n = 1; n <= L; n++) {
            p[n - 1] = 0.3 * Math.min(n, c);
        }
        Bernoulli1Result r = Dqsys_bernoulli1.dqsys_bernoulli1(new double[]{0.6}, p, L);
        assertEquals(2.064368, r.meanQueueLength, 1e-6);
        assertEquals(0.954023, r.utilization, 1e-6);
        assertEquals(0.6, r.throughput, 1e-6);
        assertEquals(0.045977, r.pmf[0], 1e-6);
        assertEquals(0.114943, r.arrivalPmf[0], 1e-6);
    }

    /** A service probability outside (0,1) has no product form. */
    @Test
    public void rejectsDegenerateServiceProbabilities() {
        assertThrows(IllegalArgumentException.class,
                () -> Dpfqn_nc.dpfqn_nc(new double[]{1.0, 0.5}, 3));
        assertThrows(IllegalArgumentException.class,
                () -> Dpfqn_nc.dpfqn_nc(new double[]{0.0, 0.5}, 3));
    }
}
