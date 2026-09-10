/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.dpfqn;

/**
 * Convolution over the discrete-time product form of a closed cycle of state
 * dependent Bernoulli servers.
 *
 * <p>Port of matlab/src/api/dpfqn/dpfqn_ncld.m. With q_j(n) = 1 - p_j(n) the
 * queue length vector has the stationary law of Daduna (2001), theorem 3.2,</p>
 *
 * <pre>
 *   pi(n_1,...,n_J) = prod_j w_j(n_j) / G(N,J)
 *   w_j(n) = prod_{h=1}^{n-1} q_j(h) / prod_{h=1}^{n} p_j(h),   w_j(0) = 1
 * </pre>
 *
 * <p>which reduces to the state independent form of {@link Dpfqn_nc} when
 * p_j(n) does not depend on n. Note that w_j misses the factor q_j(n_j) in the
 * numerator, so the extra 1/q_j on a busy node in the state independent case is
 * not tied to the node being non-empty but to its actual queue length.</p>
 *
 * <p>The constants are assembled by truncated convolution of the per-node
 * weight vectors, together with the J complement constants (the cycle with node
 * j removed) obtained from a prefix/suffix pass. Deconvolution is never used,
 * so a node with a near-unit service probability does not destroy the accuracy
 * of the other marginals.</p>
 *
 * @see Dpfqn_nc
 */
public final class Dpfqn_ncld {

    private Dpfqn_ncld() {}

    /**
     * @param P service probabilities, {@code P[j][n-1] = p_j(n)} in (0,1]
     * @param N number of customers cycling in the J nodes
     * @return the normalizing, complement and arrival constants
     */
    public static DpfqnNcLdResult dpfqn_ncld(double[][] P, int N) {
        if (N < 0) {
            throw new IllegalArgumentException("N must be a non-negative integer");
        }
        if (P == null || P.length == 0) {
            throw new IllegalArgumentException(
                    "P must be a non-empty matrix of service probabilities");
        }
        int J = P.length;
        for (int j = 0; j < J; j++) {
            if (P[j] == null || P[j].length < N) {
                throw new IllegalArgumentException(String.format(
                        "P must supply p_j(n) for every n = 1..%d", N));
            }
        }
        if (N == 0) {
            double[][] ones = new double[J][1];
            for (int j = 0; j < J; j++) {
                ones[j][0] = 1.0;
            }
            return new DpfqnNcLdResult(0.0, new double[]{1.0}, ones, ones, ones);
        }
        for (int j = 0; j < J; j++) {
            for (int n = 0; n < N; n++) {
                double v = P[j][n];
                if (Double.isNaN(v) || Double.isInfinite(v) || v <= 0 || v > 1) {
                    throw new IllegalArgumentException(
                            "service probabilities must be real and in the interval (0,1]");
                }
                // p_j(n)=1 is admissible only at the last reachable population,
                // where the missing q_j(n) never multiplies any weight.
                if (n < N - 1 && v >= 1) {
                    throw new IllegalArgumentException(
                            "service probabilities below the population bound must be strictly less than 1");
                }
            }
        }

        // Per-node log weights of theorem 3.2, then a per-node shift so that
        // every weight vector peaks at one. The shifts cancel in every ratio.
        double[][] lw = new double[J][N + 1];
        for (int j = 0; j < J; j++) {
            double cq = 0.0;
            double cp = 0.0;
            for (int n = 1; n <= N; n++) {
                cp += Math.log(P[j][n - 1]);
                lw[j][n] = cq - cp;
                cq += (P[j][n - 1] < 1.0) ? Math.log(1.0 - P[j][n - 1]) : Double.NEGATIVE_INFINITY;
            }
        }
        double[] shift = new double[J];
        double[][] W = new double[J][N + 1];
        for (int j = 0; j < J; j++) {
            double mx = Double.NEGATIVE_INFINITY;
            for (int n = 0; n <= N; n++) {
                if (lw[j][n] > mx) {
                    mx = lw[j][n];
                }
            }
            shift[j] = mx;
            for (int n = 0; n <= N; n++) {
                W[j][n] = Math.exp(lw[j][n] - mx);
            }
        }

        double[][] Wa = new double[J][N + 1];
        for (int j = 0; j < J; j++) {
            Wa[j][0] = W[j][0];
            for (int n = 1; n <= N; n++) {
                Wa[j][n] = W[j][n] * (1.0 - P[j][n - 1]);
            }
        }

        // Prefix/suffix convolution: pre[j] covers nodes 0..j-1 and suf[j]
        // covers nodes j..J-1, so the complement of node j is their product.
        double[][] pre = new double[J + 1][N + 1];
        pre[0][0] = 1.0;
        for (int j = 0; j < J; j++) {
            pre[j + 1] = conv(pre[j], W[j], N);
        }
        double[][] suf = new double[J + 1][N + 1];
        suf[J][0] = 1.0;
        for (int j = J - 1; j >= 0; j--) {
            suf[j] = conv(suf[j + 1], W[j], N);
        }
        double[] G = pre[J];
        double[][] Gc = new double[J][];
        for (int j = 0; j < J; j++) {
            Gc[j] = conv(pre[j], suf[j + 1], N);
        }

        double totalShift = 0.0;
        for (int j = 0; j < J; j++) {
            totalShift += shift[j];
        }
        return new DpfqnNcLdResult(Math.log(G[N]) + totalShift, G, W, Gc, Wa);
    }

    /** Convolution of two population tables truncated at N. */
    private static double[] conv(double[] a, double[] b, int N) {
        double[] c = new double[N + 1];
        for (int k = 0; k <= N; k++) {
            double s = 0.0;
            for (int m = 0; m <= k; m++) {
                s += a[m] * b[k - m];
            }
            c[k] = s;
        }
        return c;
    }
}
