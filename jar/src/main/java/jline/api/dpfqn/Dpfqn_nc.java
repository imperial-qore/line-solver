/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.dpfqn;

/**
 * Buzen-style recursions for the discrete-time closed cycle of Bernoulli
 * servers with state independent service probabilities.
 *
 * <p>Port of matlab/src/api/dpfqn/dpfqn_nc.m. With q_j = 1 - p_j the queue
 * length vector has the product form of Daduna (2001), corollary 3.4,</p>
 *
 * <pre>
 *   pi(n_1,...,n_J) = prod_j (q_j/p_j)^n_j (1/q_j)^{1{n_j&gt;0}} / G(N,J)
 * </pre>
 *
 * <p>whose extra factor on the busy nodes is what separates it from the
 * continuous-time Gordon-Newell form: a homogeneous cycle is uniform on the
 * state space in continuous time and is not here. G obeys proposition 3.18,</p>
 *
 * <pre>
 *   G(k,j) = G(k,j-1) + (q_j/p_j) G(k-1,j) + G(k-1,j-1)
 * </pre>
 *
 * <p>with G(0,j) = 1 and G(k,0) = 0 for k &gt;= 1. Unlike the continuous-time
 * convolution algorithm this recursion is not invariant to the numbering of the
 * nodes. The arrival constants obey proposition 3.19,</p>
 *
 * <pre>
 *   G1(k,J) = G(k-1,J-1) + (q_J/p_J) G1(k-1,J),   k &gt;= 3
 * </pre>
 *
 * <p>with G1(1,J) = 1 and G1(2,J) = q_1/p_1 + sum_{j&gt;=2} 1/p_j. By lemma 7.3
 * the arrival constant is the same at every node, so one family suffices, and
 * three consequences are what the discrete-time analyzer consumes:</p>
 *
 * <pre>
 *   throughput per slot   X   = G1(N,J) / G(N,J)   (equal at every node)
 *   utilization           U_j = X / p_j
 *   tail probability      P(X_j &gt;= k)
 *                         = (q_j/p_j)^k (1/q_j) G1(N-k+1,J) / G(N,J)
 * </pre>
 *
 * <p>The last identity is corollary 3.20(a) with its index corrected: as
 * printed there the right-hand side evaluates to P(X_j &gt;= k+1). Corollary
 * 3.20(c), which transfers the tail from node 1 to node j, holds for k &gt;= 1
 * only; at k = 0 both tails are 1 while the stated ratio is q_1/q_j.</p>
 *
 * @see Dpfqn_ncld
 */
public final class Dpfqn_nc {

    /** Table entries above this magnitude trigger a global rescale. */
    private static final double RESCALE = 1e250;

    private Dpfqn_nc() {}

    /**
     * @param p per-slot service completion probabilities p_j in (0,1)
     * @param N number of customers cycling in the J nodes
     * @return the time-stationary and arrival normalizing constants
     */
    public static DpfqnNcResult dpfqn_nc(double[] p, int N) {
        if (p == null || p.length == 0) {
            throw new IllegalArgumentException("the cycle must contain at least one node");
        }
        for (int j = 0; j < p.length; j++) {
            if (Double.isNaN(p[j]) || Double.isInfinite(p[j]) || p[j] <= 0 || p[j] >= 1) {
                throw new IllegalArgumentException(
                        "service probabilities must be real and in the open interval (0,1)");
            }
        }
        if (N < 0) {
            throw new IllegalArgumentException("N must be a non-negative integer");
        }
        int J = p.length;
        double[] x = new double[J];
        for (int j = 0; j < J; j++) {
            x[j] = (1.0 - p[j]) / p[j];
        }

        // The recursion is linear and homogeneous in the whole table, so
        // rescaling every entry at once preserves it; that is what keeps
        // (q/p)^N from overflowing for a slow node and a large population. G1
        // is rescaled in step so the two families stay in one common scale.
        double[][] Gt = new double[N + 1][J + 1];
        for (int j = 0; j <= J; j++) {
            Gt[0][j] = 1.0;
        }
        double[] G1 = new double[N + 1];
        double lscale = 0.0;
        for (int k = 1; k <= N; k++) {
            for (int j = 1; j <= J; j++) {
                Gt[k][j] = Gt[k][j - 1] + x[j - 1] * Gt[k - 1][j] + Gt[k - 1][j - 1];
            }
            if (k == 1) {
                G1[1] = Gt[0][0];
            } else if (k == 2) {
                double s = x[0];
                for (int j = 1; j < J; j++) {
                    s += 1.0 / p[j];
                }
                G1[2] = s * Gt[0][0];
            } else {
                G1[k] = Gt[k - 1][J - 1] + x[J - 1] * G1[k - 1];
            }
            double mx = 0.0;
            for (int j = 0; j <= J; j++) {
                if (Gt[k][j] > mx) {
                    mx = Gt[k][j];
                }
            }
            if (mx > RESCALE) {
                for (int a = 0; a <= N; a++) {
                    for (int b = 0; b <= J; b++) {
                        Gt[a][b] /= mx;
                    }
                }
                for (int a = 0; a <= N; a++) {
                    G1[a] /= mx;
                }
                lscale += Math.log(mx);
            }
        }
        double G = Gt[N][J];
        return new DpfqnNcResult(Math.log(G) + lscale, G, G1);
    }
}
