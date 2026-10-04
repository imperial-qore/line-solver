/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.dpfqn;

/**
 * Normalizing constants of a discrete-time closed cycle of state dependent
 * Bernoulli servers.
 *
 * <p>All four tables are in one common scale, so that the marginal law
 * {@code P(X_j = n) = W[j][n] * Gc[j][N-n] / G[N]} is exact whatever the scale;
 * {@link #lG} carries the true log constant.</p>
 */
public final class DpfqnNcLdResult {

    /** Logarithm of the normalizing constant G(N,J). */
    public final double lG;

    /** Normalizing constants G(k) for k = 0..N, in the common scale. */
    public final double[] G;

    /** Time-stationary node weights w_j(n), indexed [j][n] with n = 0..N. */
    public final double[][] W;

    /** Complement constants of the cycle with node j removed, indexed [j][k]. */
    public final double[][] Gc;

    /** Arrival node weights v_j(n) = w_j(n) q_j(n) of proposition 3.5. */
    public final double[][] Wa;

    public DpfqnNcLdResult(double lG, double[] G, double[][] W, double[][] Gc, double[][] Wa) {
        this.lG = lG;
        this.G = G;
        this.W = W;
        this.Gc = Gc;
        this.Wa = Wa;
    }

    /** Marginal queue length law of node j, {@code P(X_j = n)} for n = 0..N. */
    public double[] marginal(int j) {
        int N = G.length - 1;
        double[] marg = new double[N + 1];
        for (int n = 0; n <= N; n++) {
            marg[n] = W[j][n] * Gc[j][N - n] / G[N];
        }
        return marg;
    }
}
