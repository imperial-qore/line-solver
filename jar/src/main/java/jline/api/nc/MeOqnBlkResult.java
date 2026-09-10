/**
 * @file Maximum Entropy OQN with blocking: algorithm result
 *
 * Result of the ME algorithm for open queueing networks with finite
 * buffers, loss and transfer blocking.
 *
 * @since LINE 3.0
 */
package jline.api.nc;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Result of the ME OQN algorithm with finite buffers.
 */
public final class MeOqnBlkResult {
    private final Matrix Q;
    private final Matrix W;
    private final Matrix T;
    private final Matrix U;
    private final Matrix Ca;
    private final Matrix Cd;
    private final Matrix PBa;
    private final Matrix lambda;
    private final int iter;

    public MeOqnBlkResult(Matrix Q, Matrix W, Matrix T, Matrix U, Matrix Ca, Matrix Cd,
                          Matrix PBa, Matrix lambda, int iter) {
        this.Q = Q;
        this.W = W;
        this.T = T;
        this.U = U;
        this.Ca = Ca;
        this.Cd = Cd;
        this.PBa = PBa;
        this.lambda = lambda;
        this.iter = iter;
    }

    /** Mean number of jobs at each station, the jobs held blocked included. */
    public Matrix getQ() { return Q; }
    /** Mean response time at each station. */
    public Matrix getW() { return W; }
    /** Carried throughput at each station. */
    public Matrix getT() { return T; }
    /** Utilization at each station, a blocked server not counted as busy. */
    public Matrix getU() { return U; }
    /** Interarrival scv of the offered flow at each station. */
    public Matrix getCa() { return Ca; }
    /** Interdeparture scv at each station. */
    public Matrix getCd() { return Cd; }
    /** Probability that an arrival finds the station full. */
    public Matrix getPBa() { return PBa; }
    /** Offered arrival rate at each station, re-attempts included. */
    public Matrix getLambda() { return lambda; }
    /** Number of fixed-point iterations performed. */
    public int getIter() { return iter; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof MeOqnBlkResult)) return false;
        MeOqnBlkResult that = (MeOqnBlkResult) o;
        return iter == that.iter
                && Objects.equals(Q, that.Q)
                && Objects.equals(W, that.W)
                && Objects.equals(T, that.T)
                && Objects.equals(U, that.U)
                && Objects.equals(Ca, that.Ca)
                && Objects.equals(Cd, that.Cd)
                && Objects.equals(PBa, that.PBa)
                && Objects.equals(lambda, that.lambda);
    }

    @Override
    public int hashCode() {
        return Objects.hash(Q, W, T, U, Ca, Cd, PBa, lambda, iter);
    }
}
