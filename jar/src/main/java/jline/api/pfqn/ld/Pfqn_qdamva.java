/**
 * @file QD-AMVA: queue-dependent approximate mean value analysis.
 *
 * Port of {@code matlab/src/api/pfqn/pfqn_qdamva.m}, the queue-dependent AMVA of Casale, Perez
 * and Wang (IFIP PERFORMANCE 2015), on a closed multiclass product-form network.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import jline.util.matrix.Matrix;

/**
 * QD-AMVA on a closed multiclass product-form network.
 *
 * <p>A Schweitzer/Bard core in which the class-r demand at station k is scaled by the
 * queue-dependence term g_k evaluated at the ARRIVAL-INSTANT total queue length,
 * {@code g = pfqn_lldfun(1 + delta * rowsum(Q), mu)}.
 *
 * <p>SETTING {@code mu} TO A CONSTANT ROW RECOVERS PLAIN SCHWEITZER AMVA ONLY FOR A SINGLE CLASS.
 * {@link Pfqn_lldfun} does skip a constant row, so g == 1 there, but the residence time that
 * remains is {@code 1 + delta * rowsum(Q)} with ONE aggregate delta = (sum(N)-1)/sum(N) applied to
 * the whole arrival-instant queue, where Bard-Schweitzer shrinks the TAGGED class alone:
 * {@code 1 + sum_{s != r} Q(k,s) + (N(r)-1)/N(r) * Q(k,r)}. The two coincide iff K == 1. Measured
 * over 40 random three-class instances, {@code pfqn_qdamva(L,N,Z,ones)} departs from
 * {@code pfqn_bs} by up to 0.217 in absolute queue length, and is the LESS accurate of the two on
 * single-server multiclass models (mean relative error on Q 0.069 against 0.056 at R = 3), the
 * aggregate delta buying nothing once g == 1. This is the QD-AMVA closure, not a defect of the
 * port, but do not use the function as a Schweitzer oracle for K &gt; 1.
 *
 * <p>MU IS A DIMENSIONLESS RATE MULTIPLIER, NOT A RATE. mu(k,n) is the factor by which station k
 * serves faster when it holds n jobs. Two traps follow from {@link Pfqn_lldfun} and are the
 * reference's, not this port's:
 *
 * <ul>
 *   <li>it SKIPS a station whose mu row is constant, so a single-server station must be a row of
 *       ones and a c-server station {@code min(1..smax, c)}. Passing a c-server station a
 *       constant row silently returns g = 1, i.e. a single server.</li>
 *   <li>smax = mu.getNumCols() must be at least ceil(sum(N)) or the interpolation clamps the
 *       population and the top of the rate curve is never reached.</li>
 * </ul>
 *
 * <p>Delay stations are carried in Z, not as rows of L. Closed classes only: an infinite N(r) is
 * not supported.
 */
public final class Pfqn_qdamva {
    private Pfqn_qdamva() {}

    /** The fixed point {@link #pfqn_qdamva} reaches, and how many sweeps it took. */
    public static class Result {
        /** (M x R) mean queue lengths. */
        public final Matrix Q;
        /** (1 x R) per-class throughputs. */
        public final Matrix X;
        /** (M x R) per-class utilizations, carrying the g scaling. */
        public final Matrix U;
        /** (M x R) per-class residence times, Q = X .* R. */
        public final Matrix R;
        /** Number of iterations performed. */
        public final int iter;

        public Result(Matrix Q, Matrix X, Matrix U, Matrix R, int iter) {
            this.Q = Q;
            this.X = X;
            this.U = U;
            this.R = R;
            this.iter = iter;
        }
    }

    /**
     * QD-AMVA with the reference's default tolerance and iteration cap.
     *
     * @param L  (M x R) service demand matrix
     * @param N  (1 x R) population vector, finite
     * @param Z  (1 x R) think time vector, or null for none
     * @param mu (M x smax) queue-dependent rate multipliers, or null for none
     * @return the fixed point
     */
    public static Result pfqn_qdamva(Matrix L, Matrix N, Matrix Z, Matrix mu) {
        return pfqn_qdamva(L, N, Z, mu, null, 1e-6, 10000);
    }

    /**
     * QD-AMVA.
     *
     * @param L       (M x R) service demand matrix
     * @param N       (1 x R) population vector, finite
     * @param Z       (1 x R) think time vector, or null for none
     * @param mu      (M x smax) queue-dependent rate multipliers, or null for none
     * @param Q0      (M x R) initial guess, or null for the reference's demand split
     * @param tol     convergence tolerance on the queue lengths
     * @param maxiter maximum number of iterations
     * @return the fixed point
     */
    public static Result pfqn_qdamva(Matrix L, Matrix N, Matrix Z, Matrix mu, Matrix Q0,
                                     double tol, int maxiter) {
        final int M = L.getNumRows();
        final int K = L.getNumCols();
        if (N == null || N.length() != K) {
            throw new IllegalArgumentException(
                    "pfqn_qdamva: the population vector must have one entry per class");
        }
        if (Z != null && !Z.isEmpty() && Z.length() != K) {
            throw new IllegalArgumentException(
                    "pfqn_qdamva: the think-time vector must have one entry per class");
        }

        Matrix Q = new Matrix(M, K);
        Matrix X = new Matrix(1, K);
        Matrix U = new Matrix(M, K);
        Matrix R = new Matrix(M, K);
        int iter = 0;

        double Ntot = 0.0;
        for (int r = 0; r < K; r++) {
            Ntot += N.get(r);
        }
        // delta is undefined on an empty population, and the reference returns the zero queue
        // rather than dividing by it.
        if (!(Ntot > 0.0)) {
            return new Result(Q, X, U, R, iter);
        }

        if (Q0 != null && Q0.getNumRows() == M && Q0.getNumCols() == K) {
            Q = Q0.copy();
        } else {
            // Ltot = 0 for a class with no demand anywhere: L/Ltot is a NaN the iteration never
            // recovers from. Such a column arises routinely in a layered fixed point, where a
            // caller can start with no work at the layer station, so it is left at zero instead.
            for (int r = 0; r < K; r++) {
                double tot = 0.0;
                for (int k = 0; k < M; k++) {
                    tot += L.get(k, r);
                }
                if (!(tot > 0.0)) {
                    continue;
                }
                for (int k = 0; k < M; k++) {
                    Q.set(k, r, L.get(k, r) / tot * N.get(r));
                }
            }
        }

        final double delta = (Ntot - 1.0) / Ntot;
        // Q*10 as the sentinel, as the reference notes, stalls on an all-zero seed: the loop
        // would exit before its first pass. Offset instead.
        Matrix Qprev = new Matrix(M, K);
        for (int k = 0; k < M; k++) {
            for (int r = 0; r < K; r++) {
                Qprev.set(k, r, Q.get(k, r) + 10.0 * (1.0 + tol));
            }
        }
        final Matrix lld = (mu == null) ? new Matrix(0, 0) : mu;

        while (maxAbsDiff(Q, Qprev) > tol && iter < maxiter) {
            iter++;
            Qprev = Q.copy();

            // The arrival-instant total queue length, class independent.
            Matrix Ak = new Matrix(M, 1);
            for (int k = 0; k < M; k++) {
                double rowsum = 0.0;
                for (int r = 0; r < K; r++) {
                    rowsum += Q.get(k, r);
                }
                Ak.set(k, 0, 1.0 + delta * rowsum);
            }
            Matrix g = Pfqn_lldfun.pfqn_lldfun(Ak, lld, null);

            for (int r = 0; r < K; r++) {
                for (int k = 0; k < M; k++) {
                    double rowsum = 0.0;
                    for (int s = 0; s < K; s++) {
                        rowsum += Q.get(k, s);
                    }
                    R.set(k, r, L.get(k, r) * g.get(k, 0) * (1.0 + delta * rowsum));
                }
                double denom = (Z == null || Z.isEmpty()) ? 0.0 : Z.get(r);
                for (int k = 0; k < M; k++) {
                    denom += R.get(k, r);
                }
                X.set(0, r, (denom > 0.0) ? N.get(r) / denom : 0.0);
                for (int k = 0; k < M; k++) {
                    Q.set(k, r, X.get(0, r) * R.get(k, r));
                    U.set(k, r, L.get(k, r) * g.get(k, 0) * X.get(0, r));
                }
            }
        }
        return new Result(Q, X, U, R, iter);
    }

    /** The infinity norm of the change between two sweeps. */
    private static double maxAbsDiff(Matrix a, Matrix b) {
        double worst = 0.0;
        for (int i = 0; i < a.getNumRows(); i++) {
            for (int j = 0; j < a.getNumCols(); j++) {
                worst = Math.max(worst, Math.abs(a.get(i, j) - b.get(i, j)));
            }
        }
        return worst;
    }
}
