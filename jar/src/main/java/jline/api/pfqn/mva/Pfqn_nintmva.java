/**
 * @file Mean value analysis at a nonintegral population (fractional-base aMVA)
 *
 * Exact MVA recursion started from the FRACTIONAL base n0 = N - floor(N) instead of from
 * the empty network, giving mean performance measures at a real-valued population
 * (Dowdy and Gordon 1984, "aMVA"). Ported at parity from MATLAB pfqn_nintmva.m.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.util.matrix.Matrix;

public final class Pfqn_nintmva {
    private Pfqn_nintmva() {}

    /** Mean performance measures at a real-valued population. */
    public static final class Result {
        /** Throughput. */
        public final double X;
        /** Mean queue lengths (M x 1). */
        public final Matrix Q;
        /** Utilizations (M x 1). */
        public final Matrix U;
        /** Residence times (M x 1). */
        public final Matrix R;

        public Result(double X, Matrix Q, Matrix U, Matrix R) {
            this.X = X;
            this.Q = Q;
            this.U = U;
            this.R = R;
        }
    }

    /**
     * MVA recursion from the fractional base.
     *
     * <p>The recursion is the standard Reiser-Lavenberg one, stepped in unit increments
     * from n = N - floor(N) (where the arrival-theorem term is taken as 0, the network
     * below the base being empty) up to n = N. At integer N the base is 0 and the
     * recursion is bit-identical to exact MVA; at fractional N it interpolates smoothly
     * through the integral points, which is what a nonintegral degree of multiprogramming
     * (a time-average over a measurement window) calls for.</p>
     *
     * <p>Unlike {@link jline.api.pfqn.nc.Pfqn_dnc} this accepts a think time. It is
     * single-class: for fractional multiclass populations use Pfqn_bs, which accepts them
     * directly.</p>
     *
     * @param L service demand vector (M x 1) of the queueing stations
     * @param N population (real, nonnegative; may be fractional)
     * @param Z think time
     * @return throughput, queue lengths, utilizations and residence times
     */
    public static Result pfqn_nintmva(Matrix L, double N, double Z) {
        if (N < 0) {
            throw new IllegalArgumentException("pfqn_nintmva requires a nonnegative population.");
        }
        int M = L.length();
        Matrix Q = new Matrix(M, 1);
        Matrix R = new Matrix(M, 1);
        Matrix U = new Matrix(M, 1);
        double X = 0.0;
        if (N == 0) {
            return new Result(0.0, Q, U, R);
        }

        double n = N - Math.floor(N);
        if (n == 0) {
            n = 1.0;
        }
        while (n <= N + 1e-12) {
            double sumR = 0.0;
            for (int i = 0; i < M; i++) {
                double ri = L.get(i) * (1 + Q.get(i));
                R.set(i, 0, ri);
                sumR += ri;
            }
            X = n / (Z + sumR);
            for (int i = 0; i < M; i++) {
                Q.set(i, 0, X * R.get(i));
            }
            n += 1;
        }
        for (int i = 0; i < M; i++) {
            U.set(i, 0, X * L.get(i));
        }
        return new Result(X, Q, U, R);
    }
}
