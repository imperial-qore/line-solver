/**
 * @file Closed-network fork-join bounds and mean value analysis
 *
 * Ports of matlab/src/api/fj/fj_qgb.m, fj_amva.m and fj_respt_closed.m, the
 * closed-network group of A. Thomasian, "Analysis of Fork/Join and Related
 * Queueing Systems", ACM Computing Surveys 47(2), Article 17, 2014
 * (Eqs. (67)-(70)).
 *
 * @since LINE 3.0
 */
package jline.api.fj;

public final class FJ_closed {
    private FJ_closed() {}

    /** [Q, y] of fj_qgb: the bounded queue lengths and the geometric ratios. */
    public static final class FJQgbResult {
        public final double[] Q;
        public final double[] y;

        public FJQgbResult(double[] Q, double[] y) {
            this.Q = Q;
            this.y = y;
        }
    }

    /**
     * Geometric bound on the queue length of each fork-join subnetwork.
     *
     * y_n(M) = D_n M / (Z + sum_j D_j H_{P_j} + Dmax M),
     * Q_n(M) = H_{P_n} [ y_n/(1-y_n) - y_n^(M+1)/(1-y_n) ].
     *
     * The harmonic weights are what distinguishes this from the ordinary
     * geometric bound of Pfqn_qzgblow: a P-way fork-join subnetwork inflates its
     * own demand by H_P in the denominator and its queue length by H_P in the
     * numerator. Setting every fork degree to one recovers that bound exactly.
     */
    public static FJQgbResult fj_qgb(double[] D, int[] P, int M, double Z) {
        if (D.length != P.length) {
            throw new IllegalArgumentException("D and P must have the same number of elements.");
        }
        if (D.length == 0) {
            throw new IllegalArgumentException("At least one subnetwork is required.");
        }
        if (M < 1) {
            throw new IllegalArgumentException("M must be a positive integer. Got M=" + M + ".");
        }
        if (Z < 0) {
            throw new IllegalArgumentException("The think time must be non-negative.");
        }
        int N = D.length;
        double[] H = new double[N];
        double Dtot = 0, Dmax = 0;
        for (int n = 0; n < N; n++) {
            if (D[n] < 0) {
                throw new IllegalArgumentException("Service demands must be non-negative.");
            }
            if (P[n] < 1) {
                throw new IllegalArgumentException("Fork degrees must be positive integers.");
            }
            H[n] = FJ_harmonic.fj_harmonic(P[n]);
            Dtot += D[n] * H[n];
            if (n == 0 || D[n] > Dmax) {
                Dmax = D[n];
            }
        }
        double[] Q = new double[N];
        double[] y = new double[N];
        for (int n = 0; n < N; n++) {
            y[n] = D[n] * M / (Z + Dtot + Dmax * M);
            if (y[n] < 1.0) {
                Q[n] = H[n] * (y[n] / (1 - y[n]) - Math.pow(y[n], M + 1) / (1 - y[n]));
            } else {
                // Degenerate ratio: the bound collapses onto the full population
                Q[n] = M;
            }
        }
        return new FJQgbResult(Q, y);
    }

    public static FJQgbResult fj_qgb(double[] D, int[] P, int M) {
        return fj_qgb(D, P, M, 0.0);
    }

    /** [R, Q, X, U] of fj_amva. */
    public static final class FJAmvaResult {
        public final double[] R;
        public final double[] Q;
        public final double X;
        public final double[] U;

        public FJAmvaResult(double[] R, double[] Q, double X, double[] U) {
            this.R = R;
            this.Q = Q;
            this.X = X;
            this.U = U;
        }
    }

    /**
     * Population-by-population mean value analysis of a closed network of
     * fork-join subnetworks.
     *
     * R_n(m) = D_n [ H_{P_n} + Q_n(m-1) ], X(m) = m / (Z + sum_n R_n(m)),
     * Q_n(m) = X(m) R_n(m), started from Q_n(0) = 0.
     *
     * With every fork degree equal to one this is the exact single-class mean
     * value analysis, because H_1 = 1; above that it is an approximation whose
     * per-subnetwork residence time is an upper bound in the sense of Varki.
     */
    public static FJAmvaResult fj_amva(double[] D, int[] P, int M, double Z) {
        if (D.length != P.length) {
            throw new IllegalArgumentException("D and P must have the same number of elements.");
        }
        if (D.length == 0) {
            throw new IllegalArgumentException("At least one subnetwork is required.");
        }
        if (M < 1) {
            throw new IllegalArgumentException("M must be a positive integer. Got M=" + M + ".");
        }
        if (Z < 0) {
            throw new IllegalArgumentException("The think time must be non-negative.");
        }
        int N = D.length;
        double[] H = new double[N];
        for (int n = 0; n < N; n++) {
            if (D[n] < 0) {
                throw new IllegalArgumentException("Service demands must be non-negative.");
            }
            if (P[n] < 1) {
                throw new IllegalArgumentException("Fork degrees must be positive integers.");
            }
            H[n] = FJ_harmonic.fj_harmonic(P[n]);
        }
        double[] Q = new double[N];
        double[] R = new double[N];
        double[] U = new double[N];
        double X = 0;
        for (int m = 1; m <= M; m++) {
            double Rtot = 0;
            for (int n = 0; n < N; n++) {
                R[n] = D[n] * (H[n] + Q[n]);
                Rtot += R[n];
            }
            if (!(Rtot > 0)) {
                throw new IllegalArgumentException(
                        "The total residence time vanished; every demand is zero.");
            }
            X = m / (Z + Rtot);
            for (int n = 0; n < N; n++) {
                Q[n] = X * R[n];
            }
        }
        // Each subnetwork holds P(n) queues sharing the demand equally
        for (int n = 0; n < N; n++) {
            U[n] = X * D[n] / P[n];
        }
        return new FJAmvaResult(R, Q, X, U);
    }

    public static FJAmvaResult fj_amva(double[] D, int[] P, int M) {
        return fj_amva(D, P, M, 0.0);
    }

    /** [R, exact] of fj_respt_closed. */
    public static final class FJResptClosedResult {
        public final double R;
        public final boolean exact;

        public FJResptClosedResult(double R, boolean exact) {
            this.R = R;
            this.exact = exact;
        }
    }

    /**
     * Varki bound on the residence time of a closed fork-join subnetwork,
     * R_{P_K}(M) &lt;= x [ H_K + A ], with A the mean number of jobs an arriving
     * job finds at the subnetwork.
     */
    public static FJResptClosedResult fj_respt_closed(int K, double x, int M, double A) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }
        if (!(x > 0)) {
            throw new IllegalArgumentException("The mean service time must be positive.");
        }
        if (M < 1) {
            throw new IllegalArgumentException("M must be a positive integer. Got M=" + M + ".");
        }
        if (A < 0) {
            throw new IllegalArgumentException(
                    "The arrival-instant queue length must be non-negative.");
        }
        return new FJResptClosedResult(x * (FJ_harmonic.fj_harmonic(K) + A), false);
    }

    /**
     * The isolated parallel subsystem of Theorem 4.1, where every other job is
     * necessarily inside the subnetwork so A = M-1. The bound holds with
     * equality at K = 2.
     */
    public static FJResptClosedResult fj_respt_closed(int K, double x, int M) {
        if (M < 1) {
            throw new IllegalArgumentException("M must be a positive integer. Got M=" + M + ".");
        }
        FJResptClosedResult r = fj_respt_closed(K, x, M, (double) (M - 1));
        return new FJResptClosedResult(r.R, K == 2);
    }
}
