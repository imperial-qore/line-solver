/**
 * @file Exact interval-valued MVA for single-class closed product-form networks
 *
 * Exact output intervals of single-class MVA when the service demands, the think time and
 * the population are known only up to intervals (Luthi and Haring, Performance Evaluation
 * 32(3):185-215, 1998). Ported at parity from MATLAB pfqn_mva_interval.m.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Pfqn_mva_interval {
    private Pfqn_mva_interval() {}

    /** Zero tolerance for deciding whether a demand interval is thick. */
    private static final double ZERO_TOL = 1e-14;

    /** Interval-valued mean performance measures. Every field is [lower, upper]. */
    public static final class Result {
        /** Throughput interval (1 x 2). */
        public final Matrix X;
        /** Mean queue-length intervals per station (M x 2). */
        public final Matrix Q;
        /** Utilization enclosures per station (M x 2). */
        public final Matrix U;
        /** Residence-time intervals per station (M x 2). */
        public final Matrix R;
        /** Total response-time interval (1 x 2). */
        public final Matrix Rtot;
        /** Interval of the total number of jobs at the stations (1 x 2). */
        public final Matrix Qtot;

        public Result(Matrix X, Matrix Q, Matrix U, Matrix R, Matrix Rtot, Matrix Qtot) {
            this.X = X;
            this.Q = Q;
            this.U = U;
            this.R = R;
            this.Rtot = Rtot;
            this.Qtot = Qtot;
        }
    }

    /**
     * Exact hull of single-class MVA over an input box.
     *
     * <p>Single-class MVA is monotone in every input: the throughput decreases in each
     * demand and in the think time and increases in the population, the per-station queue
     * length and residence time increase in the own demand and in the population and
     * decrease in the other demands and in the think time, and the totals increase in every
     * demand and in the population and decrease in the think time (Luthi and Haring 1998,
     * Theorems 2-5, Table 1). By their Theorem 1 the exact range of a function monotone in
     * each argument is attained at the endpoints of the input box, so each bound below is
     * one ordinary MVA call at the corner that the sign pattern selects. This is the
     * algorithm of their Fig. 2 and it costs 2*(m+2) MVA calls, m being the number of thick
     * demand intervals; evaluating the MVA recursion in interval arithmetic instead would be
     * a valid but far wider enclosure, since every input recurs at each step (the dependency
     * problem, 14x too wide on the paper's own example).</p>
     *
     * <p>The returned interval is the exact hull of MVA over the input box, not a bound on
     * the true network: it holds conditionally on the demands lying in the box, and says
     * nothing about the accuracy of MVA itself. It must therefore not be composed with the
     * brackets of SolverBA, which bracket the exact solution of a model whose demands are
     * known.</p>
     *
     * <p>Delay stations are folded into Z, exactly as in {@link Pfqn_mva}: a delay demand
     * interval enters as a term of the think-time interval, and the hull of the sum is the
     * sum of the hulls when the delays vary independently. Load-independent single-server
     * queueing stations only, one class only; the monotonicity theorems cover no other
     * case.</p>
     *
     * @param L service demand intervals (M x 2), column 0 lower, column 1 upper
     * @param N population interval (1 x 2, or 1 x 1 for a thin population)
     * @param Z think time interval (1 x 2, or 1 x 1; null for zero)
     * @return interval-valued throughput, queue lengths, utilizations, residence times and
     *         totals
     */
    public static Result pfqn_mva_interval(Matrix L, Matrix N, Matrix Z) {
        if (L == null || L.isEmpty()) {
            throw new IllegalArgumentException("pfqn_mva_interval requires at least one queueing station.");
        }
        Matrix Lint;
        if (L.getNumCols() == 1) {
            Lint = new Matrix(L.getNumRows(), 2);
            for (int i = 0; i < L.getNumRows(); i++) {
                Lint.set(i, 0, L.get(i, 0));
                Lint.set(i, 1, L.get(i, 0));
            }
        } else if (L.getNumCols() == 2) {
            Lint = L;
        } else {
            throw new IllegalArgumentException("pfqn_mva_interval is a single-class method: L must be M x 2, "
                    + "one [lower upper] demand interval per station.");
        }

        double[] nInt = endpoints(N, "population");
        double[] zInt = Z == null || Z.isEmpty() ? new double[]{0.0, 0.0} : endpoints(Z, "think time");

        int M = Lint.getNumRows();
        Matrix Llo = new Matrix(M, 1);
        Matrix Lup = new Matrix(M, 1);
        boolean[] thick = new boolean[M];
        for (int i = 0; i < M; i++) {
            double lo = Lint.get(i, 0);
            double up = Lint.get(i, 1);
            if (lo < 0 || up < 0) {
                throw new IllegalArgumentException("demands and think times must be nonnegative.");
            }
            if (lo > up) {
                throw new IllegalArgumentException("interval lower endpoints must not exceed the upper endpoints.");
            }
            Llo.set(i, 0, lo);
            Lup.set(i, 0, up);
            thick[i] = up > lo + ZERO_TOL;
        }
        if (zInt[0] < 0) {
            throw new IllegalArgumentException("demands and think times must be nonnegative.");
        }
        if (nInt[0] > nInt[1] || zInt[0] > zInt[1]) {
            throw new IllegalArgumentException("interval lower endpoints must not exceed the upper endpoints.");
        }
        if (nInt[0] < 1) {
            throw new IllegalArgumentException("pfqn_mva_interval requires a population interval with at least one "
                    + "job; the monotonicity theorems assume n >= 1.");
        }
        if (Math.abs(nInt[0] - Math.round(nInt[0])) > ZERO_TOL
                || Math.abs(nInt[1] - Math.round(nInt[1])) > ZERO_TOL) {
            throw new IllegalArgumentException("population endpoints must be integers; use Pfqn_nintmva for a "
                    + "nonintegral population.");
        }
        double nlo = Math.round(nInt[0]);
        double nup = Math.round(nInt[1]);
        double zlo = zInt[0];
        double zup = zInt[1];

        Matrix X = new Matrix(1, 2);
        Matrix Q = new Matrix(M, 2);
        Matrix U = new Matrix(M, 2);
        Matrix R = new Matrix(M, 2);
        Matrix Rtot = new Matrix(1, 2);
        Matrix Qtot = new Matrix(1, 2);

        // S1: throughput upper bound, and the upper bounds of the stations whose demand is
        // thin (their own demand is fixed, so lowering the others maximizes them).
        Ret.pfqnMVA s1 = mva(Llo, nup, zlo);
        X.set(0, 1, s1.X.get(0));
        // S2: the same quantities at the opposite corner, giving the lower bounds.
        Ret.pfqnMVA s2 = mva(Lup, nlo, zup);
        X.set(0, 0, s2.X.get(0));
        for (int i = 0; i < M; i++) {
            if (!thick[i]) {
                Q.set(i, 1, s1.Q.get(i));
                R.set(i, 1, s1.R.get(i));
                Q.set(i, 0, s2.Q.get(i));
                R.set(i, 0, s2.R.get(i));
            }
        }

        // S3/S4: the totals increase in the demands and the population and decrease in the
        // think time, so their corners differ from those of the throughput.
        Ret.pfqnMVA s3 = mva(Llo, nlo, zup);
        Rtot.set(0, 0, sum(s3.R));
        Qtot.set(0, 0, sum(s3.Q));
        Ret.pfqnMVA s4 = mva(Lup, nup, zlo);
        Rtot.set(0, 1, sum(s4.R));
        Qtot.set(0, 1, sum(s4.Q));

        // S5/S6: one pair of calls per thick station, its own demand at the endpoint that
        // maximizes (minimizes) it and the others at the opposite endpoint.
        for (int k = 0; k < M; k++) {
            if (!thick[k]) {
                continue;
            }
            Matrix d = Llo.copy();
            d.set(k, 0, Lup.get(k));
            Ret.pfqnMVA s5 = mva(d, nup, zlo);
            Q.set(k, 1, s5.Q.get(k));
            R.set(k, 1, s5.R.get(k));

            d = Lup.copy();
            d.set(k, 0, Llo.get(k));
            Ret.pfqnMVA s6 = mva(d, nlo, zup);
            Q.set(k, 0, s6.Q.get(k));
            R.set(k, 0, s6.R.get(k));
        }

        // U = X*D is not covered by the monotonicity table, so it is enclosed by the product
        // of the two intervals, intersected with the range of a single-server utilization.
        // Where the demand is thin the product is already exact.
        for (int i = 0; i < M; i++) {
            U.set(i, 0, X.get(0, 0) * Llo.get(i));
            U.set(i, 1, Math.min(1.0, X.get(0, 1) * Lup.get(i)));
        }

        return new Result(X, Q, U, R, Rtot, Qtot);
    }

    /** Convenience overload for a thin population and think time. */
    public static Result pfqn_mva_interval(Matrix L, double N, double Z) {
        Matrix Nm = new Matrix(1, 1);
        Nm.set(0, 0, N);
        Matrix Zm = new Matrix(1, 1);
        Zm.set(0, 0, Z);
        return pfqn_mva_interval(L, Nm, Zm);
    }

    private static Ret.pfqnMVA mva(Matrix L, double n, double z) {
        Matrix N = new Matrix(1, 1);
        N.set(0, 0, n);
        Matrix Z = new Matrix(1, 1);
        Z.set(0, 0, z);
        return Pfqn_mva.pfqn_mva(L, N, Z, null);
    }

    private static double[] endpoints(Matrix v, String what) {
        if (v == null || v.isEmpty()) {
            throw new IllegalArgumentException(what + " must be given as a scalar or as a [lower upper] interval.");
        }
        int n = v.length();
        if (n == 1) {
            return new double[]{v.get(0), v.get(0)};
        }
        if (n != 2) {
            throw new IllegalArgumentException(what + " must be given as a scalar or as a [lower upper] interval.");
        }
        return new double[]{v.get(0), v.get(1)};
    }

    private static double sum(Matrix m) {
        double s = 0.0;
        for (int i = 0; i < m.length(); i++) {
            s += m.get(i);
        }
        return s;
    }
}
