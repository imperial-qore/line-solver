/**
 * @file Operational sensitivity of throughput to homogeneous-service-time violations
 *
 * Robustness certificate for a single-class closed product-form solution: how far the
 * predicted throughput can move when the homogeneous-service-time (HST) assumption fails
 * at one station (Suri 1983, eq. 3.11, Lemma 3.1, problem (P1)). Ported at parity from
 * MATLAB pfqn_hst.m.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.sens;

import org.apache.commons.math3.util.FastMath;

import jline.api.pfqn.nc.Pfqn_rgf;
import jline.util.matrix.Matrix;

public final class Pfqn_hst {
    private Pfqn_hst() {}

    /** HST sensitivity coefficients and the constrained worst case. */
    public static final class Result {
        /** Station index analysed (0-based). */
        public final int station;
        /** Product-form throughput. */
        public final double X;
        /** Utilization of that station. */
        public final double U;
        /** Mean queue length there. */
        public final double Q;
        /** P(n_i &gt;= k), k = 0..N. */
        public final double[] Pgeq;
        /** P(n_i = k), k = 0..N. */
        public final double[] p;
        /** Sensitivity coefficients c_k, k = 1..N (eq. 3.11). */
        public final double[] c;
        /** sum_k |c_k|, the unconstrained certificate per unit d. */
        public final double total;
        /** The (P1) optimum per unit d. */
        public final double worst;
        /** The worst-case deviation profile a_k/d, k = 1..N. */
        public final double[] astar;

        public Result(int station, double X, double U, double Q, double[] Pgeq, double[] p,
                      double[] c, double total, double worst, double[] astar) {
            this.station = station;
            this.X = X;
            this.U = U;
            this.Q = Q;
            this.Pgeq = Pgeq;
            this.p = p;
            this.c = c;
            this.total = total;
            this.worst = worst;
            this.astar = astar;
        }
    }

    /**
     * HST robustness certificate.
     *
     * <p>HST states that the mean service time at station i does not depend on the queue
     * length there. Suri perturbs it to S_i(n) = S_i (1 + a_n), one relative deviation per
     * queue-length level, and shows (eq. 3.11) that to first order
     * c_n = P(n_i &gt;= n+1)/u_i - P(n_i &gt;= n), with u_i = L_i X0 and
     * P(n_i &gt;= n) = L_i^n G(N-n)/G(N). The naive certificate
     * |dX0/X0| &lt;= (sum_n |c_n|) d follows from |a_n| &lt;= d alone, and by Lemma 3.1
     * that total equals Q_i(N) - Q_i(N-1).</p>
     *
     * <p>That bound is loose because an operationally consistent perturbation must leave
     * the observed mean service time unchanged, sum_n p_n a_n = 0. The constrained problem
     * (P1) max |sum_n c_n a_n| subject to |a_n| &lt;= d and sum_n p_n a_n = 0 is a
     * one-constraint linear program, solved here exactly: its optimum sets a_n = +/-d
     * according to whether c_n/p_n exceeds a threshold, with at most one fractional
     * coordinate. The feasible set is symmetric, so the optimum comes in a +/- pair and
     * only the magnitude is determined.</p>
     *
     * @param L   service demand vector (M x 1) of the queueing stations
     * @param N   population (integer, at least one job)
     * @param Z   think time
     * @param ist 0-based station index the perturbation applies to
     * @return the marginals, the sensitivity coefficients and the (P1) optimum
     */
    public static Result pfqn_hst(Matrix L, double N, double Z, int ist) {
        int M = L.length();
        if (N < 1 || N != FastMath.rint(N)) {
            throw new IllegalArgumentException("pfqn_hst requires an integer population of at least one job.");
        }
        int n = (int) FastMath.rint(N);
        if (ist < 0 || ist >= M) {
            throw new IllegalArgumentException("Station index " + ist + " is out of range (the model has "
                    + M + " queueing stations).");
        }
        if (L.get(ist) <= 0) {
            throw new IllegalArgumentException("Station " + ist
                    + " has zero demand, so its queue-length marginals are degenerate.");
        }

        double[] lg = Pfqn_rgf.pfqn_rgf(L, n, Z).lg;      // lg[k] = log G(k)
        double X = FastMath.exp(lg[n - 1] - lg[n]);
        double y = L.get(ist);

        double[] Pgeq = new double[n + 1];
        for (int k = 0; k <= n; k++) {
            Pgeq[k] = FastMath.exp(k * FastMath.log(y) + lg[n - k] - lg[n]);
        }
        double[] p = new double[n + 1];
        for (int k = 0; k <= n; k++) {
            p[k] = Pgeq[k] - (k + 1 <= n ? Pgeq[k + 1] : 0.0);
        }
        double U = y * X;
        double Q = 0.0;
        for (int k = 1; k <= n; k++) {
            Q += Pgeq[k];
        }

        // eq. (3.11)
        double[] c = new double[n];
        double total = 0.0;
        for (int i = 1; i <= n; i++) {
            double pn1 = (i + 1 <= n) ? Pgeq[i + 1] : 0.0;
            c[i - 1] = pn1 / U - Pgeq[i];
            total += FastMath.abs(c[i - 1]);
        }

        // (P1): one equality constraint plus a box. At the optimum
        // a_n = sign(c_n - lambda p_n), so sorting by c_n/p_n and sweeping the split point
        // enumerates every candidate lambda; the constraint fixes the single fractional
        // coordinate at the split.
        double[] pp = new double[n];
        System.arraycopy(p, 1, pp, 0, n);
        Integer[] order = new Integer[n];
        final double[] ratio = new double[n];
        for (int i = 0; i < n; i++) {
            order[i] = i;
            ratio[i] = c[i] / FastMath.max(pp[i], Double.MIN_VALUE);
        }
        java.util.Arrays.sort(order, new java.util.Comparator<Integer>() {
            public int compare(Integer a, Integer b) {
                return Double.compare(ratio[b], ratio[a]);
            }
        });

        double best = 0.0;
        double[] astar = new double[n];
        double[] a = new double[n];
        for (int split = 0; split <= n; split++) {
            java.util.Arrays.fill(a, -1.0);
            for (int q = 0; q < split; q++) {
                a[order[q]] = 1.0;
            }
            for (int piv = 0; piv < n; piv++) {
                if (pp[piv] <= 0) {
                    continue;
                }
                double rest = 0.0;
                for (int i = 0; i < n; i++) {
                    if (i != piv) {
                        rest += pp[i] * a[i];
                    }
                }
                double v = -rest / pp[piv];
                if (v < -1 || v > 1) {
                    continue;
                }
                double old = a[piv];
                a[piv] = v;
                double obj = 0.0;
                for (int i = 0; i < n; i++) {
                    obj += c[i] * a[i];
                }
                if (FastMath.abs(obj) > FastMath.abs(best)) {
                    best = obj;
                    astar = a.clone();
                }
                a[piv] = old;
            }
        }
        if (best < 0) {
            best = -best;
            for (int i = 0; i < n; i++) {
                astar[i] = -astar[i];   // the feasible set is symmetric
            }
        }

        return new Result(ist, X, U, Q, Pgeq, p, c, total, best, astar);
    }

    /** Convenience overload analysing the bottleneck station. */
    public static Result pfqn_hst(Matrix L, double N, double Z) {
        int best = 0;
        for (int i = 1; i < L.length(); i++) {
            if (L.get(i) > L.get(best)) {
                best = i;
            }
        }
        return pfqn_hst(L, N, Z, best);
    }
}
