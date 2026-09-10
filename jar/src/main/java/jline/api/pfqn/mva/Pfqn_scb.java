/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * @file Dowdy-Carlson-Krantz-Tripathi (1992) single-class bounds of multi-class
 *       queueing networks
 *
 * L. W. Dowdy, B. M. Carlson, A. T. Krantz, S. K. Tripathi, "Single-Class Bounds of
 * Multi-Class Queuing Networks", J. ACM 39(1):188-213, 1992. Ported at parity from
 * MATLAB pfqn_scb.m / pfqn_scbgap.m / pfqn_usumbound.m / pfqn_minclasses.m.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.util.matrix.Matrix;

public final class Pfqn_scb {
    private Pfqn_scb() {}

    /**
     * Bracket on the total throughput and the per-device utilizations of the UNKNOWN
     * multiclass system whose single-class counterpart has demand vector L at population N.
     *
     * <p>SEMANTICS DIFFER FROM EVERY OTHER pfqn_* BOUND. aba/bjb/gb/... bracket the exact
     * solution OF THE GIVEN MODEL; this brackets the multiclass system that the given
     * single-class model aggregates. The lower side is therefore the EXACT single-class
     * solution, not an approximation of it.
     *
     * <p>Theorem 2 / Corollary 2: aggregating an R-class model into its single-class
     * counterpart can only understate performance, U_k,1 &lt;= U_k,R and X_1 &lt;= X_R, and
     * Corollary 1 makes the utilization ratio uniform, U_k,R/U_k,1 = X_R/X_1 for every k.
     * Theorem 3 (their Expression 3) caps the relative throughput error at (m-1)/(N+m-1),
     * m = min(N,K), independently of the demands. The single-server capacity U_k,R &lt;= 1
     * caps the same ratio at 1/(X_1*max(L)), tight on the paper's own worst case, so both
     * are applied.
     *
     * @param L service demand vector of the single-class model (K x 1), queueing stations
     *          only; delay stations are not admitted because Theorem 3 rests on the
     *          delay-free balanced-network throughput
     * @param N total population (N &gt;= 1)
     * @return {Xlo, Xhi, Ulo(0..K-1), Uhi(0..K-1)} flattened as a length-(2+2K) array
     */
    public static double[] pfqn_scb(Matrix L, int N) {
        int K = L.length();
        if (K == 0) {
            throw new RuntimeException("pfqn_scb requires at least one queueing station.");
        }
        if (N < 1) {
            throw new RuntimeException("pfqn_scb requires N >= 1.");
        }
        // Exact single-class MVA at Z=0. This IS the lower bound (Theorem 2), so it is
        // computed exactly rather than bounded: a bounded X1 would not bracket X_R.
        double[] Q = new double[K];
        double[] Rk = new double[K];
        double X1 = 0.0;
        for (int n = 1; n <= N; n++) {
            double sumR = 0.0;
            for (int k = 0; k < K; k++) {
                Rk[k] = L.get(k) * (1.0 + Q[k]);
                sumR += Rk[k];
            }
            X1 = n / sumR;
            for (int k = 0; k < K; k++) {
                Q[k] = X1 * Rk[k];
            }
        }
        int m = Math.min(N, K);
        double ratio = ((double) (N + m - 1)) / N;      // Theorem 3, Expression (3)
        double Dmax = L.elementMax();
        if (X1 * Dmax > 0) {
            // U_k,R <= 1 with the uniform ratio of Corollary 1. Tight at the worst case.
            ratio = Math.min(ratio, 1.0 / (X1 * Dmax));
        }
        double[] out = new double[2 + 2 * K];
        out[0] = X1;
        out[1] = X1 * ratio;
        for (int k = 0; k < K; k++) {
            out[2 + k] = X1 * L.get(k);
            out[2 + K + k] = X1 * L.get(k) * ratio;
        }
        return out;
    }

    /**
     * Demand-free bound on the relative throughput error incurred when r of the N
     * single-customer classes of a closed product-form network are merged into one class.
     * With r = N this is the full single-class aggregation error of their Theorem 3, at
     * most 50%; with r &lt; N it is the partial-aggregation error of their Theorem 4.
     *
     * <p>General case, dominating classes allowed (Expression 4, and with r = N Expression
     * 3): e = (min(r,K)-1)/(r+min(r,K)-1). Undominated case, every customer placing the
     * same total demand (Theorem 5 and its comment (3), which lifts the N = R restriction):
     * e = r(r-1)/(min(N,K)(2r-1)), valid for r &lt;= K only, smaller than the general case by
     * the factor r/min(N,K) and equal to it at r = K. THE DOMAIN IS NOT COSMETIC: Theorem 5
     * gives each of its R classes a dedicated device, so r never exceeds K there, and
     * comment (3) states the generalization for r &lt; K. Evaluated at r &gt; K the expression
     * climbs past the general bound and past the 50% cap of Theorem 3, i.e. it stops being a
     * bound, so r &gt; K is refused rather than returned.
     *
     * @param N            total number of customers; one per class, so also the number of
     *                     classes before merging
     * @param K            number of queueing devices
     * @param r            number of classes merged into one (1 &lt;= r &lt;= N)
     * @param undominated  true for the tighter Theorem-5 form, valid only when no class
     *                     dominates, i.e. every customer's total device demand is equal, and
     *                     only for r &lt;= K
     * @return maximum relative throughput error, in [0,1/2]
     */
    public static double pfqn_scbgap(int N, int K, int r, boolean undominated) {
        if (N < 1 || K < 1) {
            throw new RuntimeException("pfqn_scbgap requires N >= 1 and K >= 1.");
        }
        if (r < 1 || r > N) {
            throw new RuntimeException("pfqn_scbgap requires 1 <= r <= N (r=" + r + ", N=" + N + ").");
        }
        if (r == 1) {
            return 0.0;                                     // merging one class changes nothing
        }
        if (undominated) {
            if (r > K) {
                throw new RuntimeException("The undominated (Theorem 5) form is defined for r <= K only (r="
                        + r + ", K=" + K + "); beyond it the expression exceeds the general bound and the 50% cap.");
            }
            return ((double) r * (r - 1)) / (Math.min(N, K) * (2.0 * r - 1));   // Thm 5, comment (3)
        }
        int m = Math.min(r, K);
        return ((double) (m - 1)) / (r + m - 1);            // Expression (4); r=N gives (3)
    }

    /** Full single-class aggregation, dominating classes allowed: r = N, undominated = false. */
    public static double pfqn_scbgap(int N, int K) {
        return pfqn_scbgap(N, K, N, false);
    }

    /**
     * Largest value the sum of device utilizations sum_k U_k,R can take in any closed
     * product-form network with R classes, K devices and N customers (their Theorem 6):
     * sum_k U_k,R &lt;= (H-1) + (K-H+1)(N-H+1)/(K+N-2H+1), H = min(R,K). The bound is
     * demand-free and nondecreasing in R, which is what makes it invertible into a lower
     * bound on the number of necessary classes; see {@link #pfqn_minclasses}.
     *
     * @param R number of single-customer classes (1 &lt;= R &lt;= N)
     * @param K number of devices
     * @param N total number of customers
     * @return upper bound on sum_k U_k,R
     */
    public static double pfqn_usumbound(int R, int K, int N) {
        if (N < 1 || K < 1) {
            throw new RuntimeException("pfqn_usumbound requires N >= 1 and K >= 1.");
        }
        if (R < 1 || R > N) {
            throw new RuntimeException("pfqn_usumbound requires 1 <= R <= N (R=" + R + ", N=" + N + ").");
        }
        int H = Math.min(R, K);
        return (H - 1) + ((double) (K - H + 1) * (N - H + 1)) / (K + N - 2 * H + 1);
    }

    /**
     * Smallest number of customer classes R consistent with an observed sum of device
     * utilizations, obtained by inverting the demand-free Expression (6) bound of
     * {@link #pfqn_usumbound}, which is nondecreasing in R. Only measured quantities are
     * needed, so the answer is available BEFORE any class-specific demand has been
     * characterized. An upper bound on R is meaningless and none is returned.
     *
     * <p>The paper's example: K = 2 devices, N = 3 customers, measured sum_k U_k = 1.6. A
     * single class admits at most 2N/(N+1) = 1.5, so Rmin = 2.
     *
     * @param Usum measured sum of device utilizations
     * @param K    number of devices
     * @param N    total number of customers
     * @return least R in 1..N with pfqn_usumbound(R,K,N) &gt;= Usum; -1 when Usum exceeds
     *         min(N,K) and so is unattainable by ANY class structure
     */
    public static int pfqn_minclasses(double Usum, int K, int N) {
        if (N < 1 || K < 1) {
            throw new RuntimeException("pfqn_minclasses requires N >= 1 and K >= 1.");
        }
        if (Usum < 0) {
            throw new RuntimeException("pfqn_minclasses requires a nonnegative utilization sum.");
        }
        double tol = 1e-12 * Math.max(1.0, Math.abs(Usum));
        for (int R = 1; R <= N; R++) {
            if (pfqn_usumbound(R, K, N) >= Usum - tol) {
                return R;
            }
        }
        return -1;
    }
}
