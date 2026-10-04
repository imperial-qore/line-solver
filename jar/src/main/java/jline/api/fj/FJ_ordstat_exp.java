/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.fj;

import jline.util.Maths;
import jline.util.matrix.Matrix;

import org.apache.commons.math3.util.FastMath;

/**
 * Mean of the k-th smallest of n independent EXPONENTIAL branch completion times, i.e. the
 * instant a k-of-n (quorum) join fires. k = n is the ordinary AND-join, the maximum, and k = 1
 * the minimum.
 *
 * <p>With lambda_i = 1/ri(i) and m = n-k stragglers allowed,</p>
 *
 * <pre>
 *   E[X_(k)] = sum_{j=m+1..n} (-1)^(j-m-1) C(j-1,m) e_j,
 *   e_j      = sum_{|S|=j} 1 / sum_{i in S} lambda_i
 * </pre>
 *
 * <p>the inclusion-exclusion identity for the order statistics of independent exponentials. At
 * m = 0 it collapses to sum_j (-1)^(j-1) e_j, the classical expression for the maximum, TERM BY
 * TERM: a full join therefore evaluates exactly as it did before this class existed.</p>
 *
 * <p>The sum has 2^n terms and its signs alternate, so it is evaluated exactly only while the
 * branch count is small. Beyond MAXEXACT branches a genuine quorum (k &lt; n) is evaluated by
 * {@link FJ_quorum#quorumMoments(double[], double[], int)} instead, whose Poisson-binomial
 * recurrence adds no cancellation; a full join keeps the exact path at every n so that no
 * existing result moves.</p>
 *
 * <p>Port of matlab/src/api/fj/fj_ordstat_exp.m.</p>
 *
 * <p>Reference: A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems", ACM Computing
 * Surveys 47(2), Article 17, 2014, Sec. 3 (Eq. 18-19).</p>
 */
public class FJ_ordstat_exp {

    /** Branch count above which a genuine quorum leaves the exact alternating sum. */
    public static final int MAXEXACT = 15;

    /**
     * @param ri branch completion time means
     * @param k  quorum, 1 &lt;= k &lt;= ri.length
     * @return   the mean instant the k-of-n join fires
     */
    public static double fj_ordstat_exp(double[] ri, int k) {
        int nAll = 0;
        for (int i = 0; i < ri.length; i++) {
            if (!Double.isNaN(ri[i]) && !Double.isInfinite(ri[i])) {
                nAll++;
            }
        }
        if (nAll == 0) {
            return 0.0;
        }
        double[] r = new double[nAll];
        int p = 0;
        for (int i = 0; i < ri.length; i++) {
            if (!Double.isNaN(ri[i]) && !Double.isInfinite(ri[i])) {
                r[p++] = ri[i];
            }
        }
        int n = nAll;
        if (k < 1 || k > n) {
            throw new IllegalArgumentException(
                    "fj_ordstat_exp: k must satisfy 1 <= k <= n. Got k=" + k + ", n=" + n + ".");
        }
        // A branch of zero mean completes instantly: it never delays the join and it counts
        // toward the quorum at once. Removing it here keeps the reciprocal below finite, which
        // an exact arithmetic requires and IEEE only tolerates.
        int nzero = 0;
        for (int i = 0; i < n; i++) {
            if (r[i] <= 0) nzero++;
        }
        if (nzero > 0) {
            k -= nzero;
            if (k <= 0) {
                return 0.0;
            }
            double[] pos = new double[n - nzero];
            int q = 0;
            for (int i = 0; i < n; i++) {
                if (r[i] > 0) pos[q++] = r[i];
            }
            r = pos;
            n = pos.length;
        }
        if (n == 1) {
            return r[0];
        }

        if (k < n && n > MAXEXACT) {
            // Branch times are taken as exponential, so the variance is the square of the mean.
            double[] var = new double[n];
            for (int i = 0; i < n; i++) {
                var[i] = r[i] * r[i];
            }
            return FJ_quorum.quorumMoments(r, var, k)[0];
        }

        Matrix lambdai = new Matrix(1, n);
        for (int i = 0; i < n; i++) {
            lambdai.set(0, i, 1.0 / r[i]);
        }
        int nstrag = n - k;
        double m = 0;
        for (int j = nstrag + 1; j <= n; j++) {
            Matrix nk = Maths.nCk(lambdai, j).sumRows();
            double ej = Matrix.ones(nk.getNumRows(), 1).elementDiv(nk).elementSum();
            m += FastMath.pow(-1, j - nstrag - 1) * Maths.binomialCoeff(j - 1, nstrag) * ej;
        }
        return m;
    }

    /**
     * @param ri branch completion time means, as a row or column matrix
     * @param k  quorum, 1 &lt;= k &lt;= ri.length()
     * @return   the mean instant the k-of-n join fires
     */
    public static double fj_ordstat_exp(Matrix ri, int k) {
        double[] r = new double[ri.length()];
        for (int i = 0; i < r.length; i++) {
            r[i] = ri.get(i);
        }
        return fj_ordstat_exp(r, k);
    }
}
