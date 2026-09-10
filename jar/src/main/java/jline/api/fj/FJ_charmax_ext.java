/**
 * @file Characteristic maximum for lattice laws and the Blom plotting position
 *
 * Ports of matlab/src/api/fj/fj_char_max_discrete.m and fj_char_max_blom.m,
 * Section 4.9 of A. Thomasian, "Analysis of Fork/Join and Related Queueing
 * Systems", ACM Computing Surveys 47(2), Article 17, 2014 (Eq. (48) and the
 * Blom and Kruskal-Weiss positions).
 *
 * @since LINE 3.0
 */
package jline.api.fj;

import java.util.function.DoubleUnaryOperator;

public final class FJ_charmax_ext {
    private FJ_charmax_ext() {}

    /** The lattice laws whose characteristic maximum is closed. */
    public enum DiscreteDist { GEOMETRIC, POISSON }

    /** [MK, mK, exact] of fj_char_max_discrete. */
    public static final class FJCharMaxDiscreteResult {
        public final double MK;
        public final int mK;
        public final double exact;

        public FJCharMaxDiscreteResult(double MK, int mK, double exact) {
            this.MK = MK;
            this.mK = mK;
            this.exact = exact;
        }
    }

    /**
     * Gravey's characteristic maximum for a lattice law. With m_K the smallest
     * integer at which P(X &gt; m_K) &lt;= 1/K,
     *
     * M_K = m_K + K sum_{k &gt;= m_K} P(X &gt; k),
     *
     * which upper bounds the expected maximum of K i.i.d. copies at O(1)
     * instead of the alternating binomial sum. For the geometric the tail sum
     * is K p^(m_K+1)/(1-p); for the Poisson it is the same sum rewritten
     * through E[(X-m)^+] = theta P(X&gt;m-1) - m P(X&gt;m).
     */
    public static FJCharMaxDiscreteResult fj_char_max_discrete(int K, DiscreteDist dist,
                                                               double par) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }
        if (dist == DiscreteDist.GEOMETRIC) {
            double p = par;
            if (!(p > 0) || !(p < 1)) {
                throw new IllegalArgumentException(
                        "The geometric parameter must lie in (0,1). Got p=" + p + ".");
            }
            int mK = (int) Math.max(0, Math.ceil(-Math.log(K) / Math.log(p) - 1e-12));
            double MK = mK + K * Math.pow(p, mK + 1) / (1 - p);
            // Exact maximum by inclusion-exclusion on the geometric tail
            double exact = 0;
            double pk = 1;
            for (int k = 1; k <= K; k++) {
                pk *= p;
                double term = binom(K, k) * pk / (1 - pk);
                exact += (k % 2 == 1) ? term : -term;
            }
            return new FJCharMaxDiscreteResult(MK, mK, exact);
        }
        double theta = par;
        if (!(theta > 0)) {
            throw new IllegalArgumentException(
                    "The Poisson mean must be positive. Got theta=" + theta + ".");
        }
        int kmax = (int) Math.ceil(theta + 12 * Math.sqrt(theta) + 40);
        double[] cdf = new double[kmax + 1];
        double[] tail = new double[kmax + 1];
        double term = Math.exp(-theta);
        double acc = 0;
        for (int k = 0; k <= kmax; k++) {
            if (k > 0) {
                term = term * theta / k;
            }
            acc += term;
            cdf[k] = Math.min(acc, 1.0);
            tail[k] = 1 - cdf[k];
        }
        int mK = -1;
        for (int k = 0; k <= kmax; k++) {
            if (tail[k] <= 1.0 / K) {
                mK = k;
                break;
            }
        }
        if (mK < 0) {
            throw new IllegalArgumentException(
                    "The Poisson lattice truncation never reached a tail of 1/K.");
        }
        double tailPrev = (mK == 0) ? 1.0 : tail[mK - 1];
        double MK = mK * (1 - K * tail[mK]) + K * theta * tailPrev;
        double exact = 0;
        for (int k = 0; k <= kmax; k++) {
            exact += 1 - Math.pow(cdf[k], K);
        }
        return new FJCharMaxDiscreteResult(MK, mK, exact);
    }

    /** [mK, lo, hi] of fj_char_max_blom; the bracket is the standard normal one. */
    public static final class FJCharMaxBlomResult {
        public final double mK;
        public final double lo;
        public final double hi;
        public final boolean bracketAvailable;

        public FJCharMaxBlomResult(double mK, double lo, double hi, boolean bracketAvailable) {
            this.mK = mK;
            this.lo = lo;
            this.hi = hi;
            this.bracketAvailable = bracketAvailable;
        }
    }

    /**
     * Blom-corrected plotting position for the characteristic maximum,
     * m_K = F^-1( (K - alpha) / (K - alpha - beta + 1) ), which for
     * alpha = beta = 0 falls back on the naive K/(K+1). For the standard normal
     * the position is bracketed for K &gt;= 5 by
     * sqrt(2 ln K - ln ln K - 3) &lt; m_K &lt; sqrt(2 ln K - ln ln K).
     */
    public static FJCharMaxBlomResult fj_char_max_blom(int K, DoubleUnaryOperator Finv,
                                                       double alpha, double beta) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }
        double den = K - alpha - beta + 1;
        if (!(den > 0)) {
            throw new IllegalArgumentException(
                    "The Blom offsets leave a non-positive denominator at K=" + K + ".");
        }
        double q = (K - alpha) / den;
        if (!(q > 0) || !(q < 1)) {
            throw new IllegalArgumentException(
                    "The Blom plotting position " + q + " fell outside (0,1).");
        }
        double mK = (Finv == null) ? normalQuantile(q) : Finv.applyAsDouble(q);
        double lo = Double.NaN, hi = Double.NaN;
        boolean have = false;
        if (K >= 5) {
            double z = 2 * Math.log(K) - Math.log(Math.log(K));
            if (z > 3) {
                lo = Math.sqrt(z - 3);
                hi = Math.sqrt(z);
                have = true;
            }
        }
        return new FJCharMaxBlomResult(mK, lo, hi, have);
    }

    public static FJCharMaxBlomResult fj_char_max_blom(int K) {
        return fj_char_max_blom(K, null, 0.4886, 0.3140);
    }

    /**
     * The standard normal quantile, MATLAB's sqrt(2) erfinv(2u-1), by bisection
     * on the complementary error function. It is called once per query here, so
     * the iterations are free and the result is tight to the last bit.
     */
    private static double normalQuantile(double u) {
        if (!(u > 0) || !(u < 1)) {
            throw new IllegalArgumentException("The probability must lie in (0,1).");
        }
        double lo = -40, hi = 40;
        for (int it = 0; it < 200; it++) {
            double mid = 0.5 * (lo + hi);
            if (mid == lo || mid == hi) {
                break;
            }
            if (0.5 * erfc(-mid / Math.sqrt(2.0)) < u) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        return 0.5 * (lo + hi);
    }

    /** Complementary error function, Numerical Recipes' Chebyshev form. */
    private static double erfc(double x) {
        double z = Math.abs(x);
        double t = 2.0 / (2.0 + z);
        double ty = 4.0 * t - 2.0;
        double[] cof = {-1.3026537197817094, 6.4196979235649026e-1, 1.9476473204185836e-2,
                -9.561514786808631e-3, -9.46595344482036e-4, 3.66839497852761e-4,
                4.2523324806907e-5, -2.0278578112534e-5, -1.624290004647e-6,
                1.303655835580e-6, 1.5626441722e-8, -8.5238095915e-8, 6.529054439e-9,
                5.059343495e-9, -9.91364156e-10, -2.27365122e-10, 9.6467911e-11,
                2.394038e-12, -6.886027e-12, 8.94487e-13, 3.13092e-13, -1.12708e-13,
                3.81e-16, 7.106e-15};
        double d = 0, dd = 0;
        for (int j = cof.length - 1; j > 0; j--) {
            double tmp = d;
            d = ty * d - dd + cof[j];
            dd = tmp;
        }
        double ans = t * Math.exp(-z * z + 0.5 * (cof[0] + ty * d) - dd);
        return x >= 0 ? ans : 2.0 - ans;
    }

    private static double binom(int n, int k) {
        if (k > n) {
            return 0;
        }
        int kk = Math.min(k, n - k);
        double r = 1;
        for (int i = 1; i <= kk; i++) {
            r = r * (n - kk + i) / i;
        }
        return r;
    }
}
