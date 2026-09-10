/**
 * @file Maxima of heterogeneous and general branch distributions
 *
 * Ports of matlab/src/api/fj/fj_xmax_het.m, fj_lst_max_het.m,
 * fj_xmax_moments_het.m, fj_xmax_hz.m, fj_xmax_hz_het.m, fj_cox_fit.m and
 * fj_xmax_coxian.m, the order-statistics group of A. Thomasian, "Analysis of
 * Fork/Join and Related Queueing Systems", ACM Computing Surveys 47(2),
 * Article 17, 2014 (Eqs. (27)-(30), (38)-(40), (46)).
 *
 * @since LINE 3.0
 */
package jline.api.fj;

import java.util.function.DoubleUnaryOperator;

public final class FJ_maxima {
    private FJ_maxima() {}

    /**
     * Exact n-th moment of the maximum of K heterogeneous exponentials, by
     * inclusion-exclusion on the survival function:
     *
     * E[Y^n] = sum over the nonempty subsets S of (-1)^(|S|+1) n! / (sum_S lambda_i)^n.
     *
     * At n = 1 and K = 2 this collapses to 1/l1 + 1/l2 - 1/(l1+l2), and for
     * equal rates to H_K/lambda.
     */
    public static double fj_xmax_het(double[] lambda, int n) {
        int K = lambda.length;
        if (K < 1) {
            throw new IllegalArgumentException("At least one rate is required.");
        }
        if (n < 1) {
            throw new IllegalArgumentException("The moment order must be a positive integer.");
        }
        if (K > 24) {
            throw new IllegalArgumentException(
                    "Inclusion-exclusion needs 2^K terms; use fj_xmax_moments_het instead.");
        }
        for (int i = 0; i < K; i++) {
            if (!(lambda[i] > 0)) {
                throw new IllegalArgumentException("All exponential rates must be positive.");
            }
        }
        double nfact = 1;
        for (int i = 2; i <= n; i++) {
            nfact *= i;
        }
        double acc = 0;
        int nmask = 1 << K;
        for (int mask = 1; mask < nmask; mask++) {
            double rate = 0;
            int card = 0;
            for (int i = 0; i < K; i++) {
                if ((mask & (1 << i)) != 0) {
                    rate += lambda[i];
                    card++;
                }
            }
            double term = nfact / Math.pow(rate, n);
            acc += (card % 2 == 1) ? term : -term;
        }
        return acc;
    }

    public static double fj_xmax_het(double[] lambda) {
        return fj_xmax_het(lambda, 1);
    }

    /**
     * Laplace-Stieltjes transform of the maximum of heterogeneous exponentials,
     * by the Harrison-Zertal recurrence
     *
     * ( s + sum_{j=1..m} lambda_j ) L*_m = sum_{j=1..m} lambda_j L*_{m-1}(lambda \ j),
     *
     * anchored at L*_0 = 1 because the maximum of an empty collection is zero.
     */
    public static double fj_lst_max_het(double[] lambda, double s) {
        int K = lambda.length;
        if (K < 1) {
            throw new IllegalArgumentException("At least one rate is required.");
        }
        if (K > 22) {
            throw new IllegalArgumentException("The recurrence enumerates 2^K sub-collections.");
        }
        if (s < 0) {
            throw new IllegalArgumentException("The transform argument must be non-negative.");
        }
        for (int i = 0; i < K; i++) {
            if (!(lambda[i] > 0)) {
                throw new IllegalArgumentException("All exponential rates must be positive.");
            }
        }
        int nmask = 1 << K;
        double[] tab = new double[nmask];
        tab[0] = 1;
        for (int mask = 1; mask < nmask; mask++) {
            double tot = 0, acc = 0;
            for (int j = 0; j < K; j++) {
                if ((mask & (1 << j)) != 0) {
                    tot += lambda[j];
                    acc += lambda[j] * tab[mask ^ (1 << j)];
                }
            }
            tab[mask] = acc / (s + tot);
        }
        return tab[nmask - 1];
    }

    /**
     * Moments of orders 1..n of the maximum of heterogeneous exponentials, by
     * the n-th derivative of the transform recurrence at the origin:
     *
     * M_m(n) = [ n M_m(n-1) + sum_j lambda_j M_{m-1}(lambda \ j, n) ] / sum_j lambda_j.
     *
     * Eq. (30) of the survey prints the second sum WITHOUT the lambda_j weight;
     * that form is not the derivative of Eq. (29) and misses the textbook
     * two-variable answer, so the weight is restored here. fj_xmax_het is the
     * independent inclusion-exclusion check.
     */
    public static double[] fj_xmax_moments_het(double[] lambda, int n) {
        int K = lambda.length;
        if (K < 1) {
            throw new IllegalArgumentException("At least one rate is required.");
        }
        if (n < 1) {
            throw new IllegalArgumentException("The moment order must be a positive integer.");
        }
        if (K > 22) {
            throw new IllegalArgumentException("The recurrence enumerates 2^K sub-collections.");
        }
        for (int i = 0; i < K; i++) {
            if (!(lambda[i] > 0)) {
                throw new IllegalArgumentException("All exponential rates must be positive.");
            }
        }
        int nmask = 1 << K;
        double[][] tab = new double[nmask][n + 1];
        for (int mask = 0; mask < nmask; mask++) {
            tab[mask][0] = 1;
        }
        for (int order = 1; order <= n; order++) {
            // The empty sub-collection has a zero maximum, so all its moments vanish
            tab[0][order] = 0;
            for (int mask = 1; mask < nmask; mask++) {
                double tot = 0;
                double acc = order * tab[mask][order - 1];
                for (int j = 0; j < K; j++) {
                    if ((mask & (1 << j)) != 0) {
                        tot += lambda[j];
                        acc += lambda[j] * tab[mask ^ (1 << j)][order];
                    }
                }
                tab[mask][order] = acc / tot;
            }
        }
        double[] out = new double[n];
        for (int k = 1; k <= n; k++) {
            out[k - 1] = tab[nmask - 1][k];
        }
        return out;
    }

    /**
     * Harrison-Zertal closed form for K identically distributed branches given
     * their first two moments, X_K^max ~ m1 + ( m2/(2 m1) ) ( H_K - 1 ).
     *
     * The correction is the equilibrium mean of the branch law scaled by
     * H_K - 1, which is exact for the exponential and reduces to m1 at K = 1
     * for every branch law.
     */
    public static double fj_xmax_hz(double m1, double m2, int K) {
        if (!(m1 > 0)) {
            throw new IllegalArgumentException("The branch mean must be positive.");
        }
        if (m2 < m1 * m1) {
            throw new IllegalArgumentException(
                    "The second moment is below the square of the mean.");
        }
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }
        return m1 + (m2 / (2 * m1)) * (FJ_harmonic.fj_harmonic(K) - 1);
    }

    /**
     * Harrison-Zertal recurrence for independent but not identically
     * distributed branches, each supplied through its first two moments and its
     * distribution function:
     *
     * I(S) = (1/|S|) sum_{i in S} [ I(S \ i) + (m2_i/(2 m1_i)) L*_{S\i}(1/m1_i) ],
     *
     * anchored at I({i}) = m1_i. The transform of the maximum over a
     * sub-collection is recovered from the product of the distribution
     * functions by composite Simpson quadrature.
     */
    public static double fj_xmax_hz_het(double[] m1, double[] m2, DoubleUnaryOperator[] cdf,
                                        double tol, int npanels) {
        int K = m1.length;
        if (m2.length != K || cdf.length != K) {
            throw new IllegalArgumentException("m1, m2 and cdf must have the same length.");
        }
        if (K < 1) {
            throw new IllegalArgumentException("At least one branch is required.");
        }
        if (K > 14) {
            throw new IllegalArgumentException(
                    "The recurrence enumerates 2^K sub-collections with a quadrature each.");
        }
        for (int i = 0; i < K; i++) {
            if (!(m1[i] > 0)) {
                throw new IllegalArgumentException("All branch means must be positive.");
            }
            if (m2[i] < m1[i] * m1[i]) {
                throw new IllegalArgumentException(
                        "A second moment is below the square of its mean.");
            }
        }
        if (npanels % 2 != 0) {
            npanels++;
        }
        double[] alpha = new double[K];
        double[] resid = new double[K];
        double mmax = m1[0];
        for (int i = 0; i < K; i++) {
            alpha[i] = 1.0 / m1[i];
            resid[i] = m2[i] / (2 * m1[i]);
            if (m1[i] > mmax) {
                mmax = m1[i];
            }
        }
        // Horizon: widen until every branch is essentially complete
        double U = 8 * mmax;
        for (int it = 0; it < 60; it++) {
            double prodF = 1;
            for (int j = 0; j < K; j++) {
                prodF *= cdf[j].applyAsDouble(U);
            }
            if (1 - prodF < tol) {
                break;
            }
            U *= 2;
        }
        int nmask = 1 << K;
        double[] Ival = new double[nmask];
        for (int mask = 1; mask < nmask; mask++) {
            int card = Integer.bitCount(mask);
            if (card == 1) {
                Ival[mask] = m1[Integer.numberOfTrailingZeros(mask)];
                continue;
            }
            double acc = 0;
            for (int i = 0; i < K; i++) {
                if ((mask & (1 << i)) == 0) {
                    continue;
                }
                int rest = mask ^ (1 << i);
                acc += Ival[rest] + resid[i] * lstMax(cdf, rest, K, alpha[i], U, npanels);
            }
            Ival[mask] = acc / card;
        }
        return Ival[nmask - 1];
    }

    public static double fj_xmax_hz_het(double[] m1, double[] m2, DoubleUnaryOperator[] cdf) {
        return fj_xmax_hz_het(m1, m2, cdf, 1e-10, 2000);
    }

    /** L*_T(s) for the sub-collection selected by mask, by composite Simpson. */
    private static double lstMax(DoubleUnaryOperator[] cdf, int mask, int K, double s, double U,
                                 int npanels) {
        if (mask == 0) {
            return 1;
        }
        double h = U / npanels;
        double acc = 0;
        for (int i = 0; i <= npanels; i++) {
            double t = h * i;
            double g = Math.exp(-s * t);
            for (int j = 0; j < K; j++) {
                if ((mask & (1 << j)) != 0) {
                    g *= cdf[j].applyAsDouble(t);
                }
            }
            double w = (i == 0 || i == npanels) ? 1 : ((i % 2 == 1) ? 4 : 2);
            acc += w * g;
        }
        return s * (h / 3) * acc;
    }

    /** [mu1, mu2, q, kmin, kmax] of fj_cox_fit. */
    public static final class FJCoxFitResult {
        public final double mu1;
        public final double mu2;
        public final double q;
        public final int kmin;
        public final int kmax;

        public FJCoxFitResult(double mu1, double mu2, double q, int kmin, int kmax) {
            this.mu1 = mu1;
            this.mu2 = mu2;
            this.q = q;
            this.kmin = kmin;
            this.kmax = kmax;
        }
    }

    /**
     * Marie's balanced-stage fit of a two-stage Coxian law to a target mean and
     * squared coefficient of variation: requiring 1/mu1 = q/mu2 closes the
     * system of two moment equations in three unknowns and gives
     * mu1 = 2 mu, q = 1/(2 c2), mu2 = mu/c2, which needs c2 &gt;= 0.5.
     */
    public static FJCoxFitResult fj_cox_fit(double m1, double c2) {
        if (!(m1 > 0)) {
            throw new IllegalArgumentException("The target mean must be positive.");
        }
        if (c2 < 0.5) {
            throw new IllegalArgumentException(
                    "The balanced-stage Coxian fit needs c2 >= 0.5. Got c2=" + c2 + ".");
        }
        double mu = 1.0 / m1;
        double inv = 1.0 / c2;
        int kmin = (int) Math.max(1, Math.ceil(inv - 1e-12));
        int kmax = (int) Math.floor(inv + 1e-12) + 1;
        if (kmax < kmin) {
            kmax = kmin;
        }
        return new FJCoxFitResult(2 * mu, mu / c2, 1.0 / (2 * c2), kmin, kmax);
    }

    /** [Xmax, m1, c2] of fj_xmax_coxian. */
    public static final class FJXmaxCoxianResult {
        public final double Xmax;
        public final double m1;
        public final double c2;

        public FJXmaxCoxianResult(double Xmax, double m1, double c2) {
            this.Xmax = Xmax;
            this.m1 = m1;
            this.c2 = c2;
        }
    }

    /**
     * Exact expected maximum of K i.i.d. two-stage Coxian branches. The survival
     * function is a two-term exponential mixture, so expanding 1 - (1-S)^K
     * binomially and integrating term by term is closed:
     *
     * E[Y_K] = sum_j (-1)^(j+1) C(K,j) sum_i C(j,i) A^(j-i) B^i / ((j-i) mu1 + i mu2).
     *
     * At coincident stage rates the mixture degenerates into
     * S(t) = (1 + q mu t) exp(-mu t), which is handled by the same expansion
     * with the polynomial integrals.
     */
    public static FJXmaxCoxianResult fj_xmax_coxian(int K, double mu1, double mu2, double q) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }
        if (!(mu1 > 0) || !(mu2 > 0)) {
            throw new IllegalArgumentException("Both stage rates must be positive.");
        }
        if (q < 0 || q > 1) {
            throw new IllegalArgumentException("The branching probability must lie in [0,1].");
        }
        if (K > 60) {
            throw new IllegalArgumentException(
                    "The binomial expansion loses precision past K=60. Got K=" + K + ".");
        }
        double m1 = 1.0 / mu1 + q / mu2;
        double var1 = 1.0 / (mu1 * mu1) + q * (2 - q) / (mu2 * mu2);
        double c2 = var1 / (m1 * m1);
        double acc = 0;
        if (Math.abs(mu2 - mu1) > 1e-12 * Math.max(mu1, mu2)) {
            double A = (1 - q) + q * mu2 / (mu2 - mu1);
            double B = -q * mu1 / (mu2 - mu1);
            for (int j = 1; j <= K; j++) {
                double inner = 0;
                for (int i = 0; i <= j; i++) {
                    double rate = (j - i) * mu1 + i * mu2;
                    inner += binom(j, i) * Math.pow(A, j - i) * Math.pow(B, i) / rate;
                }
                double term = binom(K, j) * inner;
                acc += (j % 2 == 1) ? term : -term;
            }
        } else {
            double mu = mu1;
            for (int j = 1; j <= K; j++) {
                double inner = 0;
                double jmu = j * mu;
                for (int i = 0; i <= j; i++) {
                    double fact = 1;
                    for (int e = 2; e <= i; e++) {
                        fact *= e;
                    }
                    inner += binom(j, i) * Math.pow(q * mu, i) * fact / Math.pow(jmu, i + 1);
                }
                double term = binom(K, j) * inner;
                acc += (j % 2 == 1) ? term : -term;
            }
        }
        return new FJXmaxCoxianResult(acc, m1, c2);
    }

    /** Binomial coefficient, built multiplicatively so every partial product is integral. */
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
