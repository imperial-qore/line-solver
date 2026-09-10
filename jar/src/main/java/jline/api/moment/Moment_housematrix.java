package jline.api.moment;

import jline.util.matrix.Matrix;

/**
 * Conversion matrix of one edge of the house of moments.
 *
 * <p>The edge is returned as a linear map on the moment subspace {m_0 = 1}.
 * Four edges (the Lah pair and the shifted-binomial pair) pin their zeroth
 * output to 1 rather than propagating element 0, so as maps of the whole space
 * they are affine. Here the offset is folded into column 0, which is empty for
 * those edges, making every edge a genuine matrix. On a moment vector, whose
 * element 0 is 1 by definition, the two agree. This is also what makes those
 * edges usable dimension by dimension in the joint conversions.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003.
 *
 * @since LINE 3.0
 */
public final class Moment_housematrix {
    private Moment_housematrix() {}

    /**
     * Returns the conversion matrix of one edge of the house of moments.
     *
     * @param edge one of "factorial_from_raw", "raw_from_factorial",
     *             "upfactorial_from_raw", "raw_from_upfactorial",
     *             "binomial_from_factorial", "factorial_from_binomial",
     *             "negbinomial_from_upfactorial", "upfactorial_from_negbinomial",
     *             "factorial_from_upfactorial", "upfactorial_from_factorial",
     *             "negbinomial_from_binomial", "binomial_from_negbinomial",
     *             "binomial_from_tail", "tail_from_binomial"
     * @param n maximum order of the mode (n &gt;= 0)
     * @return the (n+1)x(n+1) conversion matrix
     */
    public static Matrix moment_housematrix(String edge, int n) {
        if (n < 0) {
            throw new IllegalArgumentException("moment_housematrix: The maximum order n must be "
                    + "a nonnegative integer.");
        }
        if ("factorial_from_raw".equals(edge)) {
            return Moment_stirling1.moment_stirling1(n);
        }
        if ("raw_from_factorial".equals(edge)) {
            return Moment_stirling2.moment_stirling2(n);
        }
        if ("upfactorial_from_raw".equals(edge)) {
            return Moment_stirlingcycle.moment_stirlingcycle(n);
        }
        Matrix T = new Matrix(n + 1, n + 1);
        if ("raw_from_upfactorial".equals(edge)) {
            Matrix S = Moment_stirling2.moment_stirling2(n);
            for (int i = 0; i <= n; i++) {
                for (int j = 0; j <= i; j++) {
                    double sign = ((i - j) % 2 == 0) ? 1.0 : -1.0;
                    T.set(i, j, sign * S.get(i, j));
                }
            }
            return T;
        }
        if ("binomial_from_factorial".equals(edge) || "negbinomial_from_upfactorial".equals(edge)) {
            double fact = 1.0;
            for (int i = 0; i <= n; i++) {
                if (i > 0) {
                    fact *= i;
                }
                T.set(i, i, 1.0 / fact);
            }
            return T;
        }
        if ("factorial_from_binomial".equals(edge) || "upfactorial_from_negbinomial".equals(edge)) {
            double fact = 1.0;
            for (int i = 0; i <= n; i++) {
                if (i > 0) {
                    fact *= i;
                }
                T.set(i, i, fact);
            }
            return T;
        }
        if ("factorial_from_upfactorial".equals(edge) || "upfactorial_from_factorial".equals(edge)) {
            boolean signed = "factorial_from_upfactorial".equals(edge);
            Matrix L = Moment_lah.moment_lah(n);
            T.set(0, 0, 1.0);
            for (int i = 1; i <= n; i++) {
                for (int k = 1; k <= i; k++) {
                    double sign = (!signed || (i - k) % 2 == 0) ? 1.0 : -1.0;
                    T.set(i, k, sign * L.get(i, k));
                }
            }
            return T;
        }
        if ("negbinomial_from_binomial".equals(edge) || "binomial_from_negbinomial".equals(edge)) {
            boolean signed = "binomial_from_negbinomial".equals(edge);
            T.set(0, 0, 1.0);
            for (int i = 1; i <= n; i++) {
                for (int k = 1; k <= i; k++) {
                    double sign = (!signed || (i - k) % 2 == 0) ? 1.0 : -1.0;
                    T.set(i, k, sign * Moment_binotrans.nchoosek(i - 1, k - 1));
                }
            }
            return T;
        }
        if ("binomial_from_tail".equals(edge) || "tail_from_binomial".equals(edge)) {
            boolean signed = "tail_from_binomial".equals(edge);
            T.set(0, 0, 1.0);
            for (int i = 1; i <= n; i++) {
                for (int k = i; k <= n; k++) {
                    double sign = (!signed || (k - i) % 2 == 0) ? 1.0 : -1.0;
                    T.set(i, k, sign * Moment_binotrans.nchoosek(k - 1, i - 1));
                }
            }
            return T;
        }
        throw new IllegalArgumentException("moment_housematrix: Unknown edge " + edge + ".");
    }
}
