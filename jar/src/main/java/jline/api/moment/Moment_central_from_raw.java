package jline.api.moment;

import jline.util.Maths;
import jline.util.matrix.Matrix;

/**
 * Central moments from power (raw) moments.
 *
 * <p>Converts the power (raw) moments m_n = E[N^n] of a random variable N into
 * the central moments m_n^c = E[(N-m_1)^n] by means of the binomial transform
 * in the variation that involves the mean m_1,
 *
 * <pre>
 *   m_n^c = sum_{k=0}^{n} (-1)^(n-k) * binom(n,k) * m_k * m_1^(n-k)
 * </pre>
 *
 * <p>The conversion also holds for continuous random variables.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003, Section 4.
 *
 * @since LINE 3.0
 */
public final class Moment_central_from_raw {
    private Moment_central_from_raw() {}

    /**
     * Converts power (raw) moments into central moments.
     *
     * @param m column vector of length n+1 holding m_0,...,m_n, i.e. element i
     *          is the moment of order i and element 0 is m_0 = 1. At least the
     *          mean m_1 must be given, hence the length must be at least 2
     * @return column vector of length n+1 holding m_0^c,...,m_n^c, with
     *         m_0^c = 1 and m_1^c = 0 by construction
     * @throws IllegalArgumentException if m holds fewer than 2 elements
     */
    public static Matrix moment_central_from_raw(Matrix m) {
        int n = m.length() - 1;
        if (n < 1) {
            throw new IllegalArgumentException(
                    "The mean m_1 is required for this conversion, hence m must have at least 2 elements.");
        }
        double m1 = m.get(1);
        Matrix mc = new Matrix(n + 1, 1);
        for (int i = 0; i <= n; i++) {
            double acc = 0.0;
            for (int k = 0; k <= i; k++) {
                acc += (((i - k) % 2 == 0) ? 1.0 : -1.0) * Maths.nchoosek(i, k) * m.get(k)
                        * Math.pow(m1, i - k);
            }
            mc.set(i, 0, acc);
        }
        return mc;
    }
}
