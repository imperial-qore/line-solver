package jline.api.moment;

import jline.util.Maths;
import jline.util.matrix.Matrix;

/**
 * Power (raw) moments from central moments.
 *
 * <p>Converts the central moments m_n^c = E[(N-m_1)^n] of a random variable N
 * into the power (raw) moments m_n = E[N^n] by means of the inverse binomial
 * transform in the variation that involves the mean m_1,
 *
 * <pre>
 *   m_n = sum_{k=0}^{n} binom(n,k) * m_k^c * m_1^(n-k)
 * </pre>
 *
 * <p>The mean must be supplied separately since m_1^c = 0 carries no
 * information on it. The conversion also holds for continuous random variables.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003, Section 4.
 *
 * @since LINE 3.0
 */
public final class Moment_raw_from_central {
    private Moment_raw_from_central() {}

    /**
     * Converts central moments into power (raw) moments.
     *
     * @param mc column vector of length n+1 holding m_0^c,...,m_n^c, i.e.
     *           element i is the moment of order i and element 0 is m_0^c = 1
     * @param m1 mean of N
     * @return column vector of length n+1 holding m_0,...,m_n
     */
    public static Matrix moment_raw_from_central(Matrix mc, double m1) {
        int n = mc.length() - 1;
        Matrix m = new Matrix(n + 1, 1);
        for (int i = 0; i <= n; i++) {
            double acc = 0.0;
            for (int k = 0; k <= i; k++) {
                acc += Maths.nchoosek(i, k) * mc.get(k) * Math.pow(m1, i - k);
            }
            m.set(i, 0, acc);
        }
        return m;
    }
}
