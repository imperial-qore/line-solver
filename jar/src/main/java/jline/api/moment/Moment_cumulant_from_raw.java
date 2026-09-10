package jline.api.moment;

import jline.util.matrix.Matrix;

/**
 * Cumulants from power (raw) moments.
 *
 * <p>Converts the power moments m_n = E[X^n] of a random variable X into its
 * cumulants kappa_n, the coefficients of the cumulant generating function
 * log E[exp(sX)] = sum_{n&gt;=1} kappa_n s^n / n!, by inverting the
 * exponential-formula recursion
 *
 * <pre>
 *   m_n = sum_{k=1}^{n} nchoosek(n-1,k-1) * kappa_k * m_(n-k)
 * </pre>
 *
 * <p>Equivalently kappa_n is the Leonov-Shiryaev partition sum over the set
 * partitions of {1,...,n}. The first cumulants are kappa_1 = m_1,
 * kappa_2 = m_2 - m_1^2 (the variance) and kappa_3 = m_3 - 3 m_1 m_2 + 2 m_1^3
 * (the third central moment). The conversion is not restricted to discrete
 * random variables.
 *
 * <p>Reference:
 * V. P. Leonov and A. N. Shiryaev. On a method of calculation of
 * semi-invariants. Theory of Probability and its Applications,
 * 4(3):319-329, 1959.
 *
 * @since LINE 3.0
 */
public final class Moment_cumulant_from_raw {
    private Moment_cumulant_from_raw() {}

    /**
     * Converts power (raw) moments into cumulants.
     *
     * @param m column vector of length n+1 holding m_0,...,m_n, with m_0 = 1
     * @return column vector of length n+1 holding kappa_0,...,kappa_n, element
     *         0 being kappa_0 = 0 and not m_0 = 1
     */
    public static Matrix moment_cumulant_from_raw(Matrix m) {
        int n = m.length() - 1;
        Matrix kappa = new Matrix(n + 1, 1);
        for (int i = 1; i <= n; i++) {
            double acc = 0.0;
            for (int k = 1; k < i; k++) {
                acc += Moment_binotrans.nchoosek(i - 1, k - 1) * kappa.get(k) * m.get(i - k);
            }
            kappa.set(i, 0, m.get(i) - acc);
        }
        return kappa;
    }
}
