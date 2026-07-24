package jline.api.moment;

import jline.util.matrix.Matrix;

/**
 * Power (raw) moments from cumulants.
 *
 * <p>Runs the exponential-formula recursion forward,
 *
 * <pre>
 *   m_n = sum_{k=1}^{n} nchoosek(n-1,k-1) * kappa_k * m_(n-k)
 * </pre>
 *
 * <p>with m_0 = 1. Inverse of {@link Moment_cumulant_from_raw}.
 *
 * <p>Reference:
 * V. P. Leonov and A. N. Shiryaev. On a method of calculation of
 * semi-invariants. Theory of Probability and its Applications,
 * 4(3):319-329, 1959.
 *
 * @since LINE 3.0
 */
public final class Moment_raw_from_cumulant {
    private Moment_raw_from_cumulant() {}

    /**
     * Converts cumulants into power (raw) moments.
     *
     * @param kappa column vector of length n+1 holding kappa_0,...,kappa_n;
     *              element 0 is ignored
     * @return column vector of length n+1 holding m_0,...,m_n, with m_0 = 1
     */
    public static Matrix moment_raw_from_cumulant(Matrix kappa) {
        int n = kappa.length() - 1;
        Matrix m = new Matrix(n + 1, 1);
        m.set(0, 0, 1.0);
        for (int i = 1; i <= n; i++) {
            double acc = 0.0;
            for (int k = 1; k <= i; k++) {
                acc += Moment_binotrans.nchoosek(i - 1, k - 1) * kappa.get(k) * m.get(i - k);
            }
            m.set(i, 0, acc);
        }
        return m;
    }
}
