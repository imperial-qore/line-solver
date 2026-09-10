package jline.api.moment;

import jline.util.matrix.Matrix;

/**
 * Upward-factorial moments from power (raw) moments.
 *
 * <p>Converts the power (raw) moments m_n = E[N^n] of a discrete random variable
 * N into the upward-factorial moments f_n^+ = E[N(N+1)...(N+n-1)] by means of
 * the Stirling cycle numbers,
 *
 * <pre>
 *   f_n^+ = sum_{k=0}^{n} sigma(n,k) * m_k
 * </pre>
 *
 * <p>Upward-factorial moments are of use in moment-matching techniques for
 * matrix-geometric and discrete phase-type distributions.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003, Section 4.
 *
 * @since LINE 3.0
 */
public final class Moment_upfactorial_from_raw {
    private Moment_upfactorial_from_raw() {}

    /**
     * Converts power (raw) moments into upward-factorial moments.
     *
     * @param m column vector of length n+1 holding m_0,...,m_n, i.e. element i
     *          is the moment of order i and element 0 is m_0 = 1
     * @return column vector of length n+1 holding f_0^+,...,f_n^+
     */
    public static Matrix moment_upfactorial_from_raw(Matrix m) {
        int n = m.length() - 1;
        Matrix mcol = new Matrix(n + 1, 1);
        for (int i = 0; i <= n; i++) {
            mcol.set(i, 0, m.get(i));
        }
        return Moment_stirlingcycle.moment_stirlingcycle(n).mult(mcol);
    }
}
