package jline.api.moment;

import jline.util.matrix.Matrix;

/**
 * Power (raw) moments from upward-factorial moments.
 *
 * <p>Converts the upward-factorial moments f_n^+ = E[N(N+1)...(N+n-1)] of a
 * discrete random variable N into the power (raw) moments m_n = E[N^n] by means
 * of the signed Stirling numbers of the second kind,
 *
 * <pre>
 *   m_n = sum_{k=0}^{n} (-1)^(n-k) * S(n,k) * f_k^+
 * </pre>
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003, Section 4.
 *
 * @since LINE 3.0
 */
public final class Moment_raw_from_upfactorial {
    private Moment_raw_from_upfactorial() {}

    /**
     * Converts upward-factorial moments into power (raw) moments.
     *
     * @param fp column vector of length n+1 holding f_0^+,...,f_n^+, i.e.
     *           element i is the moment of order i and element 0 is f_0^+ = 1
     * @return column vector of length n+1 holding m_0,...,m_n
     */
    public static Matrix moment_raw_from_upfactorial(Matrix fp) {
        int n = fp.length() - 1;
        Matrix fpcol = new Matrix(n + 1, 1);
        for (int i = 0; i <= n; i++) {
            fpcol.set(i, 0, fp.get(i));
        }
        Matrix S = Moment_stirling2.moment_stirling2(n);
        Matrix T = new Matrix(n + 1, n + 1);
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= i; j++) {
                double sign = ((i - j) % 2 == 0) ? 1.0 : -1.0;
                T.set(i, j, sign * S.get(i, j));
            }
        }
        return T.mult(fpcol);
    }
}
