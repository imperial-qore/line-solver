package jline.api.moment;

import jline.util.matrix.Matrix;

/**
 * Power (raw) moments from factorial moments.
 *
 * <p>Converts the factorial moments f_n = E[N(N-1)...(N-n+1)] of a discrete
 * random variable N into the power (raw) moments m_n = E[N^n] by means of the
 * Stirling numbers of the second kind,
 *
 * <pre>
 *   m_n = sum_{k=0}^{n} S(n,k) * f_k
 * </pre>
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003, eq. (13).
 *
 * @since LINE 3.0
 */
public final class Moment_raw_from_factorial {
    private Moment_raw_from_factorial() {}

    /**
     * Converts factorial moments into power (raw) moments.
     *
     * @param f column vector of length n+1 holding f_0,...,f_n, i.e. element i
     *          is the moment of order i and element 0 is f_0 = 1
     * @return column vector of length n+1 holding m_0,...,m_n
     */
    public static Matrix moment_raw_from_factorial(Matrix f) {
        int n = f.length() - 1;
        Matrix fcol = new Matrix(n + 1, 1);
        for (int i = 0; i <= n; i++) {
            fcol.set(i, 0, f.get(i));
        }
        return Moment_stirling2.moment_stirling2(n).mult(fcol);
    }
}
