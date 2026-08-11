package jline.api.moment;

import jline.util.matrix.Matrix;

/**
 * Upward-factorial moments from factorial moments.
 *
 * <p>Converts the factorial moments f_n = E[N(N-1)...(N-n+1)] of a discrete
 * random variable N into the upward-factorial moments
 * f_n^+ = E[N(N+1)...(N+n-1)] by means of the Lah numbers,
 *
 * <pre>
 *   f_n^+ = sum_{k=1}^{n} L(n,k) * f_k   for n &gt;= 1
 *   f_0^+ = 1
 * </pre>
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003, Section 4.
 *
 * @since LINE 3.0
 */
public final class Moment_upfactorial_from_factorial {
    private Moment_upfactorial_from_factorial() {}

    /**
     * Converts factorial moments into upward-factorial moments.
     *
     * @param f column vector of length n+1 holding f_0,...,f_n, i.e. element i
     *          is the moment of order i and element 0 is f_0 = 1
     * @return column vector of length n+1 holding f_0^+,...,f_n^+
     */
    public static Matrix moment_upfactorial_from_factorial(Matrix f) {
        int n = f.length() - 1;
        Matrix L = Moment_lah.moment_lah(n);
        Matrix fp = new Matrix(n + 1, 1);
        fp.set(0, 0, 1.0);
        for (int i = 1; i <= n; i++) {
            double acc = 0.0;
            for (int k = 1; k <= i; k++) {
                acc += L.get(i, k) * f.get(k);
            }
            fp.set(i, 0, acc);
        }
        return fp;
    }
}
