package jline.api.moment;

import jline.util.matrix.Matrix;

/**
 * Binomial moments from factorial moments.
 *
 * <p>Converts the factorial moments f_n = E[N(N-1)...(N-n+1)] of a discrete random variable N into the binomial moments b_n = E[binom(N,n)] via the one-to-one correspondence
 *
 * <pre>
 *   b_n = f_n / n!
 * </pre>
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003, Section 4.
 *
 * @since LINE 3.0
 */
public final class Moment_binomial_from_factorial {
    private Moment_binomial_from_factorial() {}

    /**
     * Binomial moments from factorial moments.
     *
     * @param f column vector of length n+1 holding the moments of order 0..n,
     *          element 0 being the order-0 moment, which equals 1
     * @return column vector of length n+1 holding the converted moments
     */
    public static Matrix moment_binomial_from_factorial(Matrix f) {
        int n = f.length() - 1;
        Matrix out = new Matrix(n + 1, 1);
        double fact = 1.0;
        for (int i = 0; i <= n; i++) {
            if (i > 0) {
                fact = fact * i;
            }
            out.set(i, 0, f.get(i) / fact);
        }
        return out;
    }
}
