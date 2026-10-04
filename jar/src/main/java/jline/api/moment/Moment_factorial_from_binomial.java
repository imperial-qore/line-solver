package jline.api.moment;

import jline.util.matrix.Matrix;

/**
 * Factorial moments from binomial moments.
 *
 * <p>Converts the binomial moments b_n = E[binom(N,n)] of a discrete random variable N into the factorial moments f_n = E[N(N-1)...(N-n+1)] via the one-to-one correspondence
 *
 * <pre>
 *   f_n = n! * b_n
 * </pre>
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003, Section 4.
 *
 * @since LINE 3.0
 */
public final class Moment_factorial_from_binomial {
    private Moment_factorial_from_binomial() {}

    /**
     * Factorial moments from binomial moments.
     *
     * @param b column vector of length n+1 holding the moments of order 0..n,
     *          element 0 being the order-0 moment, which equals 1
     * @return column vector of length n+1 holding the converted moments
     */
    public static Matrix moment_factorial_from_binomial(Matrix b) {
        int n = b.length() - 1;
        Matrix out = new Matrix(n + 1, 1);
        double fact = 1.0;
        for (int i = 0; i <= n; i++) {
            if (i > 0) {
                fact = fact * i;
            }
            out.set(i, 0, b.get(i) * fact);
        }
        return out;
    }
}
