package jline.api.moment;

import jline.util.Maths;
import jline.util.matrix.Matrix;

/**
 * Negative-binomial moments from binomial moments.
 *
 * <p>Converts the binomial moments b_n = E[binom(N,n)] of a discrete random
 * variable N into the negative-binomial moments b_n^- = E[binom(N+n-1,n)] by
 * means of the shifted binomial transform
 *
 * <pre>
 *   b_n^- = sum_{k=1}^{n} binom(n-1,k-1) * b_k   for n &gt;= 1
 *   b_0^- = 1
 * </pre>
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003, eq. (14).
 *
 * @since LINE 3.0
 */
public final class Moment_negbinomial_from_binomial {
    private Moment_negbinomial_from_binomial() {}

    /**
     * Converts binomial moments into negative-binomial moments.
     *
     * @param b column vector of length n+1 holding b_0,...,b_n, i.e. element i
     *          is the moment of order i and element 0 is b_0 = 1
     * @return column vector of length n+1 holding b_0^-,...,b_n^-
     */
    public static Matrix moment_negbinomial_from_binomial(Matrix b) {
        int n = b.length() - 1;
        Matrix bm = new Matrix(n + 1, 1);
        bm.set(0, 0, 1.0);
        for (int i = 1; i <= n; i++) {
            double acc = 0.0;
            for (int k = 1; k <= i; k++) {
                acc += Maths.nchoosek(i - 1, k - 1) * b.get(k);
            }
            bm.set(i, 0, acc);
        }
        return bm;
    }
}
