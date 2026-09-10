package jline.api.moment;

import jline.util.Maths;
import jline.util.matrix.Matrix;

/**
 * Binomial moments from negative-binomial moments.
 *
 * <p>Converts the negative-binomial moments b_n^- = E[binom(N+n-1,n)] of a
 * discrete random variable N into the binomial moments b_n = E[binom(N,n)] by
 * means of the shifted binomial transform
 *
 * <pre>
 *   b_n = sum_{k=1}^{n} (-1)^(n-k) * binom(n-1,k-1) * b_k^-   for n &gt;= 1
 *   b_0 = 1
 * </pre>
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003, eq. (14).
 *
 * @since LINE 3.0
 */
public final class Moment_binomial_from_negbinomial {
    private Moment_binomial_from_negbinomial() {}

    /**
     * Converts negative-binomial moments into binomial moments.
     *
     * @param bm column vector of length n+1 holding b_0^-,...,b_n^-, i.e.
     *           element i is the moment of order i and element 0 is b_0^- = 1
     * @return column vector of length n+1 holding b_0,...,b_n
     */
    public static Matrix moment_binomial_from_negbinomial(Matrix bm) {
        int n = bm.length() - 1;
        Matrix b = new Matrix(n + 1, 1);
        b.set(0, 0, 1.0);
        for (int i = 1; i <= n; i++) {
            double acc = 0.0;
            for (int k = 1; k <= i; k++) {
                acc += (((i - k) % 2 == 0) ? 1.0 : -1.0) * Maths.nchoosek(i - 1, k - 1) * bm.get(k);
            }
            b.set(i, 0, acc);
        }
        return b;
    }
}
