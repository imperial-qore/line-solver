package jline.api.moment;

import jline.util.matrix.Matrix;

/**
 * Survival (tail) probabilities from binomial moments.
 *
 * <pre>
 *   t_m = sum_{j&gt;=m} (-1)^(j-m) * nchoosek(j-1,m-1) * b_j,  m &gt;= 1
 * </pre>
 *
 * <p>with t_0 = 1. Inverse of {@link Moment_binomial_from_tail}. The inversion
 * is exact on the finite box supplied, the matrix being unit upper triangular,
 * but it reconstructs the true tail only if the binomial moments were
 * themselves those of a law supported on 0,...,n.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003.
 *
 * @since LINE 3.0
 */
public final class Moment_tail_from_binomial {
    private Moment_tail_from_binomial() {}

    /**
     * Converts binomial moments into survival probabilities.
     *
     * @param b column vector of length n+1 holding b_0,...,b_n
     * @return column vector of length n+1 holding t_0,...,t_n
     */
    public static Matrix moment_tail_from_binomial(Matrix b) {
        int n = b.length() - 1;
        Matrix bcol = new Matrix(n + 1, 1);
        for (int i = 0; i <= n; i++) {
            bcol.set(i, 0, b.get(i));
        }
        return Moment_housematrix.moment_housematrix("tail_from_binomial", n).mult(bcol);
    }
}
