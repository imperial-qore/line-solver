package jline.api.moment;

import jline.util.matrix.Matrix;

/**
 * Binomial moments from survival (tail) probabilities.
 *
 * <p>For a nonnegative integer random variable N with survival sequence
 * t_m = P(N &gt;= m),
 *
 * <pre>
 *   b_j = E[nchoosek(N,j)] = sum_{m&gt;=j} nchoosek(m-1,j-1) * t_m,  j &gt;= 1
 * </pre>
 *
 * <p>with b_0 = t_0 = 1. Unlike every other edge of the house of moments this
 * transform is UPPER triangular, so it consumes the whole tail: the result is
 * exact only if the sequence covers the support, that is t_m = 0 beyond the
 * last element supplied. This is the natural entry point for a closed queueing
 * network, whose queue lengths are bounded by the population and whose joint
 * survival probabilities are ratios of normalizing constants.
 *
 * <p>Truncating the tail early yields a strict LOWER bound on every b_j, since
 * all the coefficients and all the tail values are nonnegative. The bound is
 * not inherited by the central moments downstream, whose conversion alternates
 * in sign.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003.
 *
 * @since LINE 3.0
 */
public final class Moment_binomial_from_tail {
    private Moment_binomial_from_tail() {}

    /**
     * Converts survival probabilities into binomial moments.
     *
     * @param t column vector of length n+1 holding t_0,...,t_n, element m being
     *          P(N &gt;= m) and element 0 being 1
     * @return column vector of length n+1 holding b_0,...,b_n
     */
    public static Matrix moment_binomial_from_tail(Matrix t) {
        int n = t.length() - 1;
        Matrix tcol = new Matrix(n + 1, 1);
        for (int i = 0; i <= n; i++) {
            tcol.set(i, 0, t.get(i));
        }
        return Moment_housematrix.moment_housematrix("binomial_from_tail", n).mult(tcol);
    }
}
