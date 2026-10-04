package jline.api.moment;

/**
 * Joint binomial moments from joint survival probabilities.
 *
 * <p>For a nonnegative integer random vector with joint survival array
 * t_(m_1,...,m_d) = P(N_1 &gt;= m_1, ..., N_d &gt;= m_d),
 *
 * <pre>
 *   b_(k) = E[prod_j nchoosek(N_j,k_j)]
 *         = sum_(m&gt;=k) prod_j nchoosek(m_j-1,k_j-1) * t_(m)
 * </pre>
 *
 * <p>The transform is the tensor product of the univariate one, which is what
 * makes the mode-by-mode evaluation legitimate: an entry with k_j = 0 selects
 * m_j = 0, and t_(0,m_2,...) is by construction the marginal survival array of
 * the remaining coordinates. As in the univariate case it is upper triangular,
 * so the array must cover the joint support to be exact.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003.
 *
 * @since LINE 3.0
 */
public final class Moment_joint_binomial_from_tail {
    private Moment_joint_binomial_from_tail() {}

    /**
     * Converts joint survival probabilities into joint binomial moments.
     *
     * @param t flattened joint survival array in row-major order
     * @param dims extents of the array, dims[j] = n_j+1
     * @return flattened joint binomial moments, same layout
     */
    public static double[] moment_joint_binomial_from_tail(double[] t, int[] dims) {
        return Moment_jointtrans.moment_jointtrans(t, dims, "binomial_from_tail");
    }
}
