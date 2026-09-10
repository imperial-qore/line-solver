package jline.api.moment;

/**
 * Joint cumulants from joint power (raw) moments.
 *
 * <p>The joint cumulants of a random vector (N_1,...,N_d) are the coefficients
 * of the joint cumulant generating function
 *
 * <pre>
 *   log E[exp(s_1 N_1 + ... + s_d N_d)] = sum_(a != 0) kappa_a prod_j
 *   s_j^(a_j) / a_j!
 * </pre>
 *
 * <p>They obey the multivariate exponential formula, equivalently the
 * Leonov-Shiryaev partition formula. With j the first dimension in which the
 * multi-index a is nonzero,
 *
 * <pre>
 *   m_a = sum_(0&lt;b&lt;=a) prod_l nchoosek(a_l-[l=j], b_l-[l=j]) kappa_b m_(a-b)
 * </pre>
 *
 * <p>which isolates kappa_a because the b = a term has unit coefficient and
 * m_0 = 1. Unlike every other conversion in the house, this one does not factor
 * into a product of univariate transforms: the cumulant of multi-order (1,1) is
 * the covariance, which mixes the dimensions. Multi-indices are swept in
 * row-major (lexicographic) order, under which every b &lt;= a precedes a.
 *
 * <p>Reference:
 * V. P. Leonov and A. N. Shiryaev. On a method of calculation of
 * semi-invariants. Theory of Probability and its Applications,
 * 4(3):319-329, 1959.
 *
 * @since LINE 3.0
 */
public final class Moment_joint_cumulant_from_raw {
    private Moment_joint_cumulant_from_raw() {}

    /**
     * Converts joint power moments into joint cumulants.
     *
     * @param m flattened joint power moments in row-major order, element 0
     *          being 1
     * @param dims extents of the array, dims[j] = n_j+1
     * @return flattened joint cumulants, same layout, element 0 being 0
     */
    public static double[] moment_joint_cumulant_from_raw(double[] m, int[] dims) {
        return Moment_joint_raw_from_cumulant.recur(m, dims, true);
    }
}
