package jline.api.moment;

/**
 * Joint factorial moments from joint raw moments.
 *
 * <p>Converts the joint power (raw) moments m_(i_1,...,i_d) = E[prod_j N_j^(i_j)]
 * of a random vector (N_1,...,N_d) into the joint factorial moments
 * f_(i_1,...,i_d) = E[prod_j (N_j)_(i_j)], where (N)_i = N(N-1)...(N-i+1), by
 * applying the signed Stirling numbers of the first kind separately along
 * every dimension,
 *
 *   f_(i) = sum_(k) prod_j s(i_j,k_j) * m_(k)
 *
 * The joint conversion is the Kronecker product of the univariate ones,
 * which is what makes the mode-by-mode evaluation legitimate. Only the
 * cumulant and the central conversions are not of this separable form.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003.
 *
 * @since LINE 3.0
 */
public final class Moment_joint_factorial_from_raw {
    private Moment_joint_factorial_from_raw() {}

    /**
     * Converts joint raw moments into joint factorial moments.
     *
     * @param m flattened joint raw moments in row-major order, of extent
     *        dims[0]*...*dims[d-1]
     * @param dims extents of the array, dims[j] = n_j+1
     * @return flattened joint factorial moments, same layout
     */
    public static double[] moment_joint_factorial_from_raw(double[] m, int[] dims) {
        return Moment_jointtrans.moment_jointtrans(m, dims, "factorial_from_raw");
    }
}
