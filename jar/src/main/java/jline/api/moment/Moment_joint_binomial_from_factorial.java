package jline.api.moment;

/**
 * Joint binomial moments from joint factorial moments.
 *
 * <p>Converts joint factorial moments into the joint binomial moments
 * b_(i_1,...,i_d) = E[prod_j nchoosek(N_j,i_j)] = f_(i) / prod_j (i_j!).
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003.
 *
 * @since LINE 3.0
 */
public final class Moment_joint_binomial_from_factorial {
    private Moment_joint_binomial_from_factorial() {}

    /**
     * Converts joint factorial moments into joint binomial moments.
     *
     * @param f flattened joint factorial moments in row-major order, of extent
     *        dims[0]*...*dims[d-1]
     * @param dims extents of the array, dims[j] = n_j+1
     * @return flattened joint binomial moments, same layout
     */
    public static double[] moment_joint_binomial_from_factorial(double[] f, int[] dims) {
        return Moment_jointtrans.moment_jointtrans(f, dims, "binomial_from_factorial");
    }
}
