package jline.api.moment;

/**
 * Joint factorial moments from joint binomial moments.
 *
 * <p>Converts joint binomial moments into joint factorial moments. Inverse of
 * moment_joint_binomial_from_factorial.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003.
 *
 * @since LINE 3.0
 */
public final class Moment_joint_factorial_from_binomial {
    private Moment_joint_factorial_from_binomial() {}

    /**
     * Converts joint binomial moments into joint factorial moments.
     *
     * @param b flattened joint binomial moments in row-major order, of extent
     *        dims[0]*...*dims[d-1]
     * @param dims extents of the array, dims[j] = n_j+1
     * @return flattened joint factorial moments, same layout
     */
    public static double[] moment_joint_factorial_from_binomial(double[] b, int[] dims) {
        return Moment_jointtrans.moment_jointtrans(b, dims, "factorial_from_binomial");
    }
}
