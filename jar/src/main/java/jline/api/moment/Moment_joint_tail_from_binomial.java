package jline.api.moment;

/**
 * Joint survival probabilities from joint binomial moments.
 *
 * <p>Inverse of {@link Moment_joint_binomial_from_tail}.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003.
 *
 * @since LINE 3.0
 */
public final class Moment_joint_tail_from_binomial {
    private Moment_joint_tail_from_binomial() {}

    /**
     * Converts joint binomial moments into joint survival probabilities.
     *
     * @param b flattened joint binomial moments in row-major order
     * @param dims extents of the array, dims[j] = n_j+1
     * @return flattened joint survival array, same layout
     */
    public static double[] moment_joint_tail_from_binomial(double[] b, int[] dims) {
        return Moment_jointtrans.moment_jointtrans(b, dims, "tail_from_binomial");
    }
}
