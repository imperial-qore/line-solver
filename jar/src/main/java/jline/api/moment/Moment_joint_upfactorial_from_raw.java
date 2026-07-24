package jline.api.moment;

/**
 * Joint upward-factorial moments from joint raw moments.
 *
 * <p>Converts joint power (raw) moments into the joint upward-factorial moments
 * f+_(i_1,...,i_d) = E[prod_j N_j(N_j+1)...(N_j+i_j-1)], by the Stirling
 * cycle numbers along every dimension.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003.
 *
 * @since LINE 3.0
 */
public final class Moment_joint_upfactorial_from_raw {
    private Moment_joint_upfactorial_from_raw() {}

    /**
     * Converts joint raw moments into joint upward-factorial moments.
     *
     * @param m flattened joint raw moments in row-major order, of extent
     *        dims[0]*...*dims[d-1]
     * @param dims extents of the array, dims[j] = n_j+1
     * @return flattened joint upward-factorial moments, same layout
     */
    public static double[] moment_joint_upfactorial_from_raw(double[] m, int[] dims) {
        return Moment_jointtrans.moment_jointtrans(m, dims, "upfactorial_from_raw");
    }
}
