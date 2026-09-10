package jline.api.moment;

/**
 * Joint raw moments from joint upward-factorial moments.
 *
 * <p>Converts joint upward-factorial moments into joint power (raw) moments, by
 * the signed Stirling numbers of the second kind along every dimension.
 * Inverse of moment_joint_upfactorial_from_raw.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003.
 *
 * @since LINE 3.0
 */
public final class Moment_joint_raw_from_upfactorial {
    private Moment_joint_raw_from_upfactorial() {}

    /**
     * Converts joint upward-factorial moments into joint raw moments.
     *
     * @param fp flattened joint upward-factorial moments in row-major order, of extent
     *        dims[0]*...*dims[d-1]
     * @param dims extents of the array, dims[j] = n_j+1
     * @return flattened joint raw moments, same layout
     */
    public static double[] moment_joint_raw_from_upfactorial(double[] fp, int[] dims) {
        return Moment_jointtrans.moment_jointtrans(fp, dims, "raw_from_upfactorial");
    }
}
