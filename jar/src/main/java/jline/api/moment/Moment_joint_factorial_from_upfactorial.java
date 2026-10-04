package jline.api.moment;

/**
 * Joint factorial moments from joint upward-factorial moments.
 *
 * <p>Converts joint upward-factorial moments into joint factorial moments, by
 * the signed Lah numbers along every dimension.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003.
 *
 * @since LINE 3.0
 */
public final class Moment_joint_factorial_from_upfactorial {
    private Moment_joint_factorial_from_upfactorial() {}

    /**
     * Converts joint upward-factorial moments into joint factorial moments.
     *
     * @param fp flattened joint upward-factorial moments in row-major order, of extent
     *        dims[0]*...*dims[d-1]
     * @param dims extents of the array, dims[j] = n_j+1
     * @return flattened joint factorial moments, same layout
     */
    public static double[] moment_joint_factorial_from_upfactorial(double[] fp, int[] dims) {
        return Moment_jointtrans.moment_jointtrans(fp, dims, "factorial_from_upfactorial");
    }
}
