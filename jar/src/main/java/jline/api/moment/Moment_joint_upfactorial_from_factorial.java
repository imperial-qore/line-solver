package jline.api.moment;

/**
 * Joint upward-factorial moments from joint factorial moments.
 *
 * <p>Converts joint factorial moments into joint upward-factorial moments, by
 * the Lah numbers along every dimension.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003.
 *
 * @since LINE 3.0
 */
public final class Moment_joint_upfactorial_from_factorial {
    private Moment_joint_upfactorial_from_factorial() {}

    /**
     * Converts joint factorial moments into joint upward-factorial moments.
     *
     * @param f flattened joint factorial moments in row-major order, of extent
     *        dims[0]*...*dims[d-1]
     * @param dims extents of the array, dims[j] = n_j+1
     * @return flattened joint upward-factorial moments, same layout
     */
    public static double[] moment_joint_upfactorial_from_factorial(double[] f, int[] dims) {
        return Moment_jointtrans.moment_jointtrans(f, dims, "upfactorial_from_factorial");
    }
}
