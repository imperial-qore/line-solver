package jline.api.moment;

/**
 * Joint upward-factorial moments from joint negative-binomial moments.
 *
 * <p>Converts joint negative-binomial moments into joint upward-factorial
 * moments. Inverse of moment_joint_negbinomial_from_upfactorial.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003.
 *
 * @since LINE 3.0
 */
public final class Moment_joint_upfactorial_from_negbinomial {
    private Moment_joint_upfactorial_from_negbinomial() {}

    /**
     * Converts joint negative-binomial moments into joint upward-factorial moments.
     *
     * @param bm flattened joint negative-binomial moments in row-major order, of extent
     *        dims[0]*...*dims[d-1]
     * @param dims extents of the array, dims[j] = n_j+1
     * @return flattened joint upward-factorial moments, same layout
     */
    public static double[] moment_joint_upfactorial_from_negbinomial(double[] bm, int[] dims) {
        return Moment_jointtrans.moment_jointtrans(bm, dims, "upfactorial_from_negbinomial");
    }
}
