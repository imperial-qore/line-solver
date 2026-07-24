package jline.api.moment;

/**
 * Joint negative-binomial moments from joint upward-factorial moments.
 *
 * <p>Converts joint upward-factorial moments into the joint negative-binomial
 * moments b-_(i_1,...,i_d) = E[prod_j nchoosek(N_j+i_j-1,i_j)] = f+_(i) /
 * prod_j (i_j!).
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003.
 *
 * @since LINE 3.0
 */
public final class Moment_joint_negbinomial_from_upfactorial {
    private Moment_joint_negbinomial_from_upfactorial() {}

    /**
     * Converts joint upward-factorial moments into joint negative-binomial moments.
     *
     * @param fp flattened joint upward-factorial moments in row-major order, of extent
     *        dims[0]*...*dims[d-1]
     * @param dims extents of the array, dims[j] = n_j+1
     * @return flattened joint negative-binomial moments, same layout
     */
    public static double[] moment_joint_negbinomial_from_upfactorial(double[] fp, int[] dims) {
        return Moment_jointtrans.moment_jointtrans(fp, dims, "negbinomial_from_upfactorial");
    }
}
