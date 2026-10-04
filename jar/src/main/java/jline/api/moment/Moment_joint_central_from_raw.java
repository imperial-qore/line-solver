package jline.api.moment;

/**
 * Joint central moments from joint power (raw) moments.
 *
 * <p>Converts the joint power moments of a random vector (N_1,...,N_d) into the
 * joint central moments mc_(i_1,...,i_d) = E[prod_j (N_j - E N_j)^(i_j)]. The
 * means mu_j = m_(e_j) are read off the array itself, so every dimension must
 * carry at least the first order. The entry of multi-order e_j+e_l is the
 * covariance of N_j and N_l.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003, Section 4.
 *
 * @since LINE 3.0
 */
public final class Moment_joint_central_from_raw {
    private Moment_joint_central_from_raw() {}

    /**
     * Converts joint power moments into joint central moments.
     *
     * @param m flattened joint power moments in row-major order
     * @param dims extents of the array, dims[j] = n_j+1, every dims[j] &gt;= 2
     * @return flattened joint central moments, same layout
     */
    public static double[] moment_joint_central_from_raw(double[] m, int[] dims) {
        int d = dims.length;
        for (int j = 0; j < d; j++) {
            if (dims[j] < 2) {
                throw new IllegalArgumentException("moment_joint_central_from_raw: The means "
                        + "m_(e_j) are required for this conversion, hence every dimension of m "
                        + "must have at least 2 elements.");
            }
        }
        double[] mu = new double[d];
        for (int j = 0; j < d; j++) {
            int stride = 1;
            for (int i = j + 1; i < d; i++) {
                stride *= dims[i];
            }
            mu[j] = m[stride];
        }
        return Moment_joint_central_from_raw_mean.moment_joint_central_from_raw_mean(m, dims, mu);
    }
}
