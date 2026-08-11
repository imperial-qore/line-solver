package jline.api.moment;

import jline.util.matrix.Matrix;

/**
 * Joint power (raw) moments from joint central moments.
 *
 * <p>Inverts {@link Moment_joint_central_from_raw} by the multi-index binomial
 * theorem,
 *
 * <pre>
 *   m_(i) = sum_(k&lt;=i) prod_j nchoosek(i_j,k_j) mu_j^(i_j-k_j) * mc_(k)
 * </pre>
 *
 * <p>The mean vector must be supplied separately, since the first-order central
 * moments are zero and carry no information on it.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003, Section 4.
 *
 * @since LINE 3.0
 */
public final class Moment_joint_raw_from_central {
    private Moment_joint_raw_from_central() {}

    /**
     * Converts joint central moments into joint power moments.
     *
     * @param mc flattened joint central moments in row-major order
     * @param dims extents of the array, dims[j] = n_j+1
     * @param mu means E[N_1],...,E[N_d], one per dimension
     * @return flattened joint power moments, same layout
     */
    public static double[] moment_joint_raw_from_central(double[] mc, int[] dims, double[] mu) {
        if (mu.length != dims.length) {
            throw new IllegalArgumentException("moment_joint_raw_from_central: The mean vector mu "
                    + "must have one entry per dimension of mc.");
        }
        double[] out = mc;
        for (int mode = 0; mode < dims.length; mode++) {
            int n = dims[mode] - 1;
            Matrix T = new Matrix(n + 1, n + 1);
            for (int i = 0; i <= n; i++) {
                for (int k = 0; k <= i; k++) {
                    T.set(i, k, Moment_binotrans.nchoosek(i, k) * Math.pow(mu[mode], i - k));
                }
            }
            out = Moment_tensortrans.moment_tensortrans(out, dims, T, mode);
        }
        return out;
    }
}
