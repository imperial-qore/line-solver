package jline.api.moment;

import jline.util.matrix.Matrix;

/**
 * Joint central moments about a given mean vector.
 *
 * <p>Converts the joint power (raw) moments of a random vector into the joint
 * central moments about the supplied means, by the multi-index binomial
 * theorem,
 *
 * <pre>
 *   mc_(i) = sum_(k&lt;=i) prod_j (-1)^(i_j-k_j) nchoosek(i_j,k_j)
 *            mu_j^(i_j-k_j) * m_(k)
 * </pre>
 *
 * <p>which is separable, with a different shift per dimension. Supplying the
 * means makes the conversion applicable when the array does not carry the
 * first-order entries.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003, Section 4.
 *
 * @since LINE 3.0
 */
public final class Moment_joint_central_from_raw_mean {
    private Moment_joint_central_from_raw_mean() {}

    /**
     * Converts joint power moments into joint central moments about mu.
     *
     * @param m flattened joint power moments in row-major order
     * @param dims extents of the array, dims[j] = n_j+1
     * @param mu means E[N_1],...,E[N_d], one per dimension
     * @return flattened joint central moments, same layout
     */
    public static double[] moment_joint_central_from_raw_mean(double[] m, int[] dims, double[] mu) {
        if (mu.length != dims.length) {
            throw new IllegalArgumentException("moment_joint_central_from_raw_mean: The mean "
                    + "vector mu must have one entry per dimension of m.");
        }
        double[] out = m;
        for (int mode = 0; mode < dims.length; mode++) {
            int n = dims[mode] - 1;
            Matrix T = new Matrix(n + 1, n + 1);
            for (int i = 0; i <= n; i++) {
                for (int k = 0; k <= i; k++) {
                    T.set(i, k, Moment_binotrans.nchoosek(i, k) * Math.pow(-mu[mode], i - k));
                }
            }
            out = Moment_tensortrans.moment_tensortrans(out, dims, T, mode);
        }
        return out;
    }
}
