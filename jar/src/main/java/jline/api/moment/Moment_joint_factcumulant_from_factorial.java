package jline.api.moment;

/**
 * Joint factorial cumulants from joint factorial moments.
 *
 * <p>The joint factorial cumulants are the coefficients of the logarithm of the
 * joint probability generating function expanded about z = (1,...,1),
 *
 * <pre>
 *   log E[prod_j z_j^(N_j)] = sum_(a != 0) kappa_a prod_j (z_j-1)^(a_j) / a_j!
 * </pre>
 *
 * <p>and stand to the joint factorial moments exactly as the joint cumulants
 * stand to the joint power moments, so the same recursion applies. For a
 * multivariate Poisson vector with independent components every joint factorial
 * cumulant of order two or more vanishes; for the per-class counts of a marked
 * MAP they measure the departure from independent Poisson marking.
 *
 * <p>Reference:
 * V. P. Leonov and A. N. Shiryaev. On a method of calculation of
 * semi-invariants. Theory of Probability and its Applications,
 * 4(3):319-329, 1959.
 *
 * @since LINE 3.0
 */
public final class Moment_joint_factcumulant_from_factorial {
    private Moment_joint_factcumulant_from_factorial() {}

    /**
     * Converts joint factorial moments into joint factorial cumulants.
     *
     * @param f flattened joint factorial moments in row-major order, element 0
     *          being 1
     * @param dims extents of the array, dims[j] = n_j+1
     * @return flattened joint factorial cumulants, same layout
     */
    public static double[] moment_joint_factcumulant_from_factorial(double[] f, int[] dims) {
        return Moment_joint_cumulant_from_raw.moment_joint_cumulant_from_raw(f, dims);
    }
}
