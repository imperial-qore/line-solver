package jline.api.moment;

/**
 * Joint factorial moments from joint factorial cumulants.
 *
 * <p>Inverse of {@link Moment_joint_factcumulant_from_factorial}.
 *
 * <p>Reference:
 * V. P. Leonov and A. N. Shiryaev. On a method of calculation of
 * semi-invariants. Theory of Probability and its Applications,
 * 4(3):319-329, 1959.
 *
 * @since LINE 3.0
 */
public final class Moment_joint_factorial_from_factcumulant {
    private Moment_joint_factorial_from_factcumulant() {}

    /**
     * Converts joint factorial cumulants into joint factorial moments.
     *
     * @param kappa flattened joint factorial cumulants in row-major order;
     *              element 0 is ignored
     * @param dims extents of the array, dims[j] = n_j+1
     * @return flattened joint factorial moments, same layout, element 0 being 1
     */
    public static double[] moment_joint_factorial_from_factcumulant(double[] kappa, int[] dims) {
        return Moment_joint_raw_from_cumulant.moment_joint_raw_from_cumulant(kappa, dims);
    }
}
