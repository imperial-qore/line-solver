package jline.api.moment;

import jline.util.matrix.Matrix;

/**
 * Factorial moments from factorial cumulants.
 *
 * <p>Inverse of {@link Moment_factcumulant_from_factorial}.
 *
 * <p>Reference:
 * V. P. Leonov and A. N. Shiryaev. On a method of calculation of
 * semi-invariants. Theory of Probability and its Applications,
 * 4(3):319-329, 1959.
 *
 * @since LINE 3.0
 */
public final class Moment_factorial_from_factcumulant {
    private Moment_factorial_from_factcumulant() {}

    /**
     * Converts factorial cumulants into factorial moments.
     *
     * @param kappa column vector of length n+1 holding the factorial cumulants
     *              of order 0,...,n; element 0 is ignored
     * @return column vector of length n+1 holding f_0,...,f_n, with f_0 = 1
     */
    public static Matrix moment_factorial_from_factcumulant(Matrix kappa) {
        return Moment_raw_from_cumulant.moment_raw_from_cumulant(kappa);
    }
}
