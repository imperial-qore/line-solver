package jline.api.moment;

import jline.util.matrix.Matrix;

/**
 * Factorial cumulants from factorial moments.
 *
 * <p>The factorial cumulants of a discrete random variable N are the
 * coefficients of the logarithm of the probability generating function expanded
 * about z = 1,
 *
 * <pre>
 *   log E[z^N] = sum_{n&gt;=1} kappa_n (z-1)^n / n!
 * </pre>
 *
 * <p>They stand to the factorial moments exactly as the cumulants stand to the
 * power moments, so the same recursion applies. For a Poisson variable of rate
 * lambda all factorial cumulants beyond the first vanish, which makes them the
 * natural measure of departure from Poisson behaviour in the counting process
 * of a MAP.
 *
 * <p>Reference:
 * V. P. Leonov and A. N. Shiryaev. On a method of calculation of
 * semi-invariants. Theory of Probability and its Applications,
 * 4(3):319-329, 1959.
 *
 * @since LINE 3.0
 */
public final class Moment_factcumulant_from_factorial {
    private Moment_factcumulant_from_factorial() {}

    /**
     * Converts factorial moments into factorial cumulants.
     *
     * @param f column vector of length n+1 holding f_0,...,f_n, with f_0 = 1
     * @return column vector of length n+1 holding the factorial cumulants of
     *         order 0,...,n, element 0 being 0
     */
    public static Matrix moment_factcumulant_from_factorial(Matrix f) {
        return Moment_cumulant_from_raw.moment_cumulant_from_raw(f);
    }
}
