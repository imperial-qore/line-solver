package jline.api.moment;

import jline.util.matrix.Matrix;

/**
 * Triangle of the Stirling cycle numbers.
 *
 * <p>Builds the triangle of the Stirling cycle numbers (unsigned Stirling numbers
 * of the first kind), sigma(i,j) = (-1)^(i-j) * s(i,j), obtained from the
 * recursion
 *
 * <pre>
 *   sigma(i,j) = (i-1)*sigma(i-1,j) + sigma(i-1,j-1)   for j &gt; 0
 *   sigma(0,0) = 1,   sigma(i,0) = 0 for i &gt; 0
 * </pre>
 *
 * <p>These numbers are the coefficients that convert power moments into
 * upward-factorial moments.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003, eq. (12).
 *
 * @since LINE 3.0
 */
public final class Moment_stirlingcycle {
    private Moment_stirlingcycle() {}

    /**
     * Triangle of the Stirling cycle numbers sigma(i,j).
     *
     * @param n maximum order (n &gt;= 0)
     * @return (n+1)x(n+1) matrix with element (i,j) holding sigma(i,j) in the
     *         0-based notation of the reference; entries with j &gt; i are zero
     * @throws IllegalArgumentException if n is negative
     */
    public static Matrix moment_stirlingcycle(int n) {
        if (n < 0) {
            throw new IllegalArgumentException("The maximum order n must be a nonnegative integer.");
        }
        Matrix sigma = new Matrix(n + 1, n + 1);
        sigma.set(0, 0, 1.0);
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= i; j++) {
                sigma.set(i, j, (i - 1) * sigma.get(i - 1, j) + sigma.get(i - 1, j - 1));
            }
        }
        return sigma;
    }
}
