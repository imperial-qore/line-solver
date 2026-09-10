package jline.api.moment;

import jline.util.matrix.Matrix;

/**
 * Triangle of the signed Stirling numbers of the first kind.
 *
 * <p>Builds the triangle of the signed Stirling numbers of the first kind
 * s(i,j), defined as the coefficients of x^j in the falling factorial
 *
 * <pre>
 *   sum_{j=0}^{i} s(i,j) x^j = x(x-1)(x-2)...(x-i+1)
 * </pre>
 *
 * <p>These numbers are the coefficients that convert power moments into
 * factorial moments. They relate to the Stirling cycle numbers via
 * s(i,j) = (-1)^(i-j) * sigma(i,j).
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003, eq. (10) and eq. (12).
 *
 * @since LINE 3.0
 */
public final class Moment_stirling1 {
    private Moment_stirling1() {}

    /**
     * Triangle of the signed Stirling numbers of the first kind s(i,j).
     *
     * @param n maximum order (n &gt;= 0)
     * @return (n+1)x(n+1) matrix with element (i,j) holding s(i,j) in the
     *         0-based notation of the reference; entries with j &gt; i are zero
     * @throws IllegalArgumentException if n is negative
     */
    public static Matrix moment_stirling1(int n) {
        Matrix sigma = Moment_stirlingcycle.moment_stirlingcycle(n);
        Matrix s = new Matrix(n + 1, n + 1);
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= i; j++) {
                double sign = ((i - j) % 2 == 0) ? 1.0 : -1.0;
                s.set(i, j, sign * sigma.get(i, j));
            }
        }
        return s;
    }
}
