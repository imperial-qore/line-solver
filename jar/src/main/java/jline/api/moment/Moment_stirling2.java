package jline.api.moment;

import jline.util.matrix.Matrix;

/**
 * Triangle of the Stirling numbers of the second kind.
 *
 * <p>Builds the triangle of the Stirling numbers of the second kind S(i,j),
 * implicitly defined by the expansion of a power into falling factorials
 *
 * <pre>
 *   x^i = sum_{j=0}^{i} S(i,j) x(x-1)(x-2)...(x-j+1)
 * </pre>
 *
 * <p>and computed from the recursion
 *
 * <pre>
 *   S(i,j) = j*S(i-1,j) + S(i-1,j-1)   for j &gt; 0
 *   S(0,0) = 1,   S(i,0) = 0 for i &gt; 0
 * </pre>
 *
 * <p>These numbers are the coefficients that convert factorial moments back into
 * power moments.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003, eq. (11).
 *
 * @since LINE 3.0
 */
public final class Moment_stirling2 {
    private Moment_stirling2() {}

    /**
     * Triangle of the Stirling numbers of the second kind S(i,j).
     *
     * @param n maximum order (n &gt;= 0)
     * @return (n+1)x(n+1) matrix with element (i,j) holding S(i,j) in the
     *         0-based notation of the reference; entries with j &gt; i are zero
     * @throws IllegalArgumentException if n is negative
     */
    public static Matrix moment_stirling2(int n) {
        if (n < 0) {
            throw new IllegalArgumentException("The maximum order n must be a nonnegative integer.");
        }
        Matrix S = new Matrix(n + 1, n + 1);
        S.set(0, 0, 1.0);
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= i; j++) {
                S.set(i, j, j * S.get(i - 1, j) + S.get(i - 1, j - 1));
            }
        }
        return S;
    }
}
