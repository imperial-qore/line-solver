package jline.api.moment;

import jline.util.matrix.Matrix;

/**
 * Triangle of the Lah numbers.
 *
 * <p>Builds the triangle of the Lah numbers L(i,j) = (i!/j!)*nchoosek(i-1,j-1),
 * which link the factorial moments to the upward-factorial moments. The triangle
 * is built from the equivalent recursion
 *
 * <pre>
 *   L(i,j) = L(i-1,j-1) + (i+j-1)*L(i-1,j)   for j &gt; 0
 *   L(0,0) = 1,   L(i,0) = 0 for i &gt; 0
 * </pre>
 *
 * <p>which avoids the overflow of the explicit factorial form for large orders.
 *
 * <p>References:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003, Section 4.
 * I. Lah. Eine neue Art von Zahlen, ihre Eigenschaften und Anwendung in der
 * mathematischen Statistik. Mitteilungsbl. Math. Statist., 7:203-212, 1955.
 *
 * @since LINE 3.0
 */
public final class Moment_lah {
    private Moment_lah() {}

    /**
     * Triangle of the Lah numbers L(i,j).
     *
     * @param n maximum order (n &gt;= 0)
     * @return (n+1)x(n+1) matrix with element (i,j) holding L(i,j) in the
     *         0-based notation of the reference; entries with j &gt; i are zero
     * @throws IllegalArgumentException if n is negative
     */
    public static Matrix moment_lah(int n) {
        if (n < 0) {
            throw new IllegalArgumentException("The maximum order n must be a nonnegative integer.");
        }
        Matrix L = new Matrix(n + 1, n + 1);
        L.set(0, 0, 1.0);
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= i; j++) {
                L.set(i, j, L.get(i - 1, j - 1) + (i + j - 1) * L.get(i - 1, j));
            }
        }
        return L;
    }
}
