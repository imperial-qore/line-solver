/**
 * @file Dense linear algebra over the complex field
 *
 * @since LINE 3.1.0
 */
package jline.util;

import org.apache.commons.math3.complex.Complex;

/**
 * The few complex operations the transform layer needs.
 *
 * WHY THIS EXISTS. A Laplace-Stieltjes transform is evaluated OFF the real axis
 * by everything that inverts it or locates its roots: the Abate-Whitt Euler sum
 * walks the line Re(s) = A/(2t), and the matrix transform int exp(Ut) dF(t) is
 * read off the spectrum of U, which is complex in general. {@link jline.util.matrix.ComplexMatrix}
 * carries complex data but offers no solve, and pulling in a full complex linear
 * algebra stack for an n x n system of the size a phase-type law has would be
 * out of proportion. Gaussian elimination with partial pivoting is enough and is
 * what this class provides.
 */
public final class ComplexOps {
    private ComplexOps() {}

    /**
     * Solves A x = b by Gaussian elimination with partial pivoting.
     *
     * @param A n x n coefficient matrix, consumed by value (a copy is taken)
     * @param b right-hand side of length n
     * @return the solution x
     */
    public static Complex[] solve(Complex[][] A, Complex[] b) {
        final int n = b.length;
        if (A.length != n) {
            throw new IllegalArgumentException("ComplexOps.solve: shape mismatch");
        }
        final Complex[][] a = new Complex[n][n];
        final Complex[] x = new Complex[n];
        for (int i = 0; i < n; i++) {
            if (A[i].length != n) {
                throw new IllegalArgumentException("ComplexOps.solve: matrix must be square");
            }
            System.arraycopy(A[i], 0, a[i], 0, n);
            x[i] = b[i];
        }
        for (int col = 0; col < n; col++) {
            int piv = col;
            double best = a[col][col].abs();
            for (int r = col + 1; r < n; r++) {
                final double mag = a[r][col].abs();
                if (mag > best) {
                    best = mag;
                    piv = r;
                }
            }
            if (piv != col) {
                final Complex[] t = a[piv];
                a[piv] = a[col];
                a[col] = t;
                final Complex tb = x[piv];
                x[piv] = x[col];
                x[col] = tb;
            }
            if (a[col][col].abs() == 0.0) {
                throw new ArithmeticException("ComplexOps.solve: singular matrix");
            }
            for (int r = col + 1; r < n; r++) {
                final Complex f = a[r][col].divide(a[col][col]);
                for (int c = col; c < n; c++) {
                    a[r][c] = a[r][c].subtract(f.multiply(a[col][c]));
                }
                x[r] = x[r].subtract(f.multiply(x[col]));
            }
        }
        for (int row = n - 1; row >= 0; row--) {
            Complex s = x[row];
            for (int c = row + 1; c < n; c++) {
                s = s.subtract(a[row][c].multiply(x[c]));
            }
            x[row] = s.divide(a[row][row]);
        }
        return x;
    }
}
