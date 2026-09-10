package jline.lib.butools;

import jline.util.matrix.Matrix;

public final class HankelMomsFromMoms {
    private HankelMomsFromMoms() {}

    /**
     * Constructs a Hankel matrix from the given first column and last row.
     * A Hankel matrix has constant anti-diagonals.
     */
    private static Matrix hankel(double[] col, double[] row) {
        int n = col.length;
        int m = row.length;
        Matrix H = Matrix.zeros(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                if (i + j < n) {
                    H.set(i, j, col[i + j]);
                } else {
                    H.set(i, j, row[i + j - n + 1]);
                }
            }
        }
        return H;
    }

    /**
     * Returns the Hankel moments given the raw moments.
     *
     * The raw moments are: m_i = E(X^i)
     *
     * The ith Hankel moment is the determinant of matrix Delta_{i/2}, if i is even,
     * and it is the determinant of Delta^(1)_{(i+1)/2}, if i is odd.
     *
     * @param m The list of raw moments (starting with the first moment)
     * @return The list of Hankel moments
     *
     * Reference: http://en.wikipedia.org/wiki/Stieltjes_moment_problem
     */
    public static Matrix hankelMomsFromMoms(Matrix m) {
        int n = m.length();
        Matrix hm = new Matrix(m.getNumRows(), m.getNumCols(), n);

        // Extract raw moments as array for convenience
        double[] mArr = new double[n];
        for (int k = 0; k < n; k++) {
            mArr[k] = m.get(k);
        }

        for (int i = 0; i < n; i++) {
            if (i % 2 == 0) {
                // Even: N = i/2 + 1
                int N = i / 2 + 1;
                double[] col = new double[N];
                for (int k = 0; k < N; k++) {
                    col[k] = mArr[k];
                }
                double[] row = new double[N];
                for (int k = 0; k < N; k++) {
                    row[k] = mArr[N - 1 + k];
                }
                Matrix H = hankel(col, row);
                hm.set(i, H.det());
            } else {
                // Odd: N = (i+1)/2 + 1
                int N = (i + 1) / 2 + 1;
                double[] col = new double[N];
                col[0] = 1.0;
                for (int k = 1; k < N; k++) {
                    col[k] = mArr[k - 1];
                }
                double[] row = new double[N];
                for (int k = 0; k < N; k++) {
                    row[k] = mArr[N - 2 + k];
                }
                Matrix H = hankel(col, row);
                hm.set(i, H.det());
            }
        }

        return hm;
    }
}
