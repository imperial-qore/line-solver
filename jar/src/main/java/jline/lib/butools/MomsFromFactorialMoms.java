package jline.lib.butools;

import jline.util.matrix.Matrix;

public final class MomsFromFactorialMoms {
    private MomsFromFactorialMoms() {}

    /**
     * Returns the raw moments given the factorial moments.
     *
     * The raw moments are: m_i = E(X^i)
     * The factorial moments are: f_i = E(X(X-1)...(X-i+1))
     *
     * @param fm The list of factorial moments (starting with the first moment)
     * @return The list of raw moments
     *
     * Reference: http://en.wikipedia.org/wiki/Factorial_moment
     */
    public static Matrix MomsFromFactorialMoms(Matrix fm) {
        int n = fm.length();
        Matrix m = new Matrix(1, n, n);
        m.set(0, fm.get(0));

        for (int i = 1; i < n; i++) {
            // Compute polynomial with roots [0, 1, 2, ..., i]
            Matrix rootsMatrix = new Matrix(1, i + 1, i + 1);
            for (int k = 0; k <= i; k++) {
                rootsMatrix.set(k, (double) k);
            }

            // Get polynomial coefficients and negate
            // poly([0,1,...,i]) gives i+2 coefficients
            Matrix ehMatrix = Poly.poly(rootsMatrix);
            double[] eh = new double[ehMatrix.getNumCols()];
            for (int k = 0; k < ehMatrix.getNumCols(); k++) {
                eh[k] = -ehMatrix.get(0, k);
            }

            // MATLAB: eh(end-1:-1:2) for array of length i+2
            // 1-based: from index i+1 down to 2 = i elements
            // 0-based: from index i down to 1 = i elements
            double[] eh_coeff = new double[i];
            for (int k = 0; k < i; k++) {
                eh_coeff[k] = eh[i - k];
            }

            // Compute dot product: eh_coeff[0..i-1] dot m[0..i-1]
            double sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += eh_coeff[k] * m.get(k);
            }

            m.set(i, fm.get(i) + sum);
        }

        if (fm.getNumRows() > fm.getNumCols()) {
            return m.transpose();
        }
        return m;
    }
}
