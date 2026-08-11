package jline.lib.butools;

import jline.util.matrix.Matrix;

public final class FactorialMomsFromMoms {
    private FactorialMomsFromMoms() {}

    /**
     * Returns the factorial moments given the raw moments.
     *
     * The raw moments are: m_i = E(X^i)
     * The factorial moments are: f_i = E(X(X-1)...(X-i+1))
     *
     * @param m The list of raw moments (starting with the first moment)
     * @return The list of factorial moments
     *
     * Reference: http://en.wikipedia.org/wiki/Factorial_moment
     */
    public static Matrix factorialMomsFromMoms(Matrix m) {
        int n = m.length();
        Matrix fm = new Matrix(m.getNumRows(), m.getNumCols(), m.length());

        for (int i = 0; i < n; i++) {
            // Build roots matrix [0, 1, ..., i]
            Matrix rootsMatrix = new Matrix(1, i + 1, i + 1);
            for (int k = 0; k <= i; k++) {
                rootsMatrix.set(k, (double) k);
            }

            // poly([0,1,...,i]) gives coefficients of (x-0)(x-1)...(x-i)
            // which has i+2 coefficients [c_{i+1}, c_i, ..., c_1, c_0]
            Matrix ehMatrix = Poly.poly(rootsMatrix);

            // MATLAB: eh = eh(end-1:-1:1)
            // Take from 2nd-to-last down to 1st element (1-based),
            // i.e., from index i down to 0 (0-based), giving i+1 elements
            double[] eh = new double[i + 1];
            for (int k = 0; k <= i; k++) {
                eh[k] = ehMatrix.get(0, i - k);
            }

            // dot product: eh * m(1:i+1)
            double sum = 0.0;
            for (int k = 0; k <= i; k++) {
                sum += eh[k] * m.get(k);
            }
            fm.set(i, sum);
        }

        return fm;
    }
}
