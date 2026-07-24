package jline.lib.butools;

import jline.util.matrix.Matrix;

public final class JMomsFromJFactorialMoms {
    private JMomsFromJFactorialMoms() {}

    /**
     * Returns the lag-1 joint raw moments given the lag-1 joint factorial moments.
     *
     * @param jfm The matrix of joint factorial moments. The entry in row i and column j
     *            is f_{i,j}, i&gt;=1, j&gt;=1.
     * @return The matrix of joint raw moments. The entry in row i and column j
     *         is m_{i,j}, i&gt;=1, j&gt;=1.
     *
     * Reference: http://en.wikipedia.org/wiki/Factorial_moment
     */
    public static Matrix jMomsFromJFactorialMoms(Matrix jfm) {
        int s1 = jfm.getNumRows();
        int s2 = jfm.getNumCols();
        Matrix jmoms = Matrix.zeros(s1, s2);

        for (int i = 0; i < s1; i++) {
            for (int j = 0; j < s2; j++) {
                Matrix xRoots = new Matrix(1, i + 1, i + 1);
                for (int k = 0; k <= i; k++) xRoots.set(k, (double) k);
                Matrix xPoly = Poly.poly(xRoots); // i+2 coefficients
                double[] xCoeff = new double[i + 1];
                for (int k = 0; k <= i; k++) xCoeff[k] = xPoly.get(0, i - k);

                Matrix yRoots = new Matrix(1, j + 1, j + 1);
                for (int k = 0; k <= j; k++) yRoots.set(k, (double) k);
                Matrix yPoly = Poly.poly(yRoots); // j+2 coefficients
                double[] yCoeff = new double[j + 1];
                for (int k = 0; k <= j; k++) yCoeff[k] = yPoly.get(0, j - k);

                // eh = -(xCoeff' * yCoeff) (negated outer product)
                Matrix eh = Matrix.zeros(i + 1, j + 1);
                for (int r = 0; r <= i; r++) {
                    for (int c = 0; c <= j; c++) {
                        eh.set(r, c, -xCoeff[r] * yCoeff[c]);
                    }
                }

                // jmoms(i,j) = jfmoms(i,j) + trace(jmoms(1:i, 1:j) * eh')
                double traceVal = 0.0;
                for (int r = 0; r <= i; r++) {
                    for (int k = 0; k <= j; k++) {
                        traceVal += jmoms.get(r, k) * eh.get(r, k);
                    }
                }
                jmoms.set(i, j, jfm.get(i, j) + traceVal);
            }
        }

        return jmoms;
    }
}
