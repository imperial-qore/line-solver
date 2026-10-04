package jline.lib.butools;

import jline.util.matrix.Matrix;

public final class ReducedMomsFromMoms {
    private ReducedMomsFromMoms() {}

    /**
     * Returns the reduced moments given the raw moments.
     *
     * The raw moments are: m_i = E(X^i)
     * The reduced moments are: r_i = m_i / i!
     *
     * @param m The list of raw moments (starting with the first moment)
     * @return The list of reduced moments
     */
    public static Matrix reducedMomsFromMoms(Matrix m) {
        Matrix rm = new Matrix(m.getNumRows(), m.getNumCols(), m.getNumRows() * m.getNumCols());
        double invFactorial = 1.0;

        for (int i = 0; i < m.length(); i++) {
            invFactorial /= (i + 1); // Calculate 1/(i+1)!
            rm.set(i, m.get(i) * invFactorial);
        }

        return rm;
    }

    /**
     * Returns the reduced moments given the raw moments (double[] version).
     *
     * The raw moments are: m_i = E(X^i)
     * The reduced moments are: r_i = m_i / i!
     *
     * @param m The list of raw moments (starting with the first moment)
     * @return The list of reduced moments
     */
    public static double[] ReducedMomsFromMoms(double[] m) {
        double[] rm = new double[m.length];
        double invFactorial = 1.0;

        for (int i = 0; i < m.length; i++) {
            invFactorial /= (i + 1); // Calculate 1/(i+1)!
            rm[i] = m[i] * invFactorial;
        }

        return rm;
    }
}
