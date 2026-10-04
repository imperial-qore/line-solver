package jline.lib.butools;

import jline.util.matrix.Matrix;

public final class MomsFromReducedMoms {
    private MomsFromReducedMoms() {}

    /**
     * Returns the raw moments given the reduced moments.
     *
     * The raw moments are: m_i = E(X^i)
     * The reduced moments are: r_i = m_i / i!
     *
     * @param rm The list of reduced moments (starting with the first moment)
     * @return The list of raw moments
     */
    public static Matrix momsFromReducedMoms(Matrix rm) {
        Matrix m = new Matrix(rm.getNumRows(), rm.getNumCols(), rm.getNumRows() * rm.getNumCols());
        double factorial = 1.0;

        for (int i = 0; i < rm.length(); i++) {
            factorial *= (i + 1); // Calculate (i+1)!
            m.set(i, rm.get(i) * factorial);
        }

        return m;
    }
}
