package jline.lib.butools;

import jline.util.matrix.Matrix;

public final class MomsFromNormMoms {
    private MomsFromNormMoms() {}

    /**
     * Returns the raw moments given the normalized moments.
     *
     * The raw moments are: m_i = E(X^i)
     * The normalized moments are: n_i = m_i / (m_{i-1} * m_1)
     *
     * @param nm The list of normalized moments (starting with the first moment)
     * @return The list of raw moments
     */
    public static Matrix momsFromNormMoms(Matrix nm) {
        Matrix m = new Matrix(nm.getNumRows(), nm.getNumCols(), nm.getNumRows() * nm.getNumCols());

        for (int i = 0; i < nm.length(); i++) {
            if (i == 0) {
                m.set(i, nm.get(i));
            } else {
                m.set(i, m.get(0) * nm.get(i) * m.get(i - 1));
            }
        }

        return m;
    }
}
