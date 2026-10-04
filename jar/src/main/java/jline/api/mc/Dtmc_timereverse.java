package jline.api.mc;

import jline.util.matrix.Matrix;

public final class Dtmc_timereverse {
    private Dtmc_timereverse() {}

    /**
     * Compute the infinitesimal generator of the time-reversed DTMC.
     *
     * @param P Infinitesimal generator of the DTMC
     * @return Infinitesimal generator of the time-reversed DTMC
     */
    public static Matrix dtmc_timereverse(Matrix P) {
        Matrix pie = Dtmc_solve.dtmc_solve(P);
        Matrix Prev = new Matrix(P.getNumCols(), P.getNumRows());

        for (int i = 0; i < P.getNumRows(); i++) {
            for (int j = 0; j < P.getNumCols(); j++) {
                Prev.set(i, j, P.get(i, j) * pie.get(i) / pie.get(j));
            }
        }
        return Prev.transpose();
    }
}
