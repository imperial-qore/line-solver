/**
 * @file Markovian Arrival Process asymptotic index of dispersion analysis
 *
 * Computes asymptotic index of dispersion for MAP counting processes.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_idc {
    private Map_idc() {}

    /**
     * Computes the asymptotic index of dispersion (IDC) for a Markovian Arrival Process (MAP).
     *
     * @param D0 the hidden transition matrix of the MAP
     * @param D1 the visible transition matrix of the MAP
     * @return the asymptotic index of dispersion
     */
    public static double map_idc(Matrix D0, Matrix D1) {
        Matrix e = new Matrix(D0.length(), 1, D0.length());
        for (int i = 0; i < D0.length(); i++) {
            e.set(i, 0, 1);
        }

        return 1 + 2 * (Map_lambda.map_lambda(D0, D1) - Map_pie.map_pie(D0, D1)
                .mult(Map_infgen.map_infgen(D0, D1).add(1.0, e.mult(Map_piq.map_piq(D0, D1)))
                        .inv()).mult(D1).mult(e).get(0));
    }

    /**
     * Computes the asymptotic index of dispersion (IDC) for a MAP
     * stored in a MatrixCell that contains the MAP's transition matrices.
     *
     * @param MAP a MatrixCell containing the transition matrices D0 and D1 of the MAP
     * @return the asymptotic index of dispersion
     */
    public static double map_idc(MatrixCell MAP) {
        return map_idc(MAP.get(0), MAP.get(1));
    }
}
