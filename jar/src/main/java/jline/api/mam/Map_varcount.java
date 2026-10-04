/**
 * @file Markovian Arrival Process counting variance analysis
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.Maths;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_varcount {
    private Map_varcount() {}

    /**
     * Variance of the counts in a MAP over a time period t.
     */
    public static double map_varcount(Matrix D0, Matrix D1, double t) {
        int n = D0.getNumRows();
        Matrix Q = Map_infgen.map_infgen(D0, D1);
        Matrix I = Matrix.eye(n);
        Matrix e = Matrix.ones(n, 1);
        Matrix piq = Map_piq.map_piq(D0, D1);
        double lambda = 1.0 / Map_mean.map_mean(D0, D1);
        Matrix IpiqQ = Matrix.ones(n, 1).mult(piq).add(-1.0, Q).inv();
        double PRE = (lambda - 2 * lambda * lambda + 2 * piq.mult(D1).mult(IpiqQ).mult(D1).mult(e).get(0, 0));
        double POST = piq.mult(D1).mult(I.sub(Maths.matrixExp(Q.scale(t)))).mult(IpiqQ.square()).mult(D1).mult(e).get(0, 0);
        return PRE * t - 2 * POST;
    }

    /**
     * Variance of the counts in a MAP over multiple time periods.
     */
    public static Matrix map_varcount(Matrix D0, Matrix D1, Matrix t) {
        Matrix result = new Matrix(t.getNumRows(), t.getNumCols());
        for (int i = 0; i < t.getNumElements(); i++) {
            result.set(i, map_varcount(D0, D1, t.get(i)));
        }
        return result;
    }

    public static double map_varcount(MatrixCell MAP, double t) {
        return map_varcount(MAP.get(0), MAP.get(1), t);
    }

    public static Matrix map_varcount(MatrixCell MAP, Matrix t) {
        return map_varcount(MAP.get(0), MAP.get(1), t);
    }
}
