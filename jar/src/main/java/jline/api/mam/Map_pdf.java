/**
 * @file Markovian Arrival Process probability density function computation
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.Maths;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_pdf {
    private Map_pdf() {}

    /**
     * Computes the PDF of a MAP at specified time points.
     */
    public static double[] map_pdf(MatrixCell MAP, double[] tset) {
        return map_pdf(MAP.get(0), MAP.get(1), tset);
    }

    public static double[] map_pdf(Matrix D0, Matrix D1, double[] tset) {
        Matrix pi = Map_pie.map_pie(D0, D1);
        Matrix e = Matrix.ones(D1.getNumRows(), 1);
        Matrix minusD0 = D0.scale(-1.0);

        double[] result = new double[tset.length];
        for (int i = 0; i < tset.length; i++) {
            double t = tset[i];
            if (t < 0) {
                result[i] = 0.0;
            } else {
                // t = 0 IS NOT A SPECIAL CASE. f(0) = pie(-D0)e, which is
                // positive for any law with mass at the origin (an exponential,
                // a hyperexponential) and zero only for one without (an
                // Erlang). Hardcoding 0 here made Java disagree with MATLAB,
                // Python and C++ on every such law, and the error reached a
                // caller as soon as a density was integrated from 0.
                Matrix D0t = D0.scale(t);
                Matrix expD0t = Maths.matrixExp(D0t);
                Matrix temp = pi.mult(expD0t).mult(minusD0).mult(e);
                result[i] = temp.get(0, 0);
            }
        }
        return result;
    }

    public static double map_pdf(MatrixCell MAP, double t) {
        return map_pdf(MAP, new double[]{t})[0];
    }

    public static double map_pdf(Matrix D0, Matrix D1, double t) {
        return map_pdf(D0, D1, new double[]{t})[0];
    }
}
