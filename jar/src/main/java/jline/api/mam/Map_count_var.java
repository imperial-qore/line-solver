/**
 * @file Markovian Arrival Process counting process variance analysis
 *
 * Computes variance of MAP counting processes over specified time intervals using matrix
 * exponential methods. Essential for analyzing variability in arrival patterns and burstiness.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.Maths;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

public final class Map_count_var {
    private Map_count_var() {}

    /**
     * Computes the variance of the counting process over a specified interval length for a given Markovian Arrival Process (MAP).
     *
     * @param MAP The Markovian Arrival Process stored in a MatrixCell.
     * @param t   The length of the interval over which to compute the variance.
     * @return The variance of the counting process over the interval `t`.
     */
    public static double map_count_var(MatrixCell MAP, double t) {
        return map_count_var(MAP, new double[]{t})[0];
    }

    /**
     * Computes the variance of the counting process over multiple specified interval lengths for a given Markovian Arrival Process (MAP).
     *
     * @param MAP The Markovian Arrival Process stored in a MatrixCell.
     * @param t   An array of interval lengths over which to compute the variance.
     * @return An array of doubles, where each element represents the variance of the counting process.
     */
    public static double[] map_count_var(MatrixCell MAP, double[] t) {
        double[] ret = new double[t.length];
        double l = Map_lambda.map_lambda(MAP);
        Matrix D = Map_infgen.map_infgen(MAP);
        Matrix theta = Map_piq.map_piq(MAP);
        Matrix D1 = MAP.get(1);
        int n = MAP.get(0).length();
        Matrix e = Matrix.ones(n, 1);
        Matrix tmp = e.mult(theta).sub(D).inv();
        Matrix c = theta.mult(D1).mult(tmp);
        Matrix d = tmp.mult(D1).mult(e);
        Matrix ll = theta.mult(D1).mult(e);
        double l2 = 2 * FastMath.pow(l, 2);
        for (int i = 0; i < t.length; i++) {
            double t_i = t[i];
            Matrix expmDt = Maths.matrixExp(D.scale(t_i));
            Matrix term1 = ll.sub(Matrix.createLike(ll).fill(l2)).add(c.mult(MAP.get(1)).mult(e).scale(2.0)).scale(t_i);
            Matrix term2 = c.mult(Matrix.eye(n).sub(expmDt)).mult(d).scale(2.0);
            ret[i] = term1.sub(term2).value();
        }
        return ret;
    }
}
