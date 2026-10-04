/**
 * Moment computation for discrete-time MAPs.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Dmap_moment {
    private Dmap_moment() {}

    /**
     * Computes the k-th raw moment of the inter-arrival time of a discrete MAP.
     *
     *   E[T]   = pi * (I-D0)^{-1} * e
     *   E[T^2] = 2*pi*(I-D0)^{-2}*e - E[T]
     *   E[T^3] = 6*pi*(I-D0)^{-3}*e - 6*pi*(I-D0)^{-2}*e + E[T]
     */
    public static double dmap_moment(Matrix D0, Matrix D1, int order) {
        int n = D0.getNumRows();
        Matrix I = Matrix.eye(n);
        Matrix ImD0inv = I.add(-1.0, D0).inv();
        Matrix e = Matrix.ones(n, 1);
        Matrix pi = Dmap_pie.dmap_pie(D0, D1);

        switch (order) {
            case 1:
                return pi.mult(ImD0inv).mult(e).toDouble();
            case 2: {
                double m1 = pi.mult(ImD0inv).mult(e).toDouble();
                return 2.0 * pi.mult(ImD0inv.mult(ImD0inv)).mult(e).toDouble() - m1;
            }
            case 3: {
                double m1 = pi.mult(ImD0inv).mult(e).toDouble();
                Matrix ImD0inv2 = ImD0inv.mult(ImD0inv);
                return 6.0 * pi.mult(ImD0inv2.mult(ImD0inv)).mult(e).toDouble()
                        - 6.0 * pi.mult(ImD0inv2).mult(e).toDouble() + m1;
            }
            default:
                throw new IllegalArgumentException("Moments of order > 3 not implemented for DMAP");
        }
    }

    public static double dmap_moment(MatrixCell DMAP, int order) {
        return dmap_moment(DMAP.get(0), DMAP.get(1), order);
    }
}
