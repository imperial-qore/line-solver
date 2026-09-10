/**
 * @file Cache Xi Terms via Gast-van Houdt Method
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import java.util.Collections;

import org.apache.commons.math3.util.FastMath;

import jline.util.matrix.Matrix;

public final class Cache_xi_iter {
    private Cache_xi_iter() {}

    /**
     * Computes the cache xi terms using the iterative method (Gast-van Houdt, SIGMETRICS 2015).
     */
    public static Matrix cache_xi_iter(Matrix gamma, Matrix m) {
        int n = gamma.getNumRows();
        Matrix f = m.scale(1.0 / n);
        int h = f.getNumCols();

        Matrix pp = new Matrix(h + 1, n);
        pp.setRow(0, Matrix.ones(1, n));
        for (int i = 0; i < h; i++) {
            pp.setRow(i + 1, gamma.getColumn(i).transpose());
        }

        Matrix zOld = Matrix.zeros(1, h + 1);
        Matrix z = Matrix.ones(1, h + 1);

        while (z.sub(zOld).elementMaxAbs() > FastMath.pow(10.0, -12) * zOld.elementMaxAbs()) {
            zOld.setTo(z);
            Matrix temp = z.mult(pp).scale((double) n);
            for (int i = 0; i < h; i++) {
                Matrix a = temp.sub(pp.getRow(i + 1).scale(z.scale((double) n).get(0, i + 1)));
                double Fi = pp.getRow(i + 1).elementDiv(pp.getRow(i + 1).scale((double) n).add(a)).elementSum();

                double ziMin;
                double ziMax;
                if (Fi > f.get(0, i)) {
                    ziMin = 0.0;
                    ziMax = 1.0;
                } else {
                    ziMin = 1.0;
                    ziMax = 2.0;
                    while (pp.getRow(i + 1).scale(ziMax)
                            .div(pp.getRow(i + 1).scale((double) n).scale(ziMax).add(a))
                            .elementSum() < f.get(0, i)) {
                        ziMin = ziMax;
                        ziMax = ziMax * 2;
                    }
                }

                for (int x = 0; x < 50; x++) {
                    double zi = (ziMin + ziMax) / 2;
                    z.set(0, i + 1, zi);
                    if (pp.getRow(i + 1).scale(zi)
                            .div(pp.getRow(i + 1).scale(zi).scale((double) n).add(a))
                            .elementSum() < f.get(0, i)) {
                        ziMin = zi;
                    } else {
                        ziMax = zi;
                    }
                }
            }
        }

        z.removeCols(Collections.singleton(0));
        return z;
    }
}
