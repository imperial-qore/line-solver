package jline.api.cache;

import java.util.Collections;

import org.apache.commons.math3.util.FastMath;

import jline.util.matrix.Matrix;

public final class Cache_xi_bvh {
    private Cache_xi_bvh() {}

    /**
     * Computes the cache xi terms using the iterative method used in Gast-van Houdt, SIGMETRICS 2015.
     * This method calculates the xi values, which are important for understanding the distribution of items
     * in the cache. The script assumes (like the paper) that the access factors are monotone with the list index, so
     * it may not work with arbitrary (non-monotone) access costs.
     *
     * @param gamma Matrix representing the cache access factors.
     * @param m     Matrix representing the cache capacity vector.
     * @return Matrix containing the computed xi terms.
     */
    public static Matrix cache_xi_bvh(Matrix gamma, Matrix m) {
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

        while (z.sub(zOld).elementMaxAbs() > FastMath.pow(10, -12) * zOld.elementMaxAbs()) {
            zOld.setTo(z);
            Matrix temp = z.mult(pp).scale(n);
            for (int i = 0; i < h; i++) {
                Matrix a = temp.sub(pp.getRow(i + 1).scale(z.scale(n).get(0, i + 1)));
                double Fi = pp.getRow(i + 1).elementDiv(pp.getRow(i + 1).scale(n).add(a)).elementSum();

                double ziMin, ziMax;
                if (Fi > f.get(0, i)) {
                    ziMin = 0;
                    ziMax = 1;
                } else {
                    ziMin = 1;
                    ziMax = 2;
                    while (pp.getRow(i + 1).scale(ziMax).div(pp.getRow(i + 1).scale(n).scale(ziMax).add(a)).elementSum() < f.get(0, i)) {
                        ziMin = ziMax;
                        ziMax = ziMax * 2;
                    }
                }

                for (int x = 0; x < 50; x++) {
                    double zi = (ziMin + ziMax) / 2;
                    z.set(0, i + 1, zi);
                    if (pp.getRow(i + 1).scale(zi).div(pp.getRow(i + 1).scale(zi).scale(n).add(a)).elementSum() < f.get(0, i)) {
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
