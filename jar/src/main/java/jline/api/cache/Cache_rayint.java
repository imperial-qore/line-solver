package jline.api.cache;

import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Cache_rayint {
    private Cache_rayint() {}

    /**
     * Approximate the normalizing constant of the cache steady state distribution using the ray method.
     *
     * @param gamma Matrix representing the cache access factors.
     * @param m     Matrix representing the cache capacity vector.
     * @return cacheSpm: the approximated normalizing constant (Z and its logarithm lE) and the xi terms.
     */
    public static Ret.cacheSpm cache_rayint(Matrix gamma, Matrix m) {
        boolean[] rowsToKeep = new boolean[gamma.getNumRows()];
        boolean[] colsToKeep = new boolean[gamma.getNumCols()];
        for (int i = 0; i < gamma.getNumCols(); i++) {
            colsToKeep[i] = true;
        }
        for (int i = 0; i < gamma.getNumRows(); i++) {
            rowsToKeep[i] = gamma.getRow(i).elementSum() > 0;
        }
        gamma = gamma.getSlice(rowsToKeep, colsToKeep);

        int h = m.getNumElements();
        int n = gamma.getNumRows();
        double mt = m.elementSum();

        if ((double) n == mt) {
            System.err.println("The number of items equals the cache capacity.");
        }

        Matrix xi = Cache_xi_bvh.cache_xi_bvh(gamma, m);

        Matrix S = new Matrix(n, 1);
        for (int k = 0; k < n; k++) {
            double Sk = 0.0;
            for (int l = 0; l < h; l++) {
                Sk += gamma.get(k, l) * xi.get(l);
            }
            S.set(k, 0, Sk);
        }

        // phi
        double phi = 0.0;
        for (int k = 0; k < n; k++) {
            phi += FastMath.log(1 + S.get(k));
        }
        phi -= xi.copy().log().mult(m.copy().transpose()).elementSum();

        // A
        Matrix delta = Matrix.eye(h);
        Matrix C = new Matrix(h, h);
        for (int j = 0; j < h; j++) {
            for (int l = 0; l < h; l++) {
                double C1 = 0.0;
                for (int k = 0; k < n; k++) {
                    C1 += gamma.get(k, j) / (1 + S.get(k));
                }
                double C2 = 0.0;
                for (int k = 0; k < n; k++) {
                    C2 += gamma.get(k, j) * gamma.get(k, l) / FastMath.pow(1 + S.get(k), 2);
                }
                C.set(j, l, delta.get(j, l) * C1 - xi.get(j) * C2);
            }
        }

        // Z
        double Z = FastMath.exp(phi) * FastMath.pow(Math.sqrt(2 * FastMath.PI), -h)
                * m.fact().elementMult() / xi.sqrt().elementMult()
                / FastMath.sqrt(C.det());
        double lZ = -h * FastMath.log(Math.sqrt(2 * FastMath.PI)) + phi
                + m.factln().elementSum() - xi.sqrt().log().elementSum()
                - FastMath.log(Math.sqrt(C.det()));

        return new Ret.cacheSpm(Z, lZ, xi);
    }
}
