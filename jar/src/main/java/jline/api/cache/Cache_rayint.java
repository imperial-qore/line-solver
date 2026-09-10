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
        int h = m.getNumElements();
        int ncols = gamma.getNumCols();
        boolean[] rowsToKeep = new boolean[gamma.getNumRows()];
        boolean[] allCols = new boolean[ncols];
        for (int i = 0; i < ncols; i++) {
            allCols[i] = true;
        }
        int n = 0;
        for (int i = 0; i < gamma.getNumRows(); i++) {
            rowsToKeep[i] = gamma.getRow(i).elementSum() > 0;
            if (rowsToKeep[i]) {
                n++;
            }
        }
        double mt = m.elementSum();

        if ((double) n == mt) {
            // Degenerate saddle: every item is cached, so the capacity equations force
            // every multiplier to infinity and cache_xi_bvh cannot converge. Take Z
            // from the exact recursion and report the limit, rather than iterating.
            System.err.println("The number of items equals the cache capacity.");
            double Zex = Cache_erec.cache_erec(gamma.getSlice(rowsToKeep, allCols), m).get(0);
            Matrix xiInf = new Matrix(1, h);
            for (int l = 0; l < h; l++) {
                xiInf.set(0, l, Double.POSITIVE_INFINITY);
            }
            return new Ret.cacheSpm(Zex, FastMath.log(Zex), xiInf);
        }

        // A list with no capacity has xi=0, which is a boundary of the Laplace integral
        // rather than a direction of it, so it must leave the expansion: kept, its
        // -sum_l log(sqrt(xi_l)) prefactor diverges and Z comes out far too large.
        // Dropping it is exact, since setting z_l=0 in the generating function removes
        // list l from E(m) and prod_l m_l! is unchanged because 0!=1.
        boolean[] colsToKeep = new boolean[ncols];
        int hk = 0;
        for (int l = 0; l < ncols; l++) {
            colsToKeep[l] = l < h && m.get(l) > 0;
            if (colsToKeep[l]) {
                hk++;
            }
        }
        Matrix xi = new Matrix(1, h); // dropped lists keep xi = 0
        if (hk == 0) {
            return new Ret.cacheSpm(1.0, 0.0, xi); // E(0)=1 and prod_l m_l!=1
        }

        Matrix gammaLocal = gamma.getSlice(rowsToKeep, colsToKeep);
        Matrix mk = new Matrix(1, hk);
        int[] keep = new int[hk];
        for (int l = 0, a = 0; l < ncols; l++) {
            if (colsToKeep[l]) {
                keep[a] = l;
                mk.set(0, a, m.get(l));
                a++;
            }
        }

        Matrix xik = Cache_xi_bvh.cache_xi_bvh(gammaLocal, mk);
        for (int a = 0; a < hk; a++) {
            xi.set(0, keep[a], xik.get(a));
        }

        Matrix S = new Matrix(n, 1);
        for (int k = 0; k < n; k++) {
            double Sk = 0.0;
            for (int l = 0; l < hk; l++) {
                Sk += gammaLocal.get(k, l) * xik.get(l);
            }
            S.set(k, 0, Sk);
        }

        // phi
        double phi = 0.0;
        for (int k = 0; k < n; k++) {
            phi += FastMath.log(1 + S.get(k));
        }
        phi -= xik.copy().log().mult(mk.copy().transpose()).elementSum();

        // A
        Matrix delta = Matrix.eye(hk);
        Matrix C = new Matrix(hk, hk);
        for (int j = 0; j < hk; j++) {
            for (int l = 0; l < hk; l++) {
                double C1 = 0.0;
                for (int k = 0; k < n; k++) {
                    C1 += gammaLocal.get(k, j) / (1 + S.get(k));
                }
                double C2 = 0.0;
                for (int k = 0; k < n; k++) {
                    C2 += gammaLocal.get(k, j) * gammaLocal.get(k, l) / FastMath.pow(1 + S.get(k), 2);
                }
                C.set(j, l, delta.get(j, l) * C1 - xik.get(j) * C2);
            }
        }

        // Z
        double Z = FastMath.exp(phi) * FastMath.pow(Math.sqrt(2 * FastMath.PI), -hk)
                * mk.fact().elementMult() / xik.sqrt().elementMult()
                / FastMath.sqrt(C.det());
        double lZ = -hk * FastMath.log(Math.sqrt(2 * FastMath.PI)) + phi
                + mk.factln().elementSum() - xik.sqrt().log().elementSum()
                - FastMath.log(Math.sqrt(C.det()));

        return new Ret.cacheSpm(Z, lZ, xi);
    }
}
