/**
 * @file Cache Xi Terms via Fixed Point Iteration
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class Cache_xi_fp {
    private Cache_xi_fp() {}

    public static Ret.cacheXiFp cache_xi_fp(Matrix gamma, Matrix m, Matrix xi) {
        int n = gamma.getNumRows();
        int h = gamma.getNumCols();
        double tol = 1.0e-14;
        Matrix pi0 = new Matrix(1, n);
        pi0.fill(1.0 / (h + 1));
        Matrix pij = new Matrix(n, h);
        if (xi == null || xi.isEmpty()) {
            xi = new Matrix(1, h);
            for (int l = 0; l < h; l++) {
                double mean = Matrix.extractColumn(gamma, l, null).elementSum() / gamma.getNumRows();
                xi.set(l, m.get(l) / mean / (n + m.elementSum() - 1));
            }
        }
        int it = 1;
        while (it < 1.0e4) {
            Matrix pi0_1 = new Matrix(pi0);
            Matrix mul = pi0_1.mult(gamma, null);
            xi = new Matrix(m.getNumRows(), m.getNumCols());
            for (int i = 0; i < xi.getNumRows(); i++) {
                for (int j = 0; j < xi.getNumCols(); j++) {
                    xi.set(i, j, m.get(i, j) / mul.get(i, j));
                }
            }
            Matrix intermediate = gamma.elementMult(xi.repmat(n, 1), null);
            mul = gamma.mult(xi.repmat(n, 1).transpose());
            for (int i = 0; i < pij.getNumRows(); i++) {
                for (int j = 0; j < pij.getNumCols(); j++) {
                    pij.set(i, j, FastMath.abs(intermediate.get(i, j)) / FastMath.abs(1 + mul.get(i, j)));
                }
            }
            for (int i = 0; i < n; i++) {
                pi0.set(0, i, Maths.max(tol, 1 - pij.sumRows(i)));
            }
            double DELTA = 0.0;
            for (int i = 0; i < n; i++) {
                DELTA += FastMath.abs(1 - pi0.get(i) / pi0_1.get(i));
            }
            if (DELTA < tol) {
                break;
            }
            it++;
        }
        for (int i = 0; i < xi.getNumRows(); i++) {
            for (int j = 0; j < xi.getNumCols(); j++) {
                if (xi.get(i, j) < 0) {
                    xi.set(i, j, tol);
                }
            }
        }
        return new Ret.cacheXiFp(xi, pi0, pij, it);
    }
}
