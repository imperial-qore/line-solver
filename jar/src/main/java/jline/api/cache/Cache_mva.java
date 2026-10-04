/**
 * @file Cache Analysis via Mean Value Analysis
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import jline.io.Ret;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.Collections;

public final class Cache_mva {
    private Cache_mva() {}

    /**
     * Exact recursive solution of the caching model.
     */
    public static Ret.cacheMVA cache_mva(Matrix gamma, Matrix m) {
        int n = gamma.getNumRows();
        int h = gamma.getNumCols();
        Matrix SS = new Matrix(0, 0);
        for (int l = 0; l < h; l++) {
            Matrix arg = new Matrix((int) m.get(l) + 1, 1);
            int i = 1;
            while (i <= m.get(l) + 1) {
                arg.set(i - 1, (double) i);
                i++;
            }
            SS = Matrix.cartesian(SS, arg);
        }
        for (int i = 0; i < SS.getNumRows(); i++) {
            for (int j = 0; j < SS.getNumCols(); j++) {
                SS.set(i, j, SS.get(i, j) - 1);
            }
        }
        Matrix pi = new Matrix(SS.getNumRows(), n);
        Matrix[] pij = new Matrix[SS.getNumRows()];
        for (int i = 0; i < SS.getNumRows(); i++) {
            pij[i] = new Matrix(n, h);
        }
        Matrix x = new Matrix(1, h);
        int E = 1;
        for (int s = 0; s < SS.getNumRows(); s++) {
            Matrix mcur = Matrix.extractRows(SS, s, s + 1, null);
            for (int l = 0; l < h; l++) {
                Matrix mcur_l = Matrix.oner(mcur, new ArrayList<Integer>(Collections.singletonList(l)));
                int s_l = Matrix.matchrow(SS, mcur_l);
                if (s_l >= 0) {
                    Matrix one_pi = Matrix.extractRows(pi, s_l, s_l + 1, null);
                    for (int i = 0; i < one_pi.getNumCols(); i++) {
                        one_pi.set(i, 1 - one_pi.get(i));
                    }
                    x.set(l, mcur.get(l) / Matrix.extractColumn(gamma, l, null).transpose().mult(one_pi.transpose()).get(0));
                    Matrix elMult = Matrix.extractColumn(gamma, l, null).transpose().elementMult(one_pi, null);
                    for (int i = 0; i < n; i++) {
                        pij[s].set(i, l, elMult.get(i) * x.get(l));
                        pi.set(s, i, pi.get(s, i) + pij[s].get(i, l));
                    }
                }
            }
        }
        int s = Matrix.matchrow(SS, m);
        pi = Matrix.extractRows(pi, s, s + 1, null).transpose();
        Matrix newPij = new Matrix(n, h);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < h; j++) {
                newPij.set(i, j, pij[s].get(i, j));
            }
        }
        Matrix pi0 = new Matrix(pi.getNumRows(), pi.getNumCols());
        for (int i = 0; i < pi0.getNumRows(); i++) {
            for (int j = 0; j < pi0.getNumCols(); j++) {
                pi0.set(i, j, 1 - pi.get(i, j));
            }
        }
        Matrix u = new Matrix(n, h);
        for (int l = 0; l < h; l++) {
            for (int k = 0; k < n; k++) {
                u.set(k, l, x.get(l) * gamma.get(k, l));
            }
        }
        return new Ret.cacheMVA(pi, pi0, newPij, x, u, E);
    }
}
