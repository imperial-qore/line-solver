/**
 * @file Cache Miss Analysis via SPM
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import java.util.Collections;

import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Cache_miss_spm {
    private Cache_miss_spm() {}

    public static Ret.cacheMissSpm cache_miss_spm(Matrix gamma, Matrix m, MatrixCell lambda) {
        Matrix ma = m.copy();
        ma.set(0, 0, ma.value() + 1);

        double lE = Cache_spm.cache_spm(gamma, m).lZ;
        double lEa = Cache_spm.cache_spm(gamma, ma).lZ;

        double M = FastMath.exp(lEa - lE);

        int u = lambda.size();
        int n = lambda.get(0).getNumRows();
        double[] pi0 = new double[n];
        double[] MU = new double[u];
        double[] lE1 = new double[n];

        for (int k = 0; k < n; k++) {
            if (gamma.getRow(k).elementSum() > 0) {
                Matrix subGamma = gamma.copy();
                subGamma.removeRows(Collections.singleton(k));

                lE1[k] = Cache_spm.cache_spm(subGamma, m).lZ;
                pi0[k] = FastMath.exp(lE1[k] - lE);

                for (int v = 0; v < u; v++) {
                    MU[v] += lambda.get(v).get(k, 0) * pi0[k];
                }
            }
        }

        double[] MI = new double[n];
        for (int k = 0; k < n; k++) {
            if (gamma.getRow(k).elementSum() > 0) {
                MI[k] = lambda.cellsum(k, 0) * FastMath.exp(lE1[k] - lE);
            } else {
                MI[k] = 0.0;
            }
        }

        return new Ret.cacheMissSpm(M, MU, MI, pi0, lE);
    }
}
