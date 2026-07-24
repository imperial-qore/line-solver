/**
 * @file Multi-class repairman model sampling method for normalizing constants
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import java.util.Random;

import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class Pfqn_mmsample2 {
    private Pfqn_mmsample2() {}

    public static Ret.pfqnNc pfqn_mmsample2(Matrix L, Matrix N, Matrix Z, int samples) {
        Matrix L_local = L.copy();
        Matrix Z_local = Z.copy();
        int R = N.getNumElements();
        double scaleFactor = 1e-7 + FastMath.min(L_local.elementMin(), Z_local.elementMin());
        L_local.scaleEq(1.0 / scaleFactor);
        Z_local.scaleEq(1.0 / scaleFactor);

        double c = 0.5;
        int numSamples1 = (int) FastMath.ceil(c * samples);
        int numSamples2 = (int) FastMath.ceil(samples * (1 - c));
        double[] v = new double[numSamples1 + numSamples2];
        double[] du = new double[numSamples1 + numSamples2];

        Random rand = new Random();
        for (int i = 0; i < numSamples1; i++) {
            v[i] = rand.nextDouble();
        }
        double[] lv = Maths.logSpace(0.0, 5.0, numSamples2);
        if (numSamples2 >= 0) System.arraycopy(lv, 0, v, numSamples1, numSamples2);

        // see _kb/03-api-layer.md for rationale
        java.util.Arrays.sort(v);
        for (int i = 0; i < du.length; i++) {
            du[i] = (i == 0) ? v[0] : v[i] - v[i - 1];
        }

        // see _kb/03-api-layer.md for rationale
        double[] lterms = new double[v.length];
        double[] Lr = new double[R];
        double[] Zr = new double[R];
        double[] Nr = new double[R];
        for (int r = 0; r < R; r++) {
            Lr[r] = L_local.get(0, r);
            Zr[r] = Z_local.get(0, r);
            Nr[r] = N.get(r);
        }
        double lmax = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < v.length; i++) {
            double sum = 0.0;
            for (int r = 0; r < R; r++) {
                // McKenna-Mitra integrand is (Z_r + L_r*u), NOT (Z_r + L_r)*u:
                // the latter is a different function that merely agrees at u=1.
                sum += FastMath.log(Zr[r] + Lr[r] * v[i]) * Nr[r];
            }
            // Log-domain quadrature: the panel weight enters as log(du), not as
            // du added to a log.
            double lterm = FastMath.log(du[i]) - v[i] + sum;
            lterms[i] = lterm;
            if (lterm > lmax) {
                lmax = lterm;
            }
        }

        // see _kb/03-api-layer.md for rationale
        double acc = 0.0;
        for (int i = 0; i < v.length; i++) {
            acc += FastMath.exp(lterms[i] - lmax);
        }
        double lG = lmax + FastMath.log(acc) - N.factln().elementSum();
        lG = lG + N.elementSum() * FastMath.log(scaleFactor);

        return new Ret.pfqnNc(Double.valueOf(Math.exp(lG)), Double.valueOf(lG));
    }
}
