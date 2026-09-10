/**
 * @file Logistic expansion method for normalizing constant computation
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class Pfqn_le {
    private Pfqn_le() {}

    /**
     * Logistic expansion method to compute the normalizing constant.
     */
    public static Ret.pfqnNc pfqn_le(Matrix L, Matrix N, Matrix Z) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        double lGn;
        double Gn;

        if (L.isEmpty() || N.isEmpty() || N.elementSum() == 0.0
                || L.elementSum() < GlobalConstants.CoarseTol) {
            // Z can be empty here, and an empty class contributes 0, not 0*log(0).
            lGn = -Matrix.factln(N).elementSum();
            for (int r = 0; r < N.length(); r++) {
                if (N.get(r) > 0) {
                    lGn += N.get(r) * FastMath.log(Z.isEmpty() ? 0.0 : Z.sumCols(r));
                }
            }
            Gn = FastMath.exp(lGn);
        // Branch on an all-zero Z, not merely an absent one: with no think time the
        // radial integral is exact as gamma(N+M) and only the simplex is Laplaced.
        // Pfqn_ble counts Laplaced directions off this same predicate, and pfqn_nc
        // always passes a Z vector, so testing isEmpty alone sent every delay-free
        // model down the Z>0 branch while Pfqn_ble still corrected it as if Z=0.
        } else if (Z.isEmpty() || Z.elementSum() < GlobalConstants.Zero) {
            Ret.pfqnLeFpi ret = Pfqn_le_fpi.pfqn_le_fpi(L, N);
            Matrix umax = ret.u;
            Matrix A = Pfqn_le_hessian.pfqn_le_hessian(L, N, umax.transpose());
            double S = 0.0;
            for (int r = 0; r < R; r++) {
                Matrix L_col_r = new Matrix(L.getNumRows(), 1);
                Matrix.extract(L, 0, L.getNumRows(), r, r + 1, L_col_r, 0, 0);
                S += N.get(r) * FastMath.log(umax.transpose().mult(L_col_r).get(0));
            }

            Matrix tmp = new Matrix(1, N.length() + 1);
            Matrix.extract(N, 0, 1, 0, N.length(), tmp, 0, 0);
            tmp.set(N.length(), (double) (M - 1));
            Matrix log_umax = umax.copy();
            for (int i = 0; i < log_umax.length(); i++) {
                log_umax.set(i, FastMath.log(log_umax.get(i)));
            }
            lGn = (Maths.multinomialln(tmp) + Maths.factln(M - 1)
                    // Cas17 eq.(34) as published; Pfqn_ble adds the eps->0 bias correction.
                    + (M - 1) * FastMath.log(Math.sqrt(2 * FastMath.PI))
                    - FastMath.log(Math.sqrt(A.det()))) + log_umax.elementSum() + S;
            Gn = FastMath.exp(lGn);
        } else {
            Ret.pfqnLeFpiZ ret = Pfqn_le_fpiZ.pfqn_le_fpiZ(L, N, Z);
            Matrix umax = ret.u;
            double vmax = ret.v;
            Matrix A = Pfqn_le_hessianZ.pfqn_le_hessianZ(L, N, Z, umax.transpose(), vmax);
            double S = 0.0;
            for (int r = 0; r < R; r++) {
                Matrix L_col_r = new Matrix(L.getNumRows(), 1);
                Matrix.extract(L, 0, L.getNumRows(), r, r + 1, L_col_r, 0, 0);
                S += N.get(r) * FastMath.log(Z.get(r) + vmax * umax.transpose().mult(L_col_r).get(0));
            }
            Matrix log_umax = umax.copy();
            for (int i = 0; i < log_umax.length(); i++) {
                log_umax.set(i, FastMath.log(log_umax.get(i)));
            }
            lGn = (-Matrix.factln(N).elementSum() - vmax + M * FastMath.log(vmax)
                    + M * FastMath.log(Math.sqrt(2 * FastMath.PI))
                    - FastMath.log(Math.sqrt(A.det()))) + log_umax.elementSum() + S;
            Gn = FastMath.exp(lGn);
        }
        return new Ret.pfqnNc(Double.valueOf(Gn), Double.valueOf(lGn));
    }

    public static Ret.pfqnNc pfqn_le(Matrix L, Matrix N) {
        return pfqn_le(L, N, new Matrix(1, L.getNumCols()));
    }
}
