/**
 * @file The common corrected asymptotic expansion (LE-KT), computed on the cheaper side
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Pfqn_lekt {
    private Pfqn_lekt() {}

    /** Stirling remainder of a Gamma(a) direction, r(1) = 1 - log(2 pi)/2. */
    private static double r(double a) {
        return Gamma.logGamma(a) - (a - 0.5) * FastMath.log(a) + a - 0.5 * FastMath.log(2 * FastMath.PI);
    }

    /**
     * The side the estimator is computed on: "kt" when R <= M or a class self-loops
     * (one nonzero demand and no think time, which pfqn_kt extracts exactly), "le"
     * otherwise. The KT side is an R-dimensional convex solve and an R x R
     * determinant, the LE side an M-dimensional fixed point and an (M-1) x (M-1) one.
     */
    public static String route(Matrix L, Matrix N, Matrix Z) {
        int M = L.getNumRows(), R = L.getNumCols();
        boolean selfloop = false;
        if (R > 1) {
            for (int c = 0; c < R && !selfloop; c++) {
                int nnz = 0;
                for (int k = 0; k < M; k++) {
                    if (L.get(k, c) > GlobalConstants.Zero) nnz++;
                }
                double z = (Z == null || Z.isEmpty()) ? 0.0 : Z.get(0, c);
                if (nnz == 1 && z == 0.0) selfloop = true;
            }
        }
        return (R <= M || selfloop) ? "kt" : "le";
    }

    /**
     * The corrected logistic expansion (pfqn_ble) and the corrected Knessl-Tier
     * expansion (pfqn_bkt) are ONE estimator, evaluated in M-1 and in R dimensions:
     * with a think time their stationary points are one point in dual coordinates
     * (xi_r = N_r/(Z_r + v u'L_r), the class throughputs of the LE fixed point, and
     * v u_k = 1/(1-U_k), the M/M/1 factor of the KT saddle) and Sylvester's identity
     * exchanges the R x R Hessian determinant for the M x M one, after which every
     * 2 pi cancels. They agree to the accuracy of the two saddle-point solvers. Without
     * a think time the LE branch integrates the radius exactly as Gamma(N+M) while KT
     * Laplaces it, so the two differ by the constant (1-log(2 pi)/2) - r(N+M); the
     * common estimator is defined as the KT value, and the LE side here carries
     * M(1-log(2 pi)/2) - r(N+M) rather than pfqn_ble's (M-1)(1-log(2 pi)/2).
     * See _kb/03-api-layer.md.
     */
    public static Ret.pfqnNc pfqn_lekt(Matrix L, Matrix N, Matrix Z) {
        Matrix Zc = (Z == null || Z.isEmpty()) ? new Matrix(1, L.getNumCols()) : Z;
        if ("kt".equals(route(L, N, Zc))) {
            return Pfqn_bkt.pfqn_bkt(L, N, Zc);
        }
        Ret.pfqnNc base = Pfqn_ble.pfqn_ble(L, N, Zc);
        if (L.isEmpty() || N.isEmpty() || N.elementSum() == 0.0
                || L.elementSum() < GlobalConstants.CoarseTol) {
            return base; // pfqn_ble's degenerate branch: the delay term is exact
        }
        boolean noDelay = Zc.elementSum() < GlobalConstants.Zero;
        if (!noDelay) {
            return base;
        }
        // the Z = 0 branch of pfqn_ble counts M-1 directions; the common estimator
        // carries M kappa - r(N+M)
        double eta = N.elementSum() + L.getNumRows();
        double lGn = base.lG + Pfqn_ble.BLE_CORRECTION - r(eta);
        return new Ret.pfqnNc(Double.valueOf(FastMath.exp(lGn)), Double.valueOf(lGn));
    }

    public static Ret.pfqnNc pfqn_lekt(Matrix L, Matrix N) {
        return pfqn_lekt(L, N, new Matrix(1, L.getNumCols()));
    }
}
