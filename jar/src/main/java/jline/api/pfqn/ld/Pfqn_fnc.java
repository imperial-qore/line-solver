/**
 * @file Functional server scaling factor computation for load-dependent systems
 *
 * Computes scaling factors for load-dependent functional servers in product-form queueing networks.
 * Handles the mathematical transformation of load-dependent service rates into functional scaling
 * parameters, supporting both automatic parameter selection and user-specified scaling constants.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.util.Utils;
import jline.util.matrix.Matrix;

public final class Pfqn_fnc {
    private Pfqn_fnc() {}

    /**
     * Compute scaling factor of a load-dependent functional server use to calculate the mean
     */
    public static Ret.pfqnFnc pfqn_fnc(Matrix alpha) {
        int M = alpha.getNumRows();
        if (alpha.getNumCols() == 0) {
            // see _kb/03-api-layer.md for rationale
            Matrix cEmpty = new Matrix(1, M);
            cEmpty.zero();
            return new Ret.pfqnFnc(new Matrix(M, 0), cEmpty);
        }
        Matrix c = new Matrix(1, M);
        c.zero();
        Matrix mu = pfqn_fnc(alpha, c).mu;
        if (!mu.isFinite()) {
            c = Matrix.ones(1, M);
            c.scaleEq(-0.5);
            mu = pfqn_fnc(alpha, c).mu;
        }
        double dt = 0.0;
        while (!mu.isFinite()) {
            dt += 0.05;
            double c_scalar = -0.5 + dt;
            // see _kb/03-api-layer.md for rationale
            c = Matrix.ones(1, M);
            c.scaleEq(c_scalar);
            mu = pfqn_fnc(alpha, c).mu;
            if (c_scalar >= 2) {
                break;
            }
        }
        return new Ret.pfqnFnc(mu, c);
    }

    /**
     * Compute scaling factor of a load-dependent functional server use to calculate the mean instantiated
     * with scaling constant c.
     */
    public static Ret.pfqnFnc pfqn_fnc(Matrix alpha, Matrix c) {
        int M = alpha.getNumRows();
        int N = alpha.getNumCols();
        if (N == 0) {
            // No rate columns (see one-argument overload): nothing to build.
            Matrix cEmpty = new Matrix(1, M);
            cEmpty.zero();
            return new Ret.pfqnFnc(new Matrix(M, 0), cEmpty);
        }
        Matrix mu = new Matrix(M, N);
        mu.zero();
        for (int i = 0; i < M; i++) {
            mu.set(i, 0, alpha.get(i, 0) / (1 + c.get(i)));
            Matrix alphanum = new Matrix(N, N);
            alphanum.zero();
            Matrix alphaden = alphanum.copy();
            for (int n = 1; n < N; n++) {
                alphanum.set(n, 0, alpha.get(i, n));
                alphaden.set(n, 0, alpha.get(i, n - 1));
                for (int k = 1; k < n; k++) {
                    alphanum.set(n, k, alphanum.get(n, k - 1) * alpha.get(i, n - k));
                    alphaden.set(n, k, alphaden.get(n, k - 1) * alpha.get(i, n - k - 1));
                }
            }
            for (int n = 1; n < N; n++) {
                double rho = 0.0;
                double muden = 1.0;
                for (int k = 0; k < n; k++) {
                    muden *= mu.get(i, k);
                    rho += (alphanum.get(n, k) - alphaden.get(n, k)) / muden;
                }
                mu.set(i, n, (alphanum.get(n, n - 1) * alpha.get(i, 0) / muden) / (1 - rho));
            }
        }
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                if (Double.isNaN(mu.get(i, j)) || FastMath.abs(mu.get(i, j)) > 1e15) {
                    mu.set(i, j, GlobalConstants.Inf);
                }
            }
        }
        for (int i = 0; i < M; i++) {
            if (Matrix.extractRows(mu, i, i + 1, null).isFinite()) {
                continue;
            }
            boolean replaceWithInf = false;
            for (int j = 0; j < N; j++) {
                if (replaceWithInf) {
                    mu.set(i, j, GlobalConstants.Inf);
                } else if (Utils.isInf(mu.get(i, j))) {
                    replaceWithInf = true;
                }
            }
        }
        return new Ret.pfqnFnc(mu, c);
    }
}
