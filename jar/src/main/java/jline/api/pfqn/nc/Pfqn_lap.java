/**
 * @file Laplace approximation for normalizing constant
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import org.apache.commons.math3.util.FastMath;

import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class Pfqn_lap {
    private Pfqn_lap() {}

    /**
     * Compute the Laplace approximation for the log normalizing constant.
     */
    public static double pfqn_lap(Matrix L, Matrix N, Matrix Z) {
        if (L.getNumRows() > 1 && L.getNumCols() > 1
                || L.getNumElements() != N.getNumElements()) {
            // Same contract as MATLAB pfqn_lap: per-class vectors for a single
            // queueing station (repairman models)
            throw new IllegalArgumentException(
                    "pfqn_lap expects per-class vectors for a single queueing station "
                            + "(repairman models): L, N, Z must be 1xR.");
        }
        double Ntot = N.elementSum();
        int R = N.getNumCols();

        double u0 = findZero(L, N, Z, Ntot, R);

        if (u0 < 0) {
            return Double.NaN;
        }

        if (!Double.isFinite(u0)) {
            double initSign = Math.signum(evalF(0.001, L, N, Z, Ntot, R));
            u0 = Double.NaN;
            double x = 1e-4;
            while (x <= 10.0) {
                double fx = evalF(x, L, N, Z, Ntot, R);
                if (Math.signum(fx) != initSign) {
                    u0 = x;
                    break;
                }
                x += 1e-4;
            }
            if (!Double.isFinite(u0) || u0 < 0) {
                return Double.NaN;
            }
        }

        if (u0 < 0) {
            return Double.NaN;
        }

        double logI = FastMath.log(Ntot);
        for (int r = 0; r < R; r++) {
            logI -= Maths.factln(N.get(r));
        }
        logI -= Ntot * u0;
        for (int r = 0; r < R; r++) {
            logI += N.get(r) * FastMath.log(Z.get(r) + L.get(r) * u0 * Ntot);
        }
        logI += 0.5 * FastMath.log(2 * Math.PI);

        double f2 = 0.0;
        for (int r = 0; r < R; r++) {
            if (L.get(r) > 0) {
                double term = Z.get(r) / (Ntot * L.get(r)) + u0;
                f2 += (N.get(r) / Ntot) / (term * term);
            }
        }
        logI -= 0.5 * FastMath.log(f2);
        logI -= 0.5 * FastMath.log(Ntot);

        return logI;
    }

    private static double evalF(double x, Matrix L, Matrix N, Matrix Z, double Ntot, int R) {
        double sum = 0.0;
        for (int r = 0; r < R; r++) {
            double denom = Z.get(r) + Ntot * L.get(r) * x;
            if (denom > 0) {
                sum += N.get(r) * L.get(r) / denom;
            }
        }
        return 1.0 - sum;
    }

    private static double findZero(Matrix L, Matrix N, Matrix Z, double Ntot, int R) {
        double lo = 1e-10;
        double hi = 100.0;

        double flo = evalF(lo, L, N, Z, Ntot, R);
        double fhi = evalF(hi, L, N, Z, Ntot, R);

        if (flo * fhi > 0) {
            if (flo > 0) {
                lo = 1e-15;
                if (evalF(lo, L, N, Z, Ntot, R) > 0) {
                    return -1.0;
                }
            } else {
                hi = 10000.0;
                fhi = evalF(hi, L, N, Z, Ntot, R);
                if (fhi < 0) {
                    return Double.POSITIVE_INFINITY;
                }
            }
        }

        int maxIter = 200;
        double tol = 1e-14;
        for (int iter = 0; iter < maxIter; iter++) {
            double mid = (lo + hi) / 2;
            double fmid = evalF(mid, L, N, Z, Ntot, R);
            if (FastMath.abs(fmid) < tol || (hi - lo) / 2 < tol) {
                return mid;
            }
            if (fmid * evalF(lo, L, N, Z, Ntot, R) < 0) {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        return (lo + hi) / 2;
    }
}
