/**
 * CTMC Transient Analysis via Uniformization
 *
 * Computes the transient probability distribution of CTMCs using the uniformization
 * method, which transforms the continuous-time problem into a weighted sum of DTMC
 * powers. Provides numerically stable computation of time-dependent state probabilities.
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import org.apache.commons.math3.util.FastMath;

import jline.util.matrix.Matrix;

public final class Ctmc_uniformization {
    private Ctmc_uniformization() {}

    /**
     * Return the transient probability distribution of the CTMC via the
     * uniformization method. The Poisson series truncation depth is chosen
     * adaptively from the uniformization constant so that the tail mass is
     * below tolerance for any horizon t.
     *
     * @param pi0 Initial state the CTMC
     * @param Q   Infinitesimal generator of the CTMC
     * @param t   Transient analysis period boundary [0,t]
     * @return Transient probability vector at time t
     */
    public static Matrix ctmc_uniformization(Matrix pi0, Matrix Q, double t) {
        return ctmc_uniformization(pi0, Q, t, 1e-12, -1);
    }

    /**
     * Return the transient probability distribution of the CTMC via the
     * uniformization method.
     *
     * @param pi0     Initial state the CTMC
     * @param Q       Infinitesimal generator of the CTMC
     * @param t       Transient analysis period boundary [0,t]
     * @param tol     Poisson tail-mass truncation tolerance
     * @param maxiter Maximum truncation depth; pass a nonpositive value to size
     *                it adaptively as max(100, q*t + 10*sqrt(q*t) + 20)
     * @return Transient probability vector at time t
     */
    /**
     * Largest q*t per uniformization segment: exp(-q*t) must stay above the
     * double underflow threshold (exp(-745) == 0) for the Poisson recursion.
     */
    private static final double MAX_QT_PER_SEGMENT = 500.0;

    public static Matrix ctmc_uniformization(Matrix pi0, Matrix Q, double t, double tol, int maxiter) {
        double q = 0.0;
        int n = Q.getNumCols();
        for (int i = 0; i < n; i++) {
            q = FastMath.max(q, 1.1 * FastMath.abs(Q.get(i, i)));
        }
        if (q * t > MAX_QT_PER_SEGMENT) {
            // Split the horizon so exp(-q*t) never underflows within a segment;
            // exp(Q*t) = (exp(Q*t/nSeg))^nSeg applied to the row vector
            int nSeg = (int) FastMath.ceil(q * t / MAX_QT_PER_SEGMENT);
            double tSeg = t / nSeg;
            Matrix pi = pi0;
            for (int seg = 0; seg < nSeg; seg++) {
                pi = ctmc_uniformization(pi, Q, tSeg, tol, maxiter);
            }
            return pi;
        }
        if (maxiter <= 0) {
            // The Poisson(q*t) mass concentrates around q*t with spread
            // O(sqrt(q*t)); a fixed cap silently truncates for large horizons
            maxiter = (int) FastMath.max(100, FastMath.ceil(q * t + 10 * FastMath.sqrt(q * t) + 20));
        }
        Matrix Qs = Matrix.eye(n);
        Qs = Qs.add(1.0 / q, Q);
        int k = 0;
        double s = 1.0;
        double r = 1.0;
        int iter = 0;
        int kmax = 1;
        while (iter < maxiter) {
            iter++;
            k++;
            r = r * (q * t) / k;
            s = s + r;
            if (1 - FastMath.exp(-q * t) * s <= tol) {
                kmax = k;
                break;
            }
            // see _kb/03-api-layer.md for rationale
            kmax = k;
        }

        Matrix pi = new Matrix(1, n);
        pi0.scaleEq(Math.exp(-q * t), pi);
        Matrix P = new Matrix(pi0);
        double ri = FastMath.exp(-q * t);
        for (int j = 0; j < kmax; j++) {
            P.multEq(Qs);
            ri = ri * (q * t / (j + 1));
            pi = pi.add(ri, P);
        }
        return pi;
    }
}
