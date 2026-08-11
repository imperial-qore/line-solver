/**
 * CTMC time-averaged transient distribution via uniformization
 *
 * Companion of {@link Ctmc_uniformization}, which returns only the endpoint
 * pi0*exp(Q*T); this routine additionally returns the time average
 *
 *   piTimeAvg = pi0 * (1/T) * \int_0^T exp(Q*tau) d(tau)
 *
 * as well as the endpoint piExit = pi0*exp(Q*T), both obtained from the same
 * Jensen uniformization series without forming any dense matrix exponential.
 * Used by the SolverENV state-vector analyzer (deterministic-sojourn option).
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import org.apache.commons.math3.util.FastMath;

import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class Ctmc_timeaverage {
    private Ctmc_timeaverage() {}

    /**
     * Time-averaged transient distribution of a CTMC over [0,t].
     *
     * Uniformization: with q = 1.1*max|diag(Q)| and P = I + Q/q (row-stochastic),
     *   pi0*exp(Q*t)             = sum_j w_j(qt) * (pi0*P^j)
     *   pi0*\int_0^t exp(Q*tau)  = (1/q) * sum_j (1 - W_j(qt)) * (pi0*P^j)
     * where w_j and W_j are the Poisson(qt) PMF and CDF. The time average
     * divides the integral by t (equivalently the integral sum by q*t).
     *
     * @param pi0 initial (arbitrary) distribution, row vector
     * @param Q   infinitesimal generator
     * @param t   horizon
     * @return Pair(piTimeAvg, piExit), both row vectors
     */
    public static Pair<Matrix, Matrix> ctmc_timeaverage(Matrix pi0, Matrix Q, double t) {
        return ctmc_timeaverage(pi0, Q, t, 1e-12, -1);
    }

    /**
     * Largest q*t per uniformization segment: exp(-q*t) must stay above the
     * double underflow threshold (exp(-745) == 0) for the Poisson recursion.
     */
    private static final double MAX_QT_PER_SEGMENT = 500.0;

    public static Pair<Matrix, Matrix> ctmc_timeaverage(Matrix pi0, Matrix Q, double t, double tol, int maxiter) {
        int n = Q.getNumCols();
        double q = 0.0;
        for (int i = 0; i < n; i++) {
            q = FastMath.max(q, 1.1 * FastMath.abs(Q.get(i, i)));
        }
        if (q * t > MAX_QT_PER_SEGMENT) {
            // see _kb/03-api-layer.md for rationale
            int nSeg = (int) FastMath.ceil(q * t / MAX_QT_PER_SEGMENT);
            double tSeg = t / nSeg;
            Matrix piCur = pi0;
            Matrix integral = new Matrix(1, n);
            for (int seg = 0; seg < nSeg; seg++) {
                Pair<Matrix, Matrix> segResult = ctmc_timeaverage(piCur, Q, tSeg, tol, maxiter);
                integral = integral.add(tSeg, segResult.getLeft());
                piCur = segResult.getRight();
            }
            integral.scaleEq(1.0 / t);
            return new Pair<Matrix, Matrix>(integral, piCur);
        }
        if (maxiter <= 0) {
            // The Poisson(q*t) mass concentrates around q*t with spread
            // O(sqrt(q*t)); a fixed cap silently truncates for large horizons
            maxiter = (int) FastMath.max(100, FastMath.ceil(q * t + 10 * FastMath.sqrt(q * t) + 20));
        }
        Matrix Qs = Matrix.eye(n);
        Qs = Qs.add(1.0 / q, Q);
        double qt = q * t;

        // Number of Poisson terms needed (right-tail below tol).
        int k = 0;
        double s = 1.0;
        double r = 1.0;
        int iter = 0;
        int kmax = 1;
        while (iter < maxiter) {
            iter++;
            k++;
            r = r * qt / k;
            s = s + r;
            if (1 - FastMath.exp(-qt) * s <= tol) {
                kmax = k;
                break;
            }
            // Best-effort truncation depth if the loop exhausts maxiter
            kmax = k;
        }

        // Accumulate endpoint and integral over the shared pi0*P^j sequence.
        double w = FastMath.exp(-qt);   // Poisson PMF  w_0
        double W = w;                   // Poisson CDF  W_0
        Matrix P = new Matrix(pi0);     // pi0*P^0
        Matrix piExit = new Matrix(1, n);
        P.scaleEq(w, piExit);
        Matrix piIntSum = new Matrix(1, n);
        P.scaleEq(FastMath.max(1 - W, 0.0), piIntSum);
        for (int j = 1; j <= kmax; j++) {
            P.multEq(Qs);               // pi0*P^j
            w = w * qt / j;             // w_j
            W = W + w;                  // W_j
            piExit = piExit.add(w, P);
            piIntSum = piIntSum.add(FastMath.max(1 - W, 0.0), P);
        }
        piIntSum.scaleEq(1.0 / qt);     // (1/q)*sum / t  =  sum/(q*t)
        return new Pair<Matrix, Matrix>(piIntSum, piExit);
    }
}
