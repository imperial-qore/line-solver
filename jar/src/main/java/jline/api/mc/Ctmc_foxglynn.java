/**
 * CTMC Transient Analysis via Fox-Glynn Uniformization
 *
 * Computes the transient probability distribution of CTMCs by uniformization with
 * Poisson weights and truncation points obtained by the Fox-Glynn algorithm, with
 * Jansen's correction to the right tail estimate. Unlike a direct evaluation of the
 * Poisson series, neither exp(-q*t) nor (q*t)^k/k! is ever formed, so the method is
 * free of overflow and underflow and needs no horizon splitting.
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

import jline.util.matrix.Matrix;

public final class Ctmc_foxglynn {
    private Ctmc_foxglynn() {}

    private static final double SQRT_2PI = FastMath.sqrt(2.0 * FastMath.PI);

    /**
     * Smallest rate for which the Fox-Glynn tail estimates, being asymptotic
     * normal approximations, are applicable.
     */
    private static final double FOX_GLYNN_MIN_LAMBDA = 25.0;

    /**
     * Poisson weights on the Fox-Glynn truncation window.
     */
    public static final class FoxGlynnWeights {
        public final int left;
        public final int right;
        public final double[] w;

        FoxGlynnWeights(int left, int right, double[] w) {
            this.left = left;
            this.right = right;
            this.w = w;
        }

        public double weight(int k) {
            if (k < left || k > right) {
                return 0.0;
            }
            return w[k - left];
        }
    }

    /**
     * Return the transient probability distribution of the CTMC by Fox-Glynn
     * uniformization.
     *
     * @param pi0 Initial state of the CTMC
     * @param Q   Infinitesimal generator of the CTMC
     * @param t   Transient analysis period boundary [0,t]
     * @return Transient probability vector at time t
     */
    public static Matrix ctmc_foxglynn(Matrix pi0, Matrix Q, double t) {
        return ctmc_foxglynn(pi0, Q, t, 1e-12, -1);
    }

    /**
     * Return the transient probability distribution of the CTMC by Fox-Glynn
     * uniformization.
     *
     * @param pi0     Initial state of the CTMC
     * @param Q       Infinitesimal generator of the CTMC
     * @param t       Transient analysis period boundary [0,t]
     * @param tol     Poisson tail-mass truncation tolerance
     * @param maxiter Maximum truncation depth; pass a nonpositive value to let
     *                the Fox-Glynn right truncation point size it
     * @return Transient probability vector at time t
     */
    public static Matrix ctmc_foxglynn(Matrix pi0, Matrix Q, double t, double tol, int maxiter) {
        int n = Q.getNumCols();
        double q = 0.0;
        for (int i = 0; i < n; i++) {
            q = FastMath.max(q, 1.1 * FastMath.abs(Q.get(i, i)));
        }
        double lambda = q * t;
        if (q <= 0.0 || lambda <= 0.0) {
            return new Matrix(pi0);
        }

        FoxGlynnWeights fg = ctmc_foxglynn_weights(lambda, tol, maxiter);
        Matrix Qs = Matrix.eye(n).add(1.0 / q, Q);
        Matrix pi = new Matrix(1, n);
        Matrix p = new Matrix(pi0);
        for (int k = 0; k <= fg.right; k++) {
            if (k >= fg.left) {
                pi.addEq(fg.w[k - fg.left], p);
            }
            if (k < fg.right) {
                p.multEq(Qs);
            }
        }
        return pi;
    }

    /**
     * Return the Fox-Glynn truncation window and normalized Poisson weights for
     * a Poisson(lambda) mixing distribution at tail-mass tolerance tol.
     *
     * @param lambda  Poisson rate, that is the uniformization constant times the horizon
     * @param tol     Poisson tail-mass truncation tolerance
     * @param maxiter Cap on the right truncation point; nonpositive to leave it uncapped
     * @return Truncation points and weights
     */
    public static FoxGlynnWeights ctmc_foxglynn_weights(double lambda, double tol, int maxiter) {
        return ctmc_foxglynn_weights(lambda, tol, maxiter, true);
    }

    /**
     * Return the Fox-Glynn truncation window and Poisson weights for a
     * Poisson(lambda) mixing distribution at tail-mass tolerance tol.
     *
     * <p>With normalize set the window is rescaled to sum to one, as Fox and
     * Glynn prescribe, so the truncated tails are redistributed over the
     * window. Cleared, the anchor is scaled by the true mode probability
     * through a log-gamma instead, so the returned values are the Poisson
     * probabilities themselves and 1 - sum(w) is the discarded tail rather
     * than being absorbed; {@link Ctmc_fau} needs that, its error being
     * reported as missing mass rather than as a bound.</p>
     *
     * @param lambda    Poisson rate, that is the uniformization constant times the horizon
     * @param tol       Poisson tail-mass truncation tolerance
     * @param maxiter   Cap on the right truncation point; nonpositive to leave it uncapped
     * @param normalize Rescale the window to sum to one, rather than returning
     *                  the true Poisson probabilities
     * @return Truncation points and weights
     */
    public static FoxGlynnWeights ctmc_foxglynn_weights(double lambda, double tol, int maxiter,
                                                        boolean normalize) {
        double eps = (tol > 0.0) ? tol : 1e-12;
        int left = leftTruncationPoint(lambda, eps);
        int right = rightTruncationPoint(lambda, eps);
        if (maxiter > 0 && right > maxiter) {
            right = maxiter;
            if (left > right) {
                left = right;
            }
        }
        return new FoxGlynnWeights(left, right, poissonWeights(lambda, left, right, normalize));
    }

    /**
     * Chernoff exponent of the Poisson(lambda) tail at k, that is
     * lambda*h(k/lambda) with h(u) = u*log(u) - u + 1. The bound
     * exp(-lambda*h(k/lambda)) dominates P{X &gt;= k} for k &gt; lambda and
     * P{X &lt;= k} for k &lt; lambda. It certifies the Fox-Glynn estimates
     * below, which are asymptotic and valid only above FOX_GLYNN_MIN_LAMBDA.
     */
    private static double chernoffExponent(double lambda, double k) {
        if (k <= 0.0) {
            return lambda;
        }
        return lambda - k + k * FastMath.log(k / lambda);
    }

    /**
     * Right truncation point R such that P{X &gt; R} &lt;= eps/2 for
     * X ~ Poisson(lambda).
     *
     * The starting guess is the Fox-Glynn (1988) right tail estimate. With
     * m = floor(lambda) and shift s(k) = k*sqrt(2*lambda) + 3/2, the tail
     * P{X &gt;= m + ceil(s(k))} is bounded by
     *
     *   a(lambda) * d(lambda,k) * exp(-k^2/2) / (k*sqrt(2*pi)),
     *   a(lambda) = (1 + 1/lambda) * exp(1/16) * sqrt(2).
     *
     * Jansen's correction (2011), applied here as the factor
     *
     *   d(lambda,k) = 1 / (1 - exp(-(2/9)*s(k))),
     *
     * repairs the original statement, which drops this factor for the
     * k-dependent shift and is therefore optimistic at moderate lambda; the
     * factor tends to one as lambda grows, so the corrected bound agrees with
     * Fox and Glynn's asymptotically. The guess is then tightened and, if
     * necessary, grown until the Chernoff bound above holds, so the returned R
     * is certified in any regime.
     */
    private static int rightTruncationPoint(double lambda, double eps) {
        double target = FastMath.log(2.0 / eps);
        int mode = (int) FastMath.floor(lambda);
        int r = mode;
        if (lambda >= FOX_GLYNN_MIN_LAMBDA) {
            double a = (1.0 + 1.0 / lambda) * FastMath.exp(1.0 / 16.0) * FastMath.sqrt(2.0);
            double spread = FastMath.sqrt(2.0 * lambda);
            for (int k = 1; k <= 64; k++) {
                double shift = k * spread + 1.5;
                double d = 1.0 / (1.0 - FastMath.exp(-(2.0 / 9.0) * shift));
                double bound = a * d * FastMath.exp(-0.5 * k * k) / (k * SQRT_2PI);
                if (bound <= 0.5 * eps) {
                    r = mode + (int) FastMath.ceil(shift);
                    break;
                }
            }
        }
        while (r > mode && chernoffExponent(lambda, r) >= target) {
            r--;
        }
        while (chernoffExponent(lambda, r + 1.0) < target) {
            r++;
        }
        return r;
    }

    /**
     * Left truncation point L such that P{X &lt; L} &lt;= eps/2 for
     * X ~ Poisson(lambda), zero when no truncation is admissible.
     *
     * The starting guess is the Fox-Glynn (1988) left tail estimate with
     * b(lambda) = (1 + 1/lambda)*exp(1/(8*lambda)) and shift k*sqrt(lambda) + 3/2
     * below the mode. Unlike the right tail this one needs no Jansen factor,
     * the left tail of a Poisson being lighter than its normal approximation.
     * It is certified against the Chernoff bound in the same way.
     */
    private static int leftTruncationPoint(double lambda, double eps) {
        double target = FastMath.log(2.0 / eps);
        int mode = (int) FastMath.floor(lambda);
        if (chernoffExponent(lambda, 0.0) < target) {
            return 0;
        }
        int l = 0;
        if (lambda >= FOX_GLYNN_MIN_LAMBDA) {
            double b = (1.0 + 1.0 / lambda) * FastMath.exp(1.0 / (8.0 * lambda));
            double spread = FastMath.sqrt(lambda);
            for (int k = 1; k <= 64; k++) {
                double bound = b * FastMath.exp(-0.5 * k * k) / (k * SQRT_2PI);
                if (bound <= 0.5 * eps) {
                    l = mode - (int) FastMath.floor(k * spread + 1.5);
                    break;
                }
            }
            if (l < 0) {
                l = 0;
            }
        }
        while (l > 0 && chernoffExponent(lambda, l - 1.0) < target) {
            l--;
        }
        while (l < mode && chernoffExponent(lambda, (double) l) >= target) {
            l++;
        }
        return l;
    }

    /**
     * Poisson(lambda) weights on [left, right], normalized to sum to one.
     *
     * Following Fox-Glynn, the weights are built by the two-sided recursion
     * w[k-1] = w[k]*k/lambda and w[k+1] = w[k]*lambda/(k+1) anchored at the
     * mode, so neither exp(-lambda) nor lambda^k/k! is ever evaluated and the
     * overflow and underflow that limit the direct series cannot occur.
     * Anchoring at w[mode] = 1 keeps the extreme weights near eps, far above
     * the denormal threshold, making Fox and Glynn's rescaling of the mode
     * weight unnecessary here. The normalizing sum is accumulated from the two
     * tails inwards, in increasing order of magnitude, with compensated
     * summation.
     */
    private static double[] poissonWeights(double lambda, int left, int right, boolean normalize) {
        int len = right - left + 1;
        double[] w = new double[len];
        int mode = (int) FastMath.floor(lambda);
        if (mode < left) {
            mode = left;
        }
        if (mode > right) {
            mode = right;
        }
        w[mode - left] = 1.0;
        for (int k = mode; k > left; k--) {
            w[k - 1 - left] = w[k - left] * k / lambda;
        }
        for (int k = mode; k < right; k++) {
            w[k + 1 - left] = w[k - left] * lambda / (k + 1);
        }

        if (!normalize) {
            double logMode = -lambda + mode * FastMath.log(lambda) - Gamma.logGamma(mode + 1.0);
            double scale = FastMath.exp(logMode);
            for (int i = 0; i < len; i++) {
                w[i] *= scale;
            }
            return w;
        }

        double total = 0.0;
        double comp = 0.0;
        int lo = 0;
        int hi = len - 1;
        while (lo <= hi) {
            double next;
            if (w[lo] <= w[hi]) {
                next = w[lo++];
            } else {
                next = w[hi--];
            }
            double y = next - comp;
            double s = total + y;
            comp = (s - total) - y;
            total = s;
        }
        for (int i = 0; i < len; i++) {
            w[i] /= total;
        }
        return w;
    }
}
