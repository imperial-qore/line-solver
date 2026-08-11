/**
 * @file ForkTail black-box tail-latency approximation for fork-join requests
 *
 * @since LINE 3.0
 */
package jline.api.fj;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.BrentSolver;

/**
 * ForkTail black-box tail-latency approximation for fork-join requests.
 *
 * <p>Approximates the p-th percentile of the response time of a request that
 * forks into K parallel tasks and joins on the last of them, from the mean and
 * variance of the per-branch task response times alone. Each branch is treated
 * as a black box: its task response time is fitted by a generalized exponential
 * law
 *
 * <pre>
 *   F_T(x) = (1 - exp(-x/beta))^alpha
 *   E[T]   = beta*(psi(alpha+1) - psi(1))
 *   V[T]   = beta^2*(psi'(1) - psi'(alpha+1))
 * </pre>
 *
 * and the request response time is the maximum over the branches, taken as the
 * product of the branch CDFs (exact only for independent branches):
 *
 * <pre>
 *   F_X(x) = prod_i (1 - exp(-x/beta_i))^alpha_i,   x_p = F_X^{-1}(p)
 *   homogeneous: x_p = -beta*log(1 - p^(1/(K*alpha)))
 *   random fanout: F_X(x) = sum_i P_i * (1 - exp(-x/beta))^(K_i*alpha)
 * </pre>
 *
 * <p>The approximation rests on the central limit theorem for G/G/m queues in
 * heavy traffic, so it is a HIGH-LOAD result: the reference reports errors
 * within 20% and 15% at 80% and 90% utilization respectively, and makes no
 * claim at low load, where the tail is dominated by the service law rather than
 * by queueing. Prefer the FJ_codes route (see {@link FJValidation#isHomogeneous})
 * when the model is in the homogeneous MAP/PH/1 class, which is more accurate
 * there; ForkTail covers the heterogeneous branches and mixed service laws that
 * route rejects.
 *
 * <p>Port of matlab/src/api/fj/fj_tail_forktail.m and
 * matlab/src/api/fj/fj_mg1_respt_moments.m.
 *
 * <p>Reference: M. Nguyen, S. Alesawi, N. Li, H. Che, H. Jiang, "ForkTail: A
 * Black-Box Fork-Join Tail Latency Prediction Model for User-Facing Datacenter
 * Workloads", ACM HPDC 2018, pp. 206-217.
 */
public final class FJ_tail_forktail {
    private FJ_tail_forktail() {}

    private static final double FINE_TOL = 1e-8;
    private static final double COARSE_TOL = 1e-3;

    /** Fitted generalized-exponential parameters of one branch. */
    public static class GEFit {
        public final double alpha;
        public final double beta;

        public GEFit(double alpha, double beta) {
            this.alpha = alpha;
            this.beta = beta;
        }
    }

    /** Mean and variance of a branch response time. */
    public static class ResptMoments {
        public final double mean;
        public final double variance;

        public ResptMoments(double mean, double variance) {
            this.mean = mean;
            this.variance = variance;
        }
    }

    /**
     * Mean and variance of the M/G/1 FCFS response time, the white-box inputs of
     * ForkTail, from the first three moments of the service time.
     *
     * @param lambda arrival rate at the branch
     * @param ES first moment of the service time
     * @param ES2 second moment of the service time
     * @param ES3 third moment of the service time
     * @return the response time mean and variance
     */
    public static ResptMoments fj_mg1_respt_moments(double lambda, double ES, double ES2, double ES3) {
        double rho = lambda * ES;
        if (rho >= 1.0) {
            throw new IllegalArgumentException(
                "The branch is unstable (rho = " + rho + " >= 1); the response time moments do not exist.");
        }
        if (Double.isInfinite(ES3) || Double.isNaN(ES3)) {
            throw new IllegalArgumentException(
                "The service law has no finite third moment, so the ForkTail response time "
                + "variance is undefined (a Pareto branch with shape <= 3, for instance). "
                + "Use a service law with three finite moments.");
        }
        double scvS = (ES2 - ES * ES) / (ES * ES);
        double mean = ES * (1.0 + rho / (1.0 - rho) * (1.0 + scvS) / 2.0);
        double EW = lambda * ES2 / (2.0 * (1.0 - rho));
        double variance = EW * EW + lambda * ES3 / (3.0 * (1.0 - rho)) + ES2 - ES * ES;
        return new ResptMoments(mean, variance);
    }

    /**
     * Tail latency of a request forked over K identical branches.
     *
     * @param ET mean task response time of a branch
     * @param VT variance of the task response time of a branch
     * @param K number of branches
     * @param p percentile, a fraction in (0,1) or a percentage in (0,100)
     * @return the predicted p-th percentile of the request response time
     */
    public static double fj_tail_forktail(double ET, double VT, int K, double p) {
        double pp = normalizePercentile(p);
        GEFit fit = geFit(ET, VT);
        return -fit.beta * Math.log(1.0 - Math.pow(pp, 1.0 / (K * fit.alpha)));
    }

    /**
     * Tail latency of a request forked over heterogeneous branches, one entry
     * per branch.
     *
     * @param ET mean task response time per branch
     * @param VT variance of the task response time per branch
     * @param p percentile, a fraction in (0,1) or a percentage in (0,100)
     * @return the predicted p-th percentile of the request response time
     */
    public static double fj_tail_forktail(double[] ET, double[] VT, double p) {
        if (ET.length != VT.length) {
            throw new IllegalArgumentException("ET and VT must have the same number of entries.");
        }
        double pp = normalizePercentile(p);
        final int n = ET.length;
        if (n == 1) {
            return fj_tail_forktail(ET[0], VT[0], 1, pp);
        }
        final double[] alpha = new double[n];
        final double[] beta = new double[n];
        double xlo = 0.0;
        for (int i = 0; i < n; i++) {
            GEFit fit = geFit(ET[i], VT[i]);
            alpha[i] = fit.alpha;
            beta[i] = fit.beta;
            // F_X is bounded above by any single branch CDF, so the request
            // percentile is at least the largest branch percentile
            xlo = Math.max(xlo, -beta[i] * Math.log(1.0 - Math.pow(pp, 1.0 / alpha[i])));
        }
        final double logp = Math.log(pp);
        UnivariateFunction residual = new UnivariateFunction() {
            @Override
            public double value(double x) {
                double acc = 0.0;
                for (int i = 0; i < n; i++) {
                    acc += alpha[i] * Math.log1p(-Math.exp(-x / beta[i]));
                }
                return acc - logp;
            }
        };
        double xhi = xlo;
        while (residual.value(xhi) < 0.0) {
            xhi *= 2.0;
            if (Double.isInfinite(xhi)) {
                throw new RuntimeException("Could not bracket the ForkTail percentile.");
            }
        }
        double lo = xlo;
        while (residual.value(lo) > 0.0 && lo > Double.MIN_NORMAL) {
            lo /= 2.0;
        }
        return new BrentSolver(1e-12, 1e-12).solve(200, residual, lo, xhi);
    }

    /**
     * Tail latency when the fanout itself is random: a request spawns K[i] tasks
     * with probability P[i], over identical branches.
     *
     * @param ET mean task response time of a branch
     * @param VT variance of the task response time of a branch
     * @param K distinct fanouts
     * @param P probabilities of those fanouts, summing to one
     * @param p percentile, a fraction in (0,1) or a percentage in (0,100)
     * @return the predicted p-th percentile of the request response time
     */
    public static double fj_tail_forktail(double ET, double VT, final int[] K, final double[] P, double p) {
        if (K.length != P.length) {
            throw new IllegalArgumentException("A vector of fanouts K needs a probability vector P of the same length.");
        }
        double psum = 0.0;
        for (int i = 0; i < P.length; i++) {
            if (P[i] < 0.0) {
                throw new IllegalArgumentException("The fanout probabilities P must be non-negative.");
            }
            psum += P[i];
        }
        if (Math.abs(psum - 1.0) > COARSE_TOL) {
            throw new IllegalArgumentException("The fanout probabilities P must sum to one.");
        }
        final double pp = normalizePercentile(p);
        GEFit fit = geFit(ET, VT);
        final double alpha = fit.alpha;
        final double beta = fit.beta;
        int kmin = Integer.MAX_VALUE;
        int kmax = Integer.MIN_VALUE;
        for (int k : K) {
            kmin = Math.min(kmin, k);
            kmax = Math.max(kmax, k);
        }
        UnivariateFunction residual = new UnivariateFunction() {
            @Override
            public double value(double x) {
                double acc = 0.0;
                for (int i = 0; i < K.length; i++) {
                    acc += P[i] * Math.pow(1.0 - Math.exp(-x / beta), K[i] * alpha);
                }
                return acc - pp;
            }
        };
        double xlo = -beta * Math.log(1.0 - Math.pow(pp, 1.0 / (kmin * alpha)));
        double xhi = -beta * Math.log(1.0 - Math.pow(pp, 1.0 / (kmax * alpha)));
        if (xlo == xhi) {
            return xlo;
        }
        return new BrentSolver(1e-12, 1e-12).solve(200, residual, Math.min(xlo, xhi), Math.max(xlo, xhi));
    }

    /**
     * Matches a generalized exponential law on a mean and a variance. The
     * squared coefficient of variation depends on the shape alone and decreases
     * monotonically in it, so the shape is recovered by a scalar root-find on a
     * logarithmic scale and the scale then follows in closed form. SCV = 1 is
     * the exponential case alpha = 1, kept exact.
     *
     * @param ET mean of the branch response time
     * @param VT variance of the branch response time
     * @return the fitted shape and scale
     */
    public static GEFit geFit(double ET, double VT) {
        if (ET <= 0.0 || VT <= 0.0) {
            throw new IllegalArgumentException("The task response time mean and variance must be positive.");
        }
        final double scv = VT / (ET * ET);
        double alpha;
        if (Math.abs(scv - 1.0) < FINE_TOL) {
            alpha = 1.0;
        } else {
            UnivariateFunction residual = new UnivariateFunction() {
                @Override
                public double value(double la) {
                    double a = Math.exp(la);
                    double m = Gamma.digamma(a + 1.0) - Gamma.digamma(1.0);
                    return (Gamma.trigamma(1.0) - Gamma.trigamma(a + 1.0)) / (m * m) - scv;
                }
            };
            double lo = -30.0;
            double hi = 30.0;
            while (residual.value(lo) < 0.0 && lo > -700.0) {
                lo -= 30.0;   // smaller shape -> larger SCV
            }
            while (residual.value(hi) > 0.0 && hi < 700.0) {
                hi += 30.0;   // larger shape -> smaller SCV
            }
            alpha = Math.exp(new BrentSolver(1e-12, 1e-12).solve(200, residual, lo, hi));
        }
        double beta = ET / (Gamma.digamma(alpha + 1.0) - Gamma.digamma(1.0));
        return new GEFit(alpha, beta);
    }

    private static double normalizePercentile(double p) {
        double pp = p > 1.0 ? p / 100.0 : p;
        if (pp <= 0.0 || pp >= 1.0) {
            throw new IllegalArgumentException(
                "The percentile must lie strictly between 0 and 1 (or 0 and 100). Got " + p + ".");
        }
        return pp;
    }
}
