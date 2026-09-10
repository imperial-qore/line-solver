/**
 * @file M/G/1 queueing system analysis with SRPT scheduling
 *
 * @since LINE 3.0
 */
package jline.api.qsys;

import jline.io.Ret;
import jline.util.matrix.Matrix;

/**
 * Analyzes an M/G/1 queue with SRPT (Shortest Remaining Processing Time)
 * scheduling using the Schrage-Miller class-conditional response-time formula.
 *
 * <p>SRPT is a size-based policy that always serves the job with the smallest
 * remaining processing time, preempting whenever a shorter job arrives. For a
 * job of size {@code x} the mean response time is (Bansal-Harchol-Balter,
 * SIGMETRICS 2001, Sec. 4, Eqs (1)-(3), after Schrage-Miller 1966):</p>
 *
 * <pre>
 *   E[T(x)] = E[W(x)] + E[R(x)]
 *   E[W(x)] = lambda*(m2(x) + x^2*(1-F(x))) / (2*(1-rho(x))^2)   (waiting)
 *   E[R(x)] = integral_0^x dt/(1-rho(t))                          (residence)
 * </pre>
 *
 * where {@code f} is the overall (mixture) job-size density, {@code F} its CDF,
 * {@code rho(x)=lambda*int_0^x t f(t) dt} the load from jobs of size &le; x, and
 * {@code m2(x)=int_0^x t^2 f(t) dt}. The per-class mean is
 * {@code E[T_r]=int_0^inf E[T(x)] f_r(x) dx}; because {@code E[T(x)]} depends only
 * on the job size (SRPT is size-based, not class-based) this integral is exact.
 * Integrals are evaluated by cumulative trapezoidal quadrature on a common grid.
 * Each class is matched to its (mean=1/mu, scv=cs^2): exponential for cs=1, a
 * two-phase balanced hyperexponential for cs&gt;1, and a Tijms Erlang-(k-1)/Erlang-k
 * mixture for cs&lt;1. For the fully exponential case this reproduces the exact
 * M/M/1/SRPT hyperexponential-mixture result.
 */
public final class Qsys_mg1_srpt {
    private Qsys_mg1_srpt() {}

    /** Per-class job-size representation matched to (mean, scv). */
    private static final class SizeFit {
        int type;          // 0 = exp, 1 = h2 hyperexp, 2 = Erlang mixture
        double rate;       // exp rate, or common rate of the Erlang mixture
        double p;          // hyperexp / Erlang-mixture branch probability
        double r1, r2;     // hyperexp phase rates
        int k;             // Erlang mixture upper shape
        double rateMin, rateMax;
    }

    public static Ret.qsys_prio qsys_mg1_srpt(Matrix lambda, Matrix mu, Matrix cs) {
        double[] lambdaArr = lambda.toArray1D();
        double[] muArr = mu.toArray1D();
        double[] csArr = cs.toArray1D();

        if (!(lambdaArr.length == muArr.length && lambdaArr.length == csArr.length)) {
            throw new IllegalArgumentException("lambda, mu, and cs must have the same length");
        }

        int K = lambdaArr.length;
        for (int i = 0; i < K; i++) {
            if (!(lambdaArr[i] > 0.0)) {
                throw new IllegalArgumentException("lambda[" + i + "] must be positive");
            }
            if (!(muArr[i] > 0.0)) {
                throw new IllegalArgumentException("mu[" + i + "] must be positive");
            }
            if (!(csArr[i] >= 0.0)) {
                throw new IllegalArgumentException("cs[" + i + "] must be non-negative");
            }
        }

        double lambdaTotal = 0.0;
        double rhoUtil = 0.0;
        for (int i = 0; i < K; i++) {
            lambdaTotal += lambdaArr[i];
            rhoUtil += lambdaArr[i] / muArr[i];
        }
        if (!(rhoUtil < 1.0)) {
            throw new IllegalStateException("System is unstable: utilization rho = " + rhoUtil + " >= 1");
        }

        double[] pmix = new double[K];
        SizeFit[] fits = new SizeFit[K];
        double rateMin = Double.POSITIVE_INFINITY;
        double rateMax = 0.0;
        for (int r = 0; r < K; r++) {
            pmix[r] = lambdaArr[r] / lambdaTotal;
            fits[r] = buildFit(muArr[r], csArr[r]);
            if (fits[r].rateMin < rateMin) rateMin = fits[r].rateMin;
            if (fits[r].rateMax > rateMax) rateMax = fits[r].rateMax;
        }

        // Integration grid: 40 e-foldings of the slowest phase, at least 200
        // points per rate ratio to resolve the fastest phase.
        double xmax = 40.0 / rateMin;
        int N = (int) Math.min(2000000L,
                Math.max(20000L, (long) Math.ceil(200.0 * rateMax / rateMin)));
        double dx = xmax / N;

        double[] x = new double[N + 1];
        double[] ET = new double[N + 1];

        // First pass: mixture density/tail, truncated moments rho(x) and m2(x),
        // residence integral, and E[T(x)] on the grid (cumulative trapezoid).
        double cumTf = 0.0;    // int_0^x t f(t) dt
        double cumT2f = 0.0;   // int_0^x t^2 f(t) dt
        double cumRes = 0.0;   // int_0^x dt/(1-rho(t))
        double prevTf = 0.0, prevT2f = 0.0, prevInv = 0.0;
        for (int n = 0; n <= N; n++) {
            double xn = n * dx;
            x[n] = xn;
            double fmix = 0.0;
            double fbar = 0.0;
            for (int r = 0; r < K; r++) {
                fmix += pmix[r] * pdf(fits[r], xn);
                fbar += pmix[r] * tail(fits[r], xn);
            }
            double tf = xn * fmix;
            double t2f = xn * xn * fmix;
            if (n > 0) {
                cumTf += 0.5 * (tf + prevTf) * dx;
                cumT2f += 0.5 * (t2f + prevT2f) * dx;
            }
            prevTf = tf;
            prevT2f = t2f;

            double rhoX = lambdaTotal * cumTf;
            double denom = 1.0 - rhoX;
            if (denom < 1e-12) denom = 1e-12;

            double inv = 1.0 / denom;
            if (n > 0) {
                cumRes += 0.5 * (inv + prevInv) * dx;
            }
            prevInv = inv;

            double wait = lambdaTotal * (cumT2f + xn * xn * fbar) / (2.0 * denom * denom);
            ET[n] = wait + cumRes;
        }

        // Second pass: per-class mean E[T_r] = int E[T(x)] f_r(x) dx.
        double[] WArr = new double[K];
        for (int r = 0; r < K; r++) {
            double acc = 0.0;
            double prev = ET[0] * pdf(fits[r], x[0]);
            for (int n = 1; n <= N; n++) {
                double cur = ET[n] * pdf(fits[r], x[n]);
                acc += 0.5 * (cur + prev) * dx;
                prev = cur;
            }
            WArr[r] = acc;
        }

        // Load measure returned as rhohat = Q/(1+Q) (qsys convention).
        double Q = 0.0;
        for (int i = 0; i < K; i++) {
            Q += lambdaArr[i] * WArr[i];
        }
        double rhohat = Q / (1.0 + Q);

        return new Ret.qsys_prio(new Matrix(WArr), rhohat);
    }

    /** Match a job-size distribution to mean 1/mu and scv cs^2. */
    private static SizeFit buildFit(double mu, double cs) {
        SizeFit f = new SizeFit();
        double meanX = 1.0 / mu;
        double c2 = cs * cs;
        if (Math.abs(c2 - 1.0) < 1e-9) {
            f.type = 0;
            f.rate = mu;
            f.rateMin = mu;
            f.rateMax = mu;
        } else if (c2 > 1.0) {
            double pr = 0.5 * (1.0 + Math.sqrt((c2 - 1.0) / (c2 + 1.0)));
            double r1 = 2.0 * pr * mu;
            double r2 = 2.0 * (1.0 - pr) * mu;
            f.type = 1;
            f.p = pr;
            f.r1 = r1;
            f.r2 = r2;
            f.rateMin = Math.min(r1, r2);
            f.rateMax = Math.max(r1, r2);
        } else {
            int k = (int) Math.ceil(1.0 / c2);
            double pr = (1.0 / (1.0 + c2)) * (k * c2 - Math.sqrt(k * (1.0 + c2) - (double) k * k * c2));
            double rate = (k - pr) / meanX;
            f.type = 2;
            f.k = k;
            f.p = pr;
            f.rate = rate;
            f.rateMin = rate;
            f.rateMax = rate;
        }
        return f;
    }

    /** Job-size probability density at x. */
    private static double pdf(SizeFit f, double x) {
        switch (f.type) {
            case 0:
                return f.rate * Math.exp(-f.rate * x);
            case 1:
                return f.p * f.r1 * Math.exp(-f.r1 * x)
                        + (1.0 - f.p) * f.r2 * Math.exp(-f.r2 * x);
            default:
                return f.p * erlangPdf(f.k - 1, f.rate, x)
                        + (1.0 - f.p) * erlangPdf(f.k, f.rate, x);
        }
    }

    /** Complementary CDF P(X > x). */
    private static double tail(SizeFit f, double x) {
        switch (f.type) {
            case 0:
                return Math.exp(-f.rate * x);
            case 1:
                return f.p * Math.exp(-f.r1 * x) + (1.0 - f.p) * Math.exp(-f.r2 * x);
            default:
                return f.p * erlangTail(f.k - 1, f.rate, x)
                        + (1.0 - f.p) * erlangTail(f.k, f.rate, x);
        }
    }

    /**
     * Erlang-n density in log space: f(x) = rate * Poisson(n-1; rate*x).
     * n = 0 is a point mass at 0 (density 0 for x &gt; 0).
     */
    private static double erlangPdf(int n, double rate, double x) {
        if (n <= 0) {
            return 0.0;
        }
        double t = rate * x;
        int m = n - 1;
        if (t <= 0.0) {
            return (m == 0) ? rate : 0.0;
        }
        double logp = m * Math.log(t) - t - logGamma(m + 1.0);
        return rate * Math.exp(logp);
    }

    /**
     * Erlang-n complementary CDF P(X &gt; x) = sum_{j=0}^{n-1} Poisson(j; rate*x),
     * computed in log space. n = 0 tail is 0 for x &gt; 0.
     */
    private static double erlangTail(int n, double rate, double x) {
        if (n <= 0) {
            return 0.0;
        }
        double t = rate * x;
        if (t <= 0.0) {
            return 1.0; // all mass above 0
        }
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            double logp = j * Math.log(t) - t - logGamma(j + 1.0);
            sum += Math.exp(logp);
        }
        return sum;
    }

    /** Lanczos approximation of log Gamma(x) for x &gt; 0. */
    private static double logGamma(double x) {
        double[] c = {
                676.5203681218851, -1259.1392167224028, 771.32342877765313,
                -176.61502916214059, 12.507343278686905, -0.13857109526572012,
                9.9843695780195716e-6, 1.5056327351493116e-7
        };
        if (x < 0.5) {
            // Reflection formula
            return Math.log(Math.PI / Math.sin(Math.PI * x)) - logGamma(1.0 - x);
        }
        double xx = x - 1.0;
        double a = 0.99999999999980993;
        double tt = xx + 7.5;
        for (int i = 0; i < c.length; i++) {
            a += c[i] / (xx + i + 1.0);
        }
        return 0.5 * Math.log(2.0 * Math.PI) + (xx + 0.5) * Math.log(tt) - tt + Math.log(a);
    }
}
