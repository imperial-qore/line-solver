/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.qsys;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.special.Gamma;

/**
 * Extremal two-moment bounds for the GI/GI/1 queue.
 *
 * <p>WHAT THE INTERVAL MEANS. Two moments do not determine E[W]; they determine
 * a SET of possible values, and the width of that set is the honest uncertainty
 * in any two-moment approximation. The extremal distributions attain its ends:
 *
 * <ul>
 *   <li>the lower end with deterministic interarrival times and a three-point
 *       service law on multiples of that interval, closed form
 *       rho((1+cs^2)rho-1)^+/(2(1-rho)) (eq. 2.12);</li>
 *   <li>the upper end asymptotically with TWO-POINT laws: an interarrival law
 *       with an atom at 0, and a service law whose upper atom runs off to
 *       infinity while its probability vanishes.</li>
 * </ul>
 *
 * <p>Making an interarrival time larger only empties the queue once, but making
 * a service time larger delays every customer behind it, which is why the two
 * ends look so different.
 *
 * <p>HOW THE UPPER END IS COMPUTED. Chen and Whitt reduce that limit to a
 * D(1/p)/RS(D(rho),p)/1 model with p = 1/(1+ca^2) and RS a geometric random sum,
 * then evaluate its mean waiting time by Spitzer's identity with the negative
 * binomial pmf (their Algorithm 1). The sum is truncated in both indices, so
 * this bound is a numerical limit, not a formula; the closed form of eq. (3.4)
 * is within about 1% of it.
 *
 * <p>Port of MATLAB qsys_gig1_bnds_extremal.m. All four codebases form the pmf
 * in LOG space rather than by the ratio recursion Algorithm 1 prints: the two
 * are equivalent, but the recursion accumulates rounding over thousands of
 * multiplications and would leave the ports agreeing only to about 1e-4.
 *
 * <p>Reference: Y. Chen, W. Whitt (2020). Algorithms for the upper bound mean
 * waiting time in the GI/GI/1 queue. Queueing Systems 94, 327-356.
 *
 * @since LINE 3.1.0
 */
public final class Qsys_gig1_bnds_extremal {

    private Qsys_gig1_bnds_extremal() {
    }

    /** Truncation of the negative binomial value. */
    public static final int DEFAULT_K = 4000;
    /** Truncation of the random-walk length. */
    public static final int DEFAULT_N = 2000;

    /**
     * Every bound at the default truncations.
     *
     * @param lambda arrival rate
     * @param mu     service rate
     * @param ca     coefficient of variation of the interarrival time
     * @param cs     coefficient of variation of the service time
     * @return map of times in queue, keyed as in the MATLAB struct
     */
    public static Map<String, Double> qsys_gig1_bnds_extremal(double lambda, double mu, double ca,
                                                               double cs) {
        return qsys_gig1_bnds_extremal(lambda, mu, ca, cs, DEFAULT_K, DEFAULT_N, false);
    }

    /**
     * @param lambda    arrival rate
     * @param mu        service rate
     * @param ca        coefficient of variation of the interarrival time
     * @param cs        coefficient of variation of the service time
     * @param K         truncation of the negative binomial value
     * @param N         truncation of the random-walk length
     * @param skipTight skip the O(K*N) tight bound and return the closed forms only
     * @return map with trafficIntensity, lowerBound, upperBound, upperBoundClosed,
     *         upperBoundDaley, upperBoundKingman, heavyTraffic, delta,
     *         relativeWidth and tightComputed; every waiting time is a TIME IN
     *         QUEUE, so add 1/mu for a response time
     */
    public static Map<String, Double> qsys_gig1_bnds_extremal(double lambda, double mu, double ca,
                                                               double cs, int K, int N,
                                                               boolean skipTight) {
        if (lambda <= 0 || mu <= 0) {
            throw new RuntimeException(
                    "qsys_gig1_bnds_extremal: the arrival and service rates must be positive");
        }
        double rho = lambda / mu;
        if (rho >= 1) {
            throw new RuntimeException(
                    "qsys_gig1_bnds_extremal: the bounds require a stable queue, rho < 1");
        }
        double ca2 = ca * ca;
        double cs2 = cs * cs;
        // The reference sets E[U] = 1, so every waiting time carries the factor
        // 1/lambda, that time unit expressed in the caller's units.
        double scale = 1.0 / lambda;

        double delta = delta(rho);                                                       // (3.5)
        Map<String, Double> res = new HashMap<String, Double>();
        res.put("trafficIntensity", rho);
        res.put("lowerBound", scale * rho * Math.max((1 + cs2) * rho - 1, 0.0) / (2 * (1 - rho)));
        res.put("upperBoundKingman", scale * rho * rho * (ca2 / (rho * rho) + cs2) / (2 * (1 - rho)));
        res.put("upperBoundDaley",
                scale * rho * rho * ((2 - rho) * ca2 / rho + cs2) / (2 * (1 - rho)));
        res.put("heavyTraffic", scale * rho * rho * (ca2 + cs2) / (2 * (1 - rho)));
        res.put("delta", delta);
        res.put("upperBoundClosed",
                scale * (2 * (1 - rho) * rho / (1 - delta) * ca2 + rho * rho * cs2) / (2 * (1 - rho)));
        if (skipTight) {
            res.put("upperBound", res.get("upperBoundClosed"));
            res.put("tightComputed", 0.0);
        } else {
            res.put("upperBound", scale * tight(rho, ca2, cs2, K, N));
            res.put("tightComputed", 1.0);
        }
        double ub = res.get("upperBound");
        res.put("relativeWidth", ub > 0 ? (ub - res.get("lowerBound")) / ub : 0.0);
        return res;
    }

    /**
     * The D/M/1 root of eq. (3.5), delta = exp(-(1-delta)/rho), in (0,1).
     *
     * <p>g(delta) = delta - exp(-(1-delta)/rho) is negative at 0 and positive
     * just below 1, where the second root delta = 1 sits, so bisection on [0,1)
     * finds the wanted root without landing on the trivial one.
     *
     * @param rho traffic intensity
     * @return the root in (0,1)
     */
    public static double delta(double rho) {
        double lo = 0.0;
        double hi = 1.0 - 1e-15;
        for (int i = 0; i < 200; i++) {
            double mid = 0.5 * (lo + hi);
            if (mid - Math.exp(-(1 - mid) / rho) < 0) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        return 0.5 * (lo + hi);
    }

    /**
     * Algorithm 1 of the reference: the mean waiting time of the extremal model,
     * rho*ca2 + rho^2 cs2/(2(1-rho)) + E[W(D(1/p),RS(D(rho),p))], the last term
     * by Spitzer's identity sum_n E[Sn^+]/n with Sn = rho(NB(n,1-p)+n) - n/p.
     *
     * <p>The negative binomial pmf is formed in LOG space,
     * log P(NB(n,1-p)=k) = lgamma(n+k) - lgamma(k+1) - lgamma(n) + n log p + k log(1-p),
     * from one precomputed table of log-gammas. Algorithm 1 as printed steps the
     * pmf by the ratio P(n+1)/P(n) = ((n+k)/n)p, which is equivalent but
     * accumulates rounding over thousands of multiplications; the log form keeps
     * every codebase on the identical expression, which is what makes the four
     * ports comparable at 1e-12 rather than at the 1e-4 the recursion reaches.
     */
    private static double tight(double rho, double ca2, double cs2, int K, int N) {
        double p = 1.0 / (1.0 + ca2);
        double total = rho * ca2 + rho * rho * cs2 / (2.0 * (1.0 - rho));
        double logp = Math.log(p);
        double log1mp = Math.log1p(-p);
        double[] lg = new double[N + K + 2];
        for (int i = 1; i < lg.length; i++) {
            lg[i] = Gamma.logGamma(i);
        }
        for (int k = 1; k <= K; k++) {
            double s = 0.0;
            for (int n = 1; n <= N; n++) {
                double step = (n + k) * rho - n / p;
                if (step > 0) {
                    double lpmf = lg[n + k] - lg[k + 1] - lg[n] + n * logp + k * log1mp;
                    s += Math.exp(lpmf) * step / n;
                }
            }
            total += s;
        }
        return total;
    }
}
