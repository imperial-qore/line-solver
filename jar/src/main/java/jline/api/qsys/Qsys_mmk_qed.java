/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.qsys;

import java.util.HashMap;
import java.util.Map;
import java.util.function.DoubleUnaryOperator;

import org.apache.commons.math3.special.Erf;

/**
 * Halfin-Whitt QED approximation for the M/M/s queue, and the square-root
 * staffing rule that inverts it.
 *
 * <p>THE REGIME. Let s grow with the offered load a = lambda/mu so that the
 * SERVER SLACK stays of order sqrt(s), i.e. beta = (1-rho)sqrt(s) = (s-a)/sqrt(s)
 * is held fixed. The delay probability then has the non-degenerate limit
 *
 * <pre>  alpha(beta) = [ 1 + beta Phi(beta)/phi(beta) ]^(-1)</pre>
 *
 * with phi and Phi the standard normal density and cdf. That is the point of the
 * regime: servers are busy a fraction 1 - beta/sqrt(s) of the time, so efficiency
 * tends to 1, and yet the delay probability tends to a constant strictly between
 * 0 and 1, so quality does not collapse.
 *
 * <p>Useful even though M/M/s is exactly solvable, because Erlang C needs a sum
 * of s terms a^j/j! that overflows in double precision well before the thousands
 * of servers a large contact centre or thread pool has.
 *
 * <p>Port of MATLAB qsys_mmk_qed.m, qsys_mmk_qed_alpha.m and
 * qsys_mmk_qed_staffing.m.
 *
 * <p>Reference: S. Halfin, W. Whitt (1981). Heavy-traffic limits for queues with
 * many exponential servers. Operations Research 29(3), 567-588.
 *
 * @since LINE 3.1.0
 */
public final class Qsys_mmk_qed {

    private Qsys_mmk_qed() {
    }

    /**
     * The Halfin-Whitt delay-probability function alpha(beta).
     *
     * <p>Evaluated as phi/(phi + beta*Phi) rather than as the reciprocal of
     * 1 + beta*Phi/phi: the two are the same function, but the quotient Phi/phi
     * overflows once phi underflows (beta beyond about 38), whereas this form
     * degrades to 0/(0+beta) = 0, the correct limit.
     *
     * @param beta the QED server-slack parameter
     * @return alpha(beta), and 1 for a non-positive beta
     */
    public static double qsys_mmk_qed_alpha(double beta) {
        if (beta <= 0) {
            return 1.0;
        }
        double phi = Math.exp(-beta * beta / 2.0) / Math.sqrt(2.0 * Math.PI);
        double Phi = Erf.erfc(-beta / Math.sqrt(2.0)) / 2.0;
        return phi / (phi + beta * Phi);
    }

    /**
     * QED approximation for the M/M/s queue.
     *
     * @param lambda arrival rate
     * @param mu     service rate of one server
     * @param s      number of servers
     * @return map with offeredLoad, trafficIntensity, beta, probDelay,
     *         meanWaitDelayed, meanWait, meanQueueLength, meanNumber and
     *         utilization; an overloaded model has probDelay 1 and infinite waits
     */
    public static Map<String, Double> qsys_mmk_qed(double lambda, double mu, int s) {
        if (lambda <= 0) {
            throw new RuntimeException("qsys_mmk_qed: the arrival rate lambda must be positive");
        }
        if (mu <= 0) {
            throw new RuntimeException("qsys_mmk_qed: the service rate mu must be positive");
        }
        if (s < 1) {
            throw new RuntimeException("qsys_mmk_qed: the number of servers s must be at least 1");
        }
        double a = lambda / mu;
        double rho = a / s;
        double beta = (s - a) / Math.sqrt(s);
        Map<String, Double> res = new HashMap<String, Double>();
        res.put("offeredLoad", a);
        res.put("trafficIntensity", rho);
        res.put("beta", beta);
        res.put("utilization", rho);
        if (beta <= 0) {
            res.put("probDelay", 1.0);
            res.put("meanWaitDelayed", Double.POSITIVE_INFINITY);
            res.put("meanWait", Double.POSITIVE_INFINITY);
            res.put("meanQueueLength", Double.POSITIVE_INFINITY);
            res.put("meanNumber", Double.POSITIVE_INFINITY);
            return res;
        }
        double alpha = qsys_mmk_qed_alpha(beta);
        double wDelayed = 1.0 / (s * mu - lambda);
        res.put("probDelay", alpha);
        res.put("meanWaitDelayed", wDelayed);
        res.put("meanWait", alpha * wDelayed);
        res.put("meanQueueLength", lambda * alpha * wDelayed);
        res.put("meanNumber", a + lambda * alpha * wDelayed);
        return res;
    }

    /**
     * Square-root staffing for a delay-probability target.
     *
     * @param lambda arrival rate
     * @param mu     service rate of one server
     * @param target the largest acceptable P(W &gt; 0), in (0,1)
     * @return map with numServers, beta, betaTarget, offeredLoad, probDelay and
     *         meanWait
     */
    public static Map<String, Double> qsys_mmk_qed_staffing(double lambda, double mu,
                                                            double target) {
        return qsys_mmk_qed_staffing(lambda, mu, target, "delay", Double.NaN, Double.NaN, false);
    }

    /**
     * Square-root staffing under one of the three criteria.
     *
     * <p>Invert alpha(beta) = target for the server slack and staff
     * s = ceil(a + beta sqrt(a)) with a = lambda/mu: the base a erlangs of work
     * plus a cushion that grows only as the square root of the load.
     *
     * @param lambda    arrival rate
     * @param mu        service rate of one server
     * @param target    the largest acceptable P(W &gt; 0) for "delay", the largest
     *                  acceptable E[W] for "meanwait", ignored for "servicelevel"
     * @param criterion "delay", "meanwait" or "servicelevel"
     * @param deadline  the deadline of the service-level criterion
     * @param level     the probability the deadline must be met with
     * @param exact     walk s until the EXACT Erlang C measure meets the target
     * @return map with numServers, beta, betaTarget, offeredLoad, probDelay,
     *         meanWait, exactUsed and, for "servicelevel", serviceLevel
     */
    public static Map<String, Double> qsys_mmk_qed_staffing(double lambda, double mu, double target,
                                                            String criterion, double deadline,
                                                            double level, boolean exact) {
        if (lambda <= 0) {
            throw new RuntimeException("qsys_mmk_qed_staffing: the arrival rate lambda must be positive");
        }
        if (mu <= 0) {
            throw new RuntimeException("qsys_mmk_qed_staffing: the service rate mu must be positive");
        }
        final double a = lambda / mu;
        final String crit = criterion == null ? "delay" : criterion.toLowerCase();
        final double fmu = mu;
        final double ftarget = target;
        final double fdeadline = deadline;
        final double flevel = level;
        double betaTarget;
        if ("delay".equals(crit)) {
            if (!(target > 0 && target < 1)) {
                throw new RuntimeException(
                        "qsys_mmk_qed_staffing: for the delay criterion the target must be in (0,1)");
            }
            betaTarget = solve(new DoubleUnaryOperator() {
                @Override
                public double applyAsDouble(double b) {
                    return ftarget - qsys_mmk_qed_alpha(b);
                }
            });
        } else if ("meanwait".equals(crit)) {
            if (!(target > 0)) {
                throw new RuntimeException(
                        "qsys_mmk_qed_staffing: for the meanwait criterion the target must be positive");
            }
            // The residual is written target - E[W] so that it increases in beta.
            betaTarget = solve(new DoubleUnaryOperator() {
                @Override
                public double applyAsDouble(double b) {
                    return ftarget - qsys_mmk_qed_alpha(b) / (fmu * b * Math.sqrt(a));
                }
            });
        } else if ("servicelevel".equals(crit)) {
            if (!(level > 0 && level < 1) || !(deadline > 0)) {
                throw new RuntimeException("qsys_mmk_qed_staffing: the service level must be in "
                        + "(0,1) and the deadline positive");
            }
            betaTarget = solve(new DoubleUnaryOperator() {
                @Override
                public double applyAsDouble(double b) {
                    double sApprox = a + b * Math.sqrt(a);
                    return (1.0 - qsys_mmk_qed_alpha(b)
                            * Math.exp(-fmu * b * Math.sqrt(sApprox) * fdeadline)) - flevel;
                }
            });
        } else {
            throw new RuntimeException("qsys_mmk_qed_staffing: unknown criterion " + criterion);
        }

        int s = (int) Math.max(1, Math.ceil(a + betaTarget * Math.sqrt(a)));
        if (s * mu <= lambda) {
            s = (int) Math.floor(a) + 1;
        }
        if (exact) {
            while (!meets(lambda, mu, s, target, crit, deadline, level)) {
                s++;
                if (s > 10000000) {
                    throw new RuntimeException("qsys_mmk_qed_staffing: the exact refinement passed "
                            + "10^7 servers without meeting the target");
                }
            }
            while (s > 1 && meets(lambda, mu, s - 1, target, crit, deadline, level)) {
                s--;
            }
        }

        Map<String, Double> qed = qsys_mmk_qed(lambda, mu, s);
        Map<String, Double> res = new HashMap<String, Double>();
        res.put("numServers", (double) s);
        res.put("beta", qed.get("beta"));
        res.put("betaTarget", betaTarget);
        res.put("offeredLoad", a);
        res.put("probDelay", qed.get("probDelay"));
        res.put("meanWait", qed.get("meanWait"));
        res.put("exactUsed", exact ? 1.0 : 0.0);
        if ("servicelevel".equals(crit)) {
            res.put("serviceLevel",
                    1.0 - qed.get("probDelay") * Math.exp(-(s * mu - lambda) * deadline));
        }
        return res;
    }

    /**
     * Erlang C by the recursion B_j = a B_(j-1)/(j + a B_(j-1)) on the Erlang B
     * blocking probability, which never forms a^j/j! and so never overflows.
     * That matters here: this is called at the s the staffing rule proposes,
     * routinely in the thousands, where the factorial form is already infinite.
     *
     * @param s      number of servers
     * @param lambda arrival rate
     * @param mu     service rate of one server
     * @return the probability that an arrival is delayed
     */
    public static double qsys_mmk_qed_erlangc(int s, double lambda, double mu) {
        double a = lambda / mu;
        double b = 1.0;
        for (int j = 1; j <= s; j++) {
            b = a * b / (j + a * b);
        }
        double rho = a / s;
        return rho >= 1.0 ? 1.0 : b / (1.0 - rho * (1.0 - b));
    }

    /** Bisection for a root of an increasing f on (0, hi]; the bracket grows. */
    private static double solve(DoubleUnaryOperator f) {
        double lo = 1e-9;
        double hi = 1.0;
        if (f.applyAsDouble(lo) > 0) {
            return lo;
        }
        while (f.applyAsDouble(hi) < 0) {
            hi *= 2.0;
            if (hi > 1e6) {
                throw new RuntimeException(
                        "qsys_mmk_qed_staffing: no server slack meets the target");
            }
        }
        for (int i = 0; i < 200; i++) {
            double mid = 0.5 * (lo + hi);
            if (f.applyAsDouble(mid) < 0) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        return 0.5 * (lo + hi);
    }

    /** The exact M/M/s measure against the target. */
    private static boolean meets(double lambda, double mu, int s, double target, String criterion,
                                 double deadline, double level) {
        if (s * mu <= lambda) {
            return false;
        }
        double c = qsys_mmk_qed_erlangc(s, lambda, mu);
        double wq = c / (s * mu - lambda);
        if ("delay".equals(criterion)) {
            return c <= target;
        }
        if ("meanwait".equals(criterion)) {
            return wq <= target;
        }
        if ("servicelevel".equals(criterion)) {
            return (1.0 - c * Math.exp(-(s * mu - lambda) * deadline)) >= level;
        }
        return false;
    }
}
