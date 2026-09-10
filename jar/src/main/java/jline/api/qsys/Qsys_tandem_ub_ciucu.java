/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.qsys;

import java.util.function.DoubleUnaryOperator;

/**
 * Tail bounds for a GI/Hn/1 -&gt; ./Hn/1 tandem of two FCFS single servers.
 *
 * <p>Port of matlab/src/api/qsys/qsys_tandem_ub_ciucu.m. Both stations serve the
 * same hyperexponential law {@code Y, Z ~ sum_i p_i Exp(mu_i)}, a scalar
 * {@code p = 1} giving exponential service, and the arrivals are renewal with a
 * light-tailed interarrival time supplied through its Laplace-Stieltjes
 * transform {@code E[e^{-s X}]}.
 *
 * <p>With theta the positive root of {@code E[e^{theta (Y-X)}] = 1} and
 * {@code alpha = E[X e^{-theta X}]}, the test function
 *
 * <pre>
 *   gamma(u,v) = 1{0&lt;=u&lt;=v} [1 - A e^{-theta u} - (B + C u + D v) e^{-theta v}]
 * </pre>
 *
 * <p>satisfies the integral inequality of Theorem 1(b) of the reference once the
 * five sufficient conditions of its Lemma 4 fix A, B, C and D; Corollary 2 then
 * turns gamma into the tail bounds returned here. The two exponentials mix a
 * polynomial of degree one in x, which is what lets the bound follow the concave
 * bend of the tail on a linear-log scale where a purely exponential bound cannot.
 *
 * <p>In the M/M/1 -&gt; ./M/1 case the five inequalities hold as equalities, so
 * gamma is the exact joint distribution and both bounds are exact:
 * {@code P(S > x) = (1 + theta x) e^{-theta x}} and
 * {@code P(W > x) = (1 - 2 theta^2/(mu(mu+theta)) + x(mu-theta)theta/(mu+theta)) e^{-theta x}}.
 * Away from it the bound stays sharp: against an exact CTMC reference for the
 * Erlang(2)/M/1 -&gt; ./M/1 tandem it is within 2% at P(S&gt;x) = 1e-2 and within
 * 0.6% at 5e-10, with the correct asymptotic slope
 * {@code theta^2/(mu(1-alpha mu))}. Accuracy degrades with service variability,
 * to about a factor of two at CV(Y) = 2.
 *
 * <p>Reference: F. Ciucu, S. Mehri, "On the Distribution of Sojourn Times in
 * Tandem Queues", Proc. ACM Meas. Anal. Comput. Syst. 9(2), Article 27, 2025
 * (ACM SIGMETRICS 2025). Registered in .citations() as 'tandemub'.
 *
 * @since LINE 3.1.0
 */
public final class Qsys_tandem_ub_ciucu {
    private Qsys_tandem_ub_ciucu() {}

    /**
     * Bounds with the derivative of the transform obtained numerically.
     *
     * @param x   thresholds at which the tails are bounded, nonnegative
     * @param lst interarrival transform, {@code s -> E[e^{-s X}]} for s &gt;= 0
     * @param p   service phase probabilities, nonnegative and summing to one
     * @param mu  service phase rates, positive
     * @return the tail bounds and the coefficients behind them
     */
    public static QsysTandemUbResult qsys_tandem_ub_ciucu(double[] x, DoubleUnaryOperator lst,
                                                          double[] p, double[] mu) {
        return qsys_tandem_ub_ciucu(x, lst, p, mu, null);
    }

    /**
     * Bounds with an analytic derivative of the transform.
     *
     * @param x    thresholds at which the tails are bounded, nonnegative
     * @param lst  interarrival transform, {@code s -> E[e^{-s X}]} for s &gt;= 0
     * @param p    service phase probabilities, nonnegative and summing to one
     * @param mu   service phase rates, positive
     * @param dlst {@code s -> E[X e^{-s X}]}, minus the derivative of lst; null
     *             to obtain it by a Richardson-extrapolated central difference,
     *             which costs four extra transform evaluations and loses roughly
     *             four digits
     * @return the tail bounds and the coefficients behind them
     */
    public static QsysTandemUbResult qsys_tandem_ub_ciucu(double[] x, DoubleUnaryOperator lst,
                                                          double[] p, double[] mu,
                                                          DoubleUnaryOperator dlst) {
        if (x == null || x.length < 1) {
            throw new IllegalArgumentException("x must hold at least one threshold");
        }
        if (lst == null) {
            throw new IllegalArgumentException("lst must be supplied");
        }
        if (p == null || mu == null || p.length != mu.length || p.length < 1) {
            throw new IllegalArgumentException("p and mu must have the same number of phases");
        }
        double psum = 0.0;
        for (int i = 0; i < p.length; i++) {
            if (!(p[i] >= 0.0)) {
                throw new IllegalArgumentException("The phase probabilities p must be nonnegative");
            }
            if (!(mu[i] > 0.0)) {
                throw new IllegalArgumentException("The service rates mu must be positive");
            }
            psum += p[i];
        }
        if (Math.abs(psum - 1.0) > 1e-10) {
            throw new IllegalArgumentException("The phase probabilities p must sum to one");
        }
        for (int k = 0; k < x.length; k++) {
            if (!(x[k] >= 0.0)) {
                throw new IllegalArgumentException("The thresholds x must be nonnegative");
            }
        }

        double mu1 = Double.POSITIVE_INFINITY;
        for (int i = 0; i < mu.length; i++) {
            mu1 = Math.min(mu1, mu[i]);
        }
        // Stability: E[X] > E[Y] is what makes E[e^{t(Y-X)}] - 1 cross zero on (0,mu1).
        double hi = mu1 * (1.0 - 1e-12);
        if (residual(hi, lst, p, mu) <= 0.0) {
            throw new IllegalArgumentException("No positive root of E[e^{theta(Y-X)}]=1 below "
                    + "min(mu): the tandem is unstable or the service is not the lighter tail");
        }
        double lo = mu1 * 1e-12;
        while (residual(lo, lst, p, mu) >= 0.0 && lo > mu1 * 1e-16) {
            lo = lo / 10.0;                         // walk below the root at zero
        }
        if (residual(lo, lst, p, mu) >= 0.0) {
            // E[e^{t(Y-X)}]-1 is convex and vanishes at t=0, so it stays positive on
            // the whole of (0,mu1) exactly when its slope E[Y]-E[X] there is nonnegative.
            throw new IllegalArgumentException("The tandem is unstable, E[X] <= E[Y]: "
                    + "theta = 0 is the only root of E[e^{theta(Y-X)}]=1");
        }
        double theta = bisect(lo, hi, lst, p, mu);
        double alpha = dlst != null ? dlst.applyAsDouble(theta) : numericalDlst(lst, theta);

        double EexpZ = mgfY(theta, p, mu);          // E[e^{theta Z}]
        double EZexp = 0.0;                         // E[Z e^{theta Z}]
        double EY = 0.0;
        double sumPOverMuMinusTheta = 0.0;
        for (int i = 0; i < p.length; i++) {
            double dm = mu[i] - theta;
            EZexp += p[i] * mu[i] / (dm * dm);
            EY += p[i] / mu[i];
            sumPOverMuMinusTheta += p[i] / dm;
        }
        double A = 1.0;
        double C = theta * sumPOverMuMinusTheta / EZexp;
        double EUeV = EY - alpha * EexpZ;           // E[U e^{theta V}]
        double EVeV = EZexp / EexpZ - alpha * EexpZ; // E[V e^{theta V}] > 0
        double D = -C * EUeV / EVeV;
        if (!(D > 0.0)) {
            D = 0.0;
        }
        double B = D > 0.0 ? (C + D) / (mu1 - theta) - theta / mu1
                           : C * (1.0 / mu1 - alpha * EexpZ);

        double[] S = new double[x.length];
        for (int k = 0; k < x.length; k++) {
            double acc = 0.0;
            for (int i = 0; i < p.length; i++) {
                double m = mu[i];
                double dm = m - theta;
                double em = Math.exp(-m * x[k]);
                double et = Math.exp(-theta * x[k]);
                acc += p[i] * (em + m / dm * (A + B) * (et - em)
                        + m / (dm * dm) * (C + D) * ((dm * x[k] - 1.0) * et + em));
            }
            S[k] = Math.min(acc, 1.0);
        }

        double[] W = new double[x.length];
        if (p.length == 1) {
            double beta = lst.applyAsDouble(mu1);   // E[e^{-mu X}]
            for (int k = 0; k < x.length; k++) {
                double et = Math.exp(-theta * x[k]);
                double val;
                if (D == 0.0) {
                    val = (1.0 - 2.0 * theta * theta / (mu1 * (mu1 + theta))
                            + theta * (mu1 - theta) / (mu1 + theta) * x[k]) * et
                            + beta * (theta * mu1 * alpha / (2.0 * (mu1 - theta))
                            - theta / (2.0 * mu1)) * Math.exp(-mu1 * x[k]);
                } else {
                    val = (1.0 - 2.0 * theta / mu1
                            + 2.0 * theta * theta * (2.0 - alpha * mu1)
                              / ((mu1 + theta) * (mu1 + theta) * (1.0 - alpha * mu1))
                            + theta * theta * (mu1 - theta)
                              / (mu1 * (mu1 + theta) * (1.0 - alpha * mu1)) * x[k]) * et;
                }
                W[k] = Math.min(val, 1.0);
            }
        } else {
            for (int k = 0; k < x.length; k++) {
                W[k] = Double.NaN;                  // the W form is Exp-service only
            }
        }
        return new QsysTandemUbResult(S, W, theta, alpha, A, B, C, D);
    }

    /** E[e^{t Y}] of the hyperexponential service law, for t below every rate. */
    private static double mgfY(double t, double[] p, double[] mu) {
        double acc = 0.0;
        for (int i = 0; i < p.length; i++) {
            acc += p[i] * mu[i] / (mu[i] - t);
        }
        return acc;
    }

    /** The function whose positive root is theta. */
    private static double residual(double t, DoubleUnaryOperator lst, double[] p, double[] mu) {
        return mgfY(t, p, mu) * lst.applyAsDouble(t) - 1.0;
    }

    /** Bisection on a bracket with a sign change; deterministic and toolbox free. */
    private static double bisect(double lo, double hi, DoubleUnaryOperator lst,
                                 double[] p, double[] mu) {
        double a = lo;
        double b = hi;
        for (int it = 0; it < 200 && (b - a) > 1e-14 * Math.max(1.0, b); it++) {
            double m = 0.5 * (a + b);
            if (residual(m, lst, p, mu) > 0.0) {
                b = m;
            } else {
                a = m;
            }
        }
        return 0.5 * (a + b);
    }

    /** Richardson-extrapolated central difference of -lst at s, i.e. E[X e^{-s X}]. */
    private static double numericalDlst(DoubleUnaryOperator lst, double s) {
        double h = 1e-3 * (1.0 + s);
        if (h > s) {
            h = s / 2.0;
        }
        double d1 = (lst.applyAsDouble(s - h) - lst.applyAsDouble(s + h)) / (2.0 * h);
        double d2 = (lst.applyAsDouble(s - h / 2.0) - lst.applyAsDouble(s + h / 2.0)) / h;
        return (4.0 * d2 - d1) / 3.0;
    }
}
