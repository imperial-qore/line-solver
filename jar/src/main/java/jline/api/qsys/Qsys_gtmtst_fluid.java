/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.qsys;

import java.util.function.DoubleUnaryOperator;

/**
 * The Gt/Mt/st+GI many-server fluid queue.
 *
 * <p>Time-varying arrival rate lambda(t), time-varying staffing s(t),
 * exponential service at the time-varying rate mu(t), general patience with
 * complementary cdf F^c, and unlimited waiting room.
 *
 * <p>THE MODEL ALTERNATES BETWEEN TWO REGIMES, and the whole algorithm is the
 * bookkeeping of that alternation:
 *
 * <ul>
 *   <li>UNDERLOADED: the queue is empty and every arrival enters service at
 *       once, so the system is the infinite-server fluid model and B obeys
 *       B'(t) = lambda(t) - mu(t)B(t) (eq. 18 in its Mt form). It ends when B
 *       reaches s while lambda exceeds the rate Gamma(t) = s'(t) + s(t)mu(t) at
 *       which capacity frees up (eq. 15).</li>
 *   <li>OVERLOADED: every server is busy, B(t) = s(t), fluid enters service at
 *       exactly Gamma(t), and the queue is described by its BOUNDARY WAITING
 *       TIME w(t), the age of the oldest fluid still waiting. Content of age x
 *       is what arrived x ago and has not abandoned, q(t,x) = lambda(t-x)F^c(x),
 *       and the boundary moves by the delay differential equation (eq. 21)
 *       w'(t) = 1 - Gamma(t)/[lambda(t-w(t)) F^c(w(t))]. It ends when w returns
 *       to 0 with lambda no longer above Gamma (eq. 14).</li>
 * </ul>
 *
 * <p>WHY w AND NOT Q. The queue content is a functional of w, but not the other
 * way round: two systems with the same Q and different age profiles abandon at
 * different rates. Tracking the boundary keeps the age profile exact, which is
 * what makes a general patience law admissible at all.
 *
 * <p>Port of MATLAB qsys_gtmtst_fluid.m.
 *
 * <p>Reference: Y. Liu, W. Whitt (2012). The Gt/GI/st+GI many-server fluid
 * queue. Queueing Systems 71, 405-444; Y. Liu, W. Whitt (2014). Algorithms for
 * time-varying networks of many-server fluid queues. INFORMS Journal on
 * Computing 26(1), 59-73.
 *
 * @since LINE 3.1.0
 */
public final class Qsys_gtmtst_fluid {

    private Qsys_gtmtst_fluid() {
    }

    /**
     * The fluid queue started empty, on a grid of T/2000.
     *
     * @param lambdaFun     arrival rate lambda(t)
     * @param sFun          staffing s(t)
     * @param muFun         service rate mu(t)
     * @param patienceCcdf  F^c(x) = P(patience &gt; x)
     * @param T             horizon
     * @return the trajectory
     */
    public static QsysTvFluidResult qsys_gtmtst_fluid(DoubleUnaryOperator lambdaFun,
                                                       DoubleUnaryOperator sFun,
                                                       DoubleUnaryOperator muFun,
                                                       DoubleUnaryOperator patienceCcdf, double T) {
        return qsys_gtmtst_fluid(lambdaFun, sFun, muFun, patienceCcdf, T, T / 2000.0, 0.0, 0.0,
                null, null, null);
    }

    /**
     * @param lambdaFun    arrival rate lambda(t)
     * @param sFun         staffing s(t)
     * @param muFun        service rate mu(t)
     * @param patienceCcdf F^c(x) = P(patience &gt; x)
     * @param T            horizon
     * @param dt           grid step
     * @param B0           fluid in service at time 0
     * @param w0           boundary waiting time at time 0
     * @param sPrimeFun    s'(t), or null to differentiate sFun numerically
     * @param patiencePdf  the patience density, or null to difference the ccdf
     * @param lambdaPast   the arrival rate before time 0, or null for lambdaFun
     * @return the trajectory
     */
    public static QsysTvFluidResult qsys_gtmtst_fluid(DoubleUnaryOperator lambdaFun,
                                                       DoubleUnaryOperator sFun,
                                                       DoubleUnaryOperator muFun,
                                                       DoubleUnaryOperator patienceCcdf, double T,
                                                       double dt, double B0, double w0,
                                                       DoubleUnaryOperator sPrimeFun,
                                                       DoubleUnaryOperator patiencePdf,
                                                       DoubleUnaryOperator lambdaPast) {
        if (T <= 0) {
            throw new RuntimeException("qsys_gtmtst_fluid: the horizon T must be positive");
        }
        if (!(dt > 0)) {
            dt = T / 2000.0;
        }
        int n = (int) Math.round(T / dt) + 1;
        double[] t = new double[n];
        for (int i = 0; i < n; i++) {
            t[i] = T * i / (n - 1.0);
        }
        double step = t[1] - t[0];

        double[] lam = new double[n];
        double[] s = new double[n];
        double[] mu = new double[n];
        for (int i = 0; i < n; i++) {
            lam[i] = lambdaFun.applyAsDouble(t[i]);
            s[i] = sFun.applyAsDouble(t[i]);
            mu[i] = muFun.applyAsDouble(t[i]);
            if (s[i] <= 0) {
                throw new RuntimeException("qsys_gtmtst_fluid: the staffing function must be positive");
            }
            if (mu[i] <= 0) {
                throw new RuntimeException("qsys_gtmtst_fluid: the service rate must be positive");
            }
        }
        double[] sp = new double[n];
        if (sPrimeFun != null) {
            for (int i = 0; i < n; i++) {
                sp[i] = sPrimeFun.applyAsDouble(t[i]);
            }
        } else {
            for (int i = 0; i < n; i++) {
                if (i == 0) {
                    sp[i] = (s[1] - s[0]) / step;
                } else if (i == n - 1) {
                    sp[i] = (s[n - 1] - s[n - 2]) / step;
                } else {
                    sp[i] = (s[i + 1] - s[i - 1]) / (2 * step);
                }
            }
        }

        final DoubleUnaryOperator past = lambdaPast != null ? lambdaPast : lambdaFun;
        final DoubleUnaryOperator lamFun = lambdaFun;
        DoubleUnaryOperator lamOf = new DoubleUnaryOperator() {
            @Override
            public double applyAsDouble(double u) {
                return u < 0 ? past.applyAsDouble(u) : lamFun.applyAsDouble(u);
            }
        };
        final DoubleUnaryOperator ccdf = patienceCcdf;
        DoubleUnaryOperator pdf = patiencePdf;
        if (pdf == null) {
            pdf = new DoubleUnaryOperator() {
                @Override
                public double applyAsDouble(double x) {
                    double h = 1e-6;
                    return Math.max(0.0,
                            (ccdf.applyAsDouble(Math.max(0.0, x - h)) - ccdf.applyAsDouble(x + h))
                                    / (2 * h));
                }
            };
        }

        double[] B = new double[n];
        double[] Q = new double[n];
        double[] w = new double[n];
        double[] alpha = new double[n];
        int[] regime = new int[n];
        double[] gamma = new double[n];
        for (int i = 0; i < n; i++) {
            gamma[i] = sp[i] + s[i] * mu[i];             // Gamma(t), eq. (13)
        }
        B[0] = B0;
        w[0] = w0;
        boolean over = w0 > 0 || (B0 >= s[0] - 1e-12 && lam[0] > gamma[0]);
        regime[0] = over ? 1 : 0;
        if (over) {
            B[0] = s[0];
        }
        Q[0] = integrate(lamOf, ccdf, t[0], w[0], step);
        alpha[0] = integrate(lamOf, pdf, t[0], w[0], step);

        for (int i = 0; i < n - 1; i++) {
            double bNext;
            double wNext;
            if (regime[i] == 0) {
                // Underloaded: B' = lambda - mu B, by RK4 on the grid step.
                double k1 = interp(t, lam, t[i]) - interp(t, mu, t[i]) * B[i];
                double k2 = interp(t, lam, t[i] + step / 2)
                        - interp(t, mu, t[i] + step / 2) * (B[i] + step * k1 / 2);
                double k3 = interp(t, lam, t[i] + step / 2)
                        - interp(t, mu, t[i] + step / 2) * (B[i] + step * k2 / 2);
                double k4 = interp(t, lam, t[i] + step)
                        - interp(t, mu, t[i] + step) * (B[i] + step * k3);
                bNext = B[i] + step * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
                wNext = 0.0;
                if (bNext >= s[i + 1] && lam[i + 1] > gamma[i + 1]) {
                    // The servers just filled and the input outruns the freed
                    // capacity: eq. (15), the underloaded interval ends here.
                    bNext = s[i + 1];
                    regime[i + 1] = 1;
                } else {
                    regime[i + 1] = 0;
                    bNext = Math.min(bNext, s[i + 1]);
                }
            } else {
                // Overloaded: B = s and the boundary moves by eq. (21).
                double k1 = wdot(t[i], w[i], lamOf, ccdf, t, gamma);
                double k2 = wdot(t[i] + step / 2, Math.max(0.0, w[i] + step * k1 / 2), lamOf, ccdf,
                        t, gamma);
                double k3 = wdot(t[i] + step / 2, Math.max(0.0, w[i] + step * k2 / 2), lamOf, ccdf,
                        t, gamma);
                double k4 = wdot(t[i] + step, Math.max(0.0, w[i] + step * k3), lamOf, ccdf, t,
                        gamma);
                wNext = w[i] + step * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
                bNext = s[i + 1];
                if (wNext <= 0 && lam[i + 1] <= gamma[i + 1]) {
                    // The queue has drained and the input no longer outruns the
                    // freed capacity: eq. (14), the overloaded interval ends here.
                    wNext = 0.0;
                    regime[i + 1] = 0;
                } else {
                    wNext = Math.max(wNext, 0.0);
                    regime[i + 1] = 1;
                }
            }
            B[i + 1] = bNext;
            w[i + 1] = wNext;
            if (regime[i + 1] == 1) {
                Q[i + 1] = integrate(lamOf, ccdf, t[i + 1], wNext, step);
                alpha[i + 1] = integrate(lamOf, pdf, t[i + 1], wNext, step);
            }
        }

        double[] sigma = new double[n];
        double[] util = new double[n];
        double[] x = new double[n];
        double[] entry = new double[n];
        for (int i = 0; i < n; i++) {
            sigma[i] = mu[i] * B[i];                    // service completion rate, eq. (3)
            util[i] = B[i] / s[i];
            x[i] = B[i] + Q[i];
            entry[i] = t[i] - w[i];
        }
        // The potential waiting time of an arrival at t is the u-t at which the
        // boundary reaches it, i.e. the solution of u - w(u) = t. That map is
        // non-decreasing, so one interpolation inverts it.
        double[] v = new double[n];
        for (int i = 0; i < n; i++) {
            v[i] = Math.max(0.0, interp(entry, t, t[i]) - t[i]);
        }

        return new QsysTvFluidResult(t, regime, B, Q, x, w, v, sigma, alpha, util, lam, s, gamma);
    }

    /**
     * int_0^w lambda(t-x) WEIGHT(x) dx by Simpson: the ccdf gives the queue
     * content that has not abandoned, the density gives the abandonment rate.
     */
    private static double integrate(DoubleUnaryOperator lamOf, DoubleUnaryOperator weight,
                                    double ti, double wi, double dt) {
        if (wi <= 0) {
            return 0.0;
        }
        int m = Math.max(8, (int) Math.ceil(wi / dt) + 1);
        if (m % 2 == 1) {
            m++;
        }
        double h = wi / m;
        double sum = lamOf.applyAsDouble(ti) * weight.applyAsDouble(0.0)
                + lamOf.applyAsDouble(ti - wi) * weight.applyAsDouble(wi);
        for (int j = 1; j < m; j++) {
            double xx = j * h;
            sum += (j % 2 == 1 ? 4.0 : 2.0) * lamOf.applyAsDouble(ti - xx)
                    * weight.applyAsDouble(xx);
        }
        return h / 3.0 * sum;
    }

    /** Eq. (21): w' = 1 - Gamma(t)/q~(t,w) with q~(t,w) = lambda(t-w)F^c(w). */
    private static double wdot(double tt, double ww, DoubleUnaryOperator lamOf,
                               DoubleUnaryOperator ccdf, double[] t, double[] gamma) {
        double den = lamOf.applyAsDouble(tt - ww) * ccdf.applyAsDouble(ww);
        if (den <= 0) {
            // No fluid of that age survives, so the boundary can only advance
            // with the clock.
            return 1.0;
        }
        return 1.0 - interp(t, gamma, tt) / den;
    }

    /**
     * Linear interpolation on an increasing grid, clamped at both ends. Public
     * because the network solver feeds each queue its interpolated arrival rate.
     *
     * @param xs the grid, increasing
     * @param ys the values on that grid
     * @param x  where to interpolate
     * @return the interpolated value
     */
    public static double interp(double[] xs, double[] ys, double x) {
        int n = xs.length;
        if (x <= xs[0]) {
            return ys[0];
        }
        if (x >= xs[n - 1]) {
            return ys[n - 1];
        }
        int lo = 0;
        int hi = n - 1;
        while (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            if (xs[mid] <= x) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        double f = (x - xs[lo]) / (xs[hi] - xs[lo]);
        return ys[lo] + f * (ys[hi] - ys[lo]);
    }
}
