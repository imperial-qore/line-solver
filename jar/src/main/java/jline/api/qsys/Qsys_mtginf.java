/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.qsys;

import java.util.function.DoubleUnaryOperator;

/**
 * Exact time-varying analysis of the Mt/G/infinity queue.
 *
 * <p>THE RESULT IS EXACT, not an approximation. With infinitely many servers
 * customers never interact, so the model is a Poisson random measure and the
 * number in system at time t is POISSON with mean
 *
 * <pre>  m(t) = E[ int_{t-S}^{t} lambda(u) du ] = ES E[lambda(t - Se)]
 *       = int_0^Inf lambda(t-x) P(S &gt; x) dx</pre>
 *
 * where Se is the STATIONARY-EXCESS (equilibrium) law of the service time, with
 * density P(S&gt;x)/ES. Because the law is Poisson the variance equals the mean.
 *
 * <p>THE PHYSICS. Reading m(t) as ES E[lambda(t-Se)] says the time-varying load
 * is the stationary load ES lambda(t) subjected to a TIME LAG and a SPACE SHIFT:
 * to first order m(t) ~ ES lambda(t - E[Se]) with E[Se] = E[S^2]/(2 ES), so peak
 * congestion LAGS peak arrival rate, and by more than the mean service time when
 * the service law is variable. The pointwise stationary approximation
 * ES lambda(t) is the zeroth-order term of the same expansion.
 *
 * <p>Port of MATLAB qsys_mtginf.m.
 *
 * <p>Reference: S. G. Eick, W. A. Massey, W. Whitt (1993). The physics of the
 * Mt/G/infinity queue. Operations Research 41(4), 731-742.
 *
 * @since LINE 3.1.0
 */
public final class Qsys_mtginf {

    private Qsys_mtginf() {
    }

    /** Service-tail cut for the age integral. */
    public static final double DEFAULT_TOL = 1e-12;
    /** Simpson panels for the age integral. */
    public static final int DEFAULT_PANELS = 4000;

    /**
     * The queue with an infinite past and the departure rate by flow balance.
     *
     * @param lambdaFun    the arrival rate; must accept arguments in the past
     * @param serviceCcdf  G^c(x) = P(S &gt; x)
     * @param ES           the mean service time
     * @param tvals        the times at which to evaluate
     * @return the time-varying measures
     */
    public static QsysMtginfResult qsys_mtginf(DoubleUnaryOperator lambdaFun,
                                                DoubleUnaryOperator serviceCcdf, double ES,
                                                double[] tvals) {
        return qsys_mtginf(lambdaFun, serviceCcdf, ES, tvals, Double.NEGATIVE_INFINITY, Double.NaN,
                null, DEFAULT_TOL, DEFAULT_PANELS, 1e12);
    }

    /**
     * @param lambdaFun   the arrival rate
     * @param serviceCcdf G^c(x) = P(S &gt; x)
     * @param ES          the mean service time
     * @param tvals       the times at which to evaluate
     * @param startTime   time the system started empty; -Inf assumes an infinite past
     * @param ES2         the second moment of the service time, NaN to skip the lag
     * @param servicePdf  the service density for the exact departure rate, or null
     * @param tol         service-tail cut for the age integral
     * @param panels      Simpson panels for that integral
     * @param maxAge      cap on the age integrated over
     * @return the time-varying measures
     */
    public static QsysMtginfResult qsys_mtginf(DoubleUnaryOperator lambdaFun,
                                                DoubleUnaryOperator serviceCcdf, double ES,
                                                double[] tvals, double startTime, double ES2,
                                                DoubleUnaryOperator servicePdf, double tol,
                                                int panels, double maxAge) {
        if (ES <= 0) {
            throw new RuntimeException("qsys_mtginf: the mean service time ES must be positive");
        }
        if (lambdaFun == null || serviceCcdf == null) {
            throw new RuntimeException("qsys_mtginf: the arrival rate and the service ccdf must not be null");
        }
        double cut = tailCut(serviceCcdf, tol, maxAge);
        boolean unbounded = Double.isInfinite(startTime);
        double[] t = tvals.clone();

        double[] mean = meanCurve(lambdaFun, serviceCcdf, t, startTime, cut, panels, unbounded);
        double[] arrival = new double[t.length];
        double[] psa = new double[t.length];
        for (int i = 0; i < t.length; i++) {
            arrival[i] = lambdaFun.applyAsDouble(t[i]);
            psa[i] = ES * arrival[i];
        }

        double[] departure = new double[t.length];
        if (servicePdf != null) {
            for (int i = 0; i < t.length; i++) {
                double hi = unbounded ? cut : Math.min(cut, Math.max(0.0, t[i] - startTime));
                double[][] xw = simpson(0.0, hi, panels);
                double acc = 0.0;
                for (int j = 0; j < xw[0].length; j++) {
                    acc += xw[1][j] * lambdaFun.applyAsDouble(t[i] - xw[0][j])
                            * servicePdf.applyAsDouble(xw[0][j]);
                }
                departure[i] = acc;
            }
        } else {
            // Flow balance m'(t) = lambda(t) - delta(t), differentiated centrally.
            double tmax = 1.0;
            for (int i = 0; i < t.length; i++) {
                tmax = Math.max(tmax, Math.abs(t[i]));
            }
            double h = 1e-5 * tmax;
            double[] tu = new double[t.length];
            double[] td = new double[t.length];
            for (int i = 0; i < t.length; i++) {
                tu[i] = t[i] + h;
                td[i] = t[i] - h;
            }
            double[] up = meanCurve(lambdaFun, serviceCcdf, tu, startTime, cut, panels, unbounded);
            double[] dn = meanCurve(lambdaFun, serviceCcdf, td, startTime, cut, panels, unbounded);
            for (int i = 0; i < t.length; i++) {
                departure[i] = arrival[i] - (up[i] - dn[i]) / (2.0 * h);
            }
        }

        double lag = Double.NaN;
        double[] lagApprox = null;
        if (!Double.isNaN(ES2)) {
            lag = ES2 / (2.0 * ES);                 // E[Se], the time lag
            lagApprox = new double[t.length];
            for (int i = 0; i < t.length; i++) {
                lagApprox[i] = ES * lambdaFun.applyAsDouble(t[i] - lag);
            }
        }

        return new QsysMtginfResult(t, mean, mean.clone(), arrival, departure, psa, lag, lagApprox);
    }

    /**
     * The Poisson mean m(t), shared by the public entry point and the
     * finite-difference departure rate so that neither re-derives the other.
     */
    private static double[] meanCurve(DoubleUnaryOperator lambdaFun, DoubleUnaryOperator serviceCcdf,
                                      double[] t, double startTime, double cut, int panels,
                                      boolean unbounded) {
        double[][] fixed = unbounded ? simpson(0.0, cut, panels) : null;
        double[] fixedGc = null;
        if (fixed != null) {
            fixedGc = new double[fixed[0].length];
            for (int j = 0; j < fixed[0].length; j++) {
                fixedGc[j] = serviceCcdf.applyAsDouble(fixed[0][j]);
            }
        }
        double[] m = new double[t.length];
        for (int i = 0; i < t.length; i++) {
            double[] x;
            double[] w;
            double[] gc;
            if (fixed != null) {
                x = fixed[0];
                w = fixed[1];
                gc = fixedGc;
            } else {
                double hi = Math.min(cut, Math.max(0.0, t[i] - startTime));
                double[][] xw = simpson(0.0, hi, panels);
                x = xw[0];
                w = xw[1];
                gc = new double[x.length];
                for (int j = 0; j < x.length; j++) {
                    gc[j] = serviceCcdf.applyAsDouble(x[j]);
                }
            }
            double acc = 0.0;
            for (int j = 0; j < x.length; j++) {
                // m(t) = int lambda(t-x) P(S>x) dx: arrivals of age x still in service.
                acc += w[j] * lambdaFun.applyAsDouble(t[i] - x[j]) * gc[j];
            }
            m[i] = acc;
        }
        return m;
    }

    /** Nodes and weights of the composite Simpson rule on an even panel count. */
    private static double[][] simpson(double a, double b, int n) {
        if (n % 2 == 1) {
            n++;
        }
        if (b <= a) {
            return new double[][]{{a}, {0.0}};
        }
        double[] x = new double[n + 1];
        double[] w = new double[n + 1];
        double h = (b - a) / n;
        for (int i = 0; i <= n; i++) {
            x[i] = a + i * h;
            w[i] = (i == 0 || i == n) ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0);
            w[i] *= h / 3.0;
        }
        return new double[][]{x, w};
    }

    /** Smallest doubling point at which the service ccdf is below tol. */
    private static double tailCut(DoubleUnaryOperator ccdf, double tol, double cap) {
        double x = 1.0;
        while (ccdf.applyAsDouble(x) > tol) {
            x *= 2.0;
            if (x > cap) {
                return cap;
            }
        }
        return x;
    }
}
