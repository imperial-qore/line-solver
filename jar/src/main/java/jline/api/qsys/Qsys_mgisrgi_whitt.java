/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.qsys;

import java.util.ArrayList;
import java.util.List;
import java.util.function.UnaryOperator;

import org.apache.commons.math3.complex.Complex;

import jline.api.lti.Laplace_invert;

/**
 * Engineering solution of the call-center model M/GI/s/r+GI.
 *
 * <p>Poisson arrivals at rate lambda, iid general service times of mean 1/mu, s
 * servers, r extra waiting spaces and iid patience times with a general
 * distribution.
 *
 * <p>TWO APPROXIMATIONS. The general patience law is replaced by STATE-DEPENDENT
 * Markovian abandonment: a customer who is jth from the end of a queue abandons
 * at rate delta_j = h(j/lambda) for the patience hazard h = f/(1-F) (eq. 3.3),
 * because such a customer has been waiting for about j/lambda; the total
 * abandonment rate with k waiting is Delta_k = sum_{j&lt;=k} delta_j (eq. 3.4).
 * The general service law is replaced by an exponential of the same mean
 * (Section 5), accurate here because with many servers and non-negligible
 * abandonment the model behaves like a loss system, where the service law is
 * insensitive beyond its mean. What is left is the Markovian M/M/s/r+M(n) model,
 * solved as a birth-and-death process.
 *
 * <p>WHAT THE PATIENCE LAW CONTRIBUTES is only its hazard NEAR THE ORIGIN, not
 * its mean and not its tail: waits are O(1/sqrt(s)) in the many-server regime,
 * so a customer either abandons early or never.
 *
 * <p>Exact for M/M/s/r+M, which is {@link Qsys_erlanga}.
 *
 * <p>Port of MATLAB qsys_mgisrgi_whitt.m.
 *
 * <p>Reference: W. Whitt (2005). Engineering solution of a basic call-center
 * model. Management Science 51(2), 221-235.
 *
 * @since LINE 3.1.0
 */
public final class Qsys_mgisrgi_whitt {

    private Qsys_mgisrgi_whitt() {
    }

    /** Truncation level used when the waiting room is infinite. */
    public static final int DEFAULT_MAX_QUEUE = 100000;
    /** Relative tail tolerance for that truncation. */
    public static final double DEFAULT_TOL = 1e-14;

    /**
     * Steady-state measures with no waiting-time cdfs and the default truncation.
     *
     * @param lambda   arrival rate
     * @param mu       service rate of one server, the reciprocal of the mean service time
     * @param s        number of servers
     * @param r        extra waiting spaces, {@code Double.POSITIVE_INFINITY} if unbounded
     * @param patience the patience law
     * @return the steady-state measures
     */
    public static QsysAbandonResult qsys_mgisrgi_whitt(double lambda, double mu, int s, double r,
                                                       Patience patience) {
        return qsys_mgisrgi_whitt(lambda, mu, s, r, patience, null, DEFAULT_MAX_QUEUE, DEFAULT_TOL,
                "euler", 41);
    }

    /**
     * Steady-state measures, with the waiting-time cdfs when times are supplied.
     *
     * @param lambda    arrival rate
     * @param mu        service rate of one server
     * @param s         number of servers
     * @param r         extra waiting spaces, {@code Double.POSITIVE_INFINITY} if unbounded
     * @param patience  the patience law
     * @param wPoints   times at which to evaluate the waiting-time cdfs, or null
     * @param maxQueue  truncation level used when r is infinite
     * @param tol       relative tail tolerance for that truncation
     * @param invMethod Laplace inversion method for the cdfs
     * @param invN      number of inversion nodes
     * @return the steady-state measures
     */
    public static QsysAbandonResult qsys_mgisrgi_whitt(double lambda, double mu, int s, double r,
                                                       Patience patience, double[] wPoints,
                                                       int maxQueue, double tol, String invMethod,
                                                       int invN) {
        if (lambda <= 0) {
            throw new RuntimeException("qsys_mgisrgi_whitt: the arrival rate lambda must be positive");
        }
        if (mu <= 0) {
            throw new RuntimeException("qsys_mgisrgi_whitt: the service rate mu must be positive");
        }
        if (s < 1) {
            throw new RuntimeException("qsys_mgisrgi_whitt: the number of servers s must be at least 1");
        }
        if (r < 0) {
            throw new RuntimeException(
                    "qsys_mgisrgi_whitt: the number of extra waiting spaces r must be non-negative");
        }
        if (patience == null) {
            throw new RuntimeException("qsys_mgisrgi_whitt: the patience law must not be null");
        }

        boolean finiteR = !Double.isInfinite(r);
        int rr = finiteR ? (int) Math.round(r) : maxQueue;

        // The birth-death recursion of eqs. (7.4)-(7.7), unnormalized with x_s = 1.
        double[] xUp = new double[rr + 1];
        double[] dlt = new double[rr + 1];
        double[] delta = new double[Math.max(rr, 1)];
        xUp[0] = 1.0;
        double smu = s * mu;
        int kUsed = rr;
        double peak = 1.0;
        for (int k = 0; k < rr; k++) {
            int j = k + 1;
            double[] step = rateStep(j, lambda, dlt[j - 1], patience);
            delta[j - 1] = step[0];
            dlt[j] = step[1];
            xUp[k + 1] = lambda * xUp[k] / (smu + dlt[j]);
            peak = Math.max(peak, xUp[k + 1]);
            if (!finiteR && xUp[k + 1] < tol * peak && k >= 1) {
                kUsed = j;
                break;
            }
        }
        if (!finiteR) {
            if (kUsed == rr && rr > 0) {
                throw new RuntimeException("qsys_mgisrgi_whitt: the queue-length tail is still "
                        + (xUp[rr] / peak) + " of its peak at the truncation level " + rr
                        + "; with r = Inf the patience law must make the chain ergodic");
            }
            double[] xUpCut = new double[kUsed + 1];
            double[] dltCut = new double[kUsed + 1];
            double[] deltaCut = new double[kUsed];
            System.arraycopy(xUp, 0, xUpCut, 0, kUsed + 1);
            System.arraycopy(dlt, 0, dltCut, 0, kUsed + 1);
            System.arraycopy(delta, 0, deltaCut, 0, kUsed);
            xUp = xUpCut;
            dlt = dltCut;
            delta = deltaCut;
            rr = kUsed;
        }

        // The downward leg, eq. (7.5), over the states where not all servers are busy.
        double[] xDown = new double[s];
        double xk = 1.0;
        for (int k = s; k >= 1; k--) {
            xk = k * mu * xk / lambda;
            xDown[k - 1] = xk;
        }

        double[] x = new double[s + rr + 1];
        System.arraycopy(xDown, 0, x, 0, s);
        System.arraycopy(xUp, 0, x, s, rr + 1);
        double total = 0.0;
        for (int i = 0; i < x.length; i++) {
            total += x[i];
        }
        double[] p = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            p[i] = x[i] / total;
        }
        double probLoss = finiteR ? p[p.length - 1] : 0.0;
        double[] pa = new double[p.length];
        for (int i = 0; i < p.length; i++) {
            pa[i] = p[i] / (1.0 - probLoss);       // eq. (7.8), the state seen by an ENTERING customer
        }

        double meanNumber = 0.0;
        double meanQueue = 0.0;
        double busy = 0.0;
        for (int k = 0; k < p.length; k++) {
            meanNumber += k * p[k];
            meanQueue += Math.max(0, k - s) * p[k];
            busy += Math.min(k, s) * p[k];
        }
        double varNumber = 0.0;
        double varQueue = 0.0;
        for (int k = 0; k < p.length; k++) {
            varNumber += (k - meanNumber) * (k - meanNumber) * p[k];
            double q = Math.max(0, k - s);
            varQueue += (q - meanQueue) * (q - meanQueue) * p[k];
        }

        double probNoWait = 0.0;                    // eq. (7.9), states 0..s-1
        for (int k = 0; k < s; k++) {
            probNoWait += pa[k];
        }

        double[] sigma = new double[rr];
        double[] mSum = new double[rr];
        double[] vSum = new double[rr];
        double[] ewa1 = new double[rr];
        double[] ewa2 = new double[rr];
        List<double[][]> kernels = new ArrayList<double[][]>(rr);
        for (int k = 1; k <= rr; k++) {
            double[][] kern = kernel(k, smu, dlt, delta);
            double[] rate = kern[0];
            double[] phik = kern[1];
            double[] surv = kern[2];
            double prod = 1.0;
            double sm = 0.0;
            double sv = 0.0;
            double cumM = 0.0;
            double cumV = 0.0;
            double e1 = 0.0;
            double e2 = 0.0;
            for (int j = 0; j < k; j++) {
                double m = 1.0 / rate[j];
                prod *= (1.0 - phik[j]);
                sm += m;
                sv += m * m;
                // Eqs. (7.28)-(7.29): abandoning at the jth departure epoch costs
                // the sum of the first j interdeparture times.
                cumM += m;
                cumV += m * m;
                double w = surv[j] * phik[j];
                e1 += w * cumM;
                e2 += w * (cumV + cumM * cumM);
            }
            sigma[k - 1] = prod;
            mSum[k - 1] = sm;
            vSum[k - 1] = sv;
            ewa1[k - 1] = e1;
            ewa2[k - 1] = e2;
            kernels.add(kern);
        }

        // Finding s+k in system puts the arrival in position k+1, so the weights
        // are pa_{s+k} for k = 0..r-1.
        double[] wArr = new double[rr];
        System.arraycopy(pa, s, wArr, 0, rr);
        double probServed = probNoWait;
        double ews1 = 0.0;
        double ews2 = 0.0;
        double ewa1Tot = 0.0;
        double ewa2Tot = 0.0;
        for (int k = 0; k < rr; k++) {
            probServed += wArr[k] * sigma[k];
            ews1 += wArr[k] * sigma[k] * mSum[k];                          // eq. (7.16)
            ews2 += wArr[k] * sigma[k] * (vSum[k] + mSum[k] * mSum[k]);    // eq. (7.17)
            ewa1Tot += wArr[k] * ewa1[k];                                  // eq. (7.26)
            ewa2Tot += wArr[k] * ewa2[k];                                  // eq. (7.27)
        }
        double probAbandon = 1.0 - probServed;
        double meanWaitServed = ratio(ews1, probServed);
        double meanWaitAbandon = ratio(ewa1Tot, probAbandon);

        double[] waitPoints = null;
        double[] cdfWaitServed = null;
        double[] cdfWaitAbandon = null;
        double[] cdfWait = null;
        if (wPoints != null && wPoints.length > 0) {
            waitPoints = wPoints.clone();
            cdfWaitServed = new double[wPoints.length];
            cdfWaitAbandon = new double[wPoints.length];
            cdfWait = new double[wPoints.length];
            final double[] fwArr = wArr;
            final double[] fsigma = sigma;
            final List<double[][]> fkernels = kernels;
            UnaryOperator<Complex> served = new UnaryOperator<Complex>() {
                @Override
                public Complex apply(Complex z) {
                    return transform(z, fwArr, fsigma, fkernels, true).divide(z);
                }
            };
            UnaryOperator<Complex> abandoned = new UnaryOperator<Complex>() {
                @Override
                public Complex apply(Complex z) {
                    return transform(z, fwArr, fsigma, fkernels, false).divide(z);
                }
            };
            double capS = Math.max(probServed - probNoWait, 0.0);
            double capA = Math.max(probAbandon, 0.0);
            for (int i = 0; i < wPoints.length; i++) {
                double fs = Laplace_invert.laplace_invert(served, wPoints[i], invMethod, invN);
                double fa = Laplace_invert.laplace_invert(abandoned, wPoints[i], invMethod, invN);
                fs = Math.min(Math.max(fs, 0.0), capS);
                fa = Math.min(Math.max(fa, 0.0), capA);
                cdfWaitServed[i] = (probNoWait + fs) / Math.max(probServed, Double.MIN_NORMAL);
                cdfWaitAbandon[i] = fa / Math.max(probAbandon, Double.MIN_NORMAL);
                cdfWait[i] = probNoWait + fs + fa;
            }
        }

        return new QsysAbandonResult(p, probLoss, probNoWait, probServed, probAbandon, meanNumber,
                varNumber, meanQueue, varQueue, busy / s,
                lambda * (1.0 - probLoss) * probServed, lambda * (1.0 - probLoss) * probAbandon,
                meanWaitServed,
                Math.max(0.0, ratio(ews2, probServed) - meanWaitServed * meanWaitServed),
                meanWaitAbandon,
                Math.max(0.0, ratio(ewa2Tot, probAbandon) - meanWaitAbandon * meanWaitAbandon),
                ews1 + ewa1Tot, ews2 + ewa2Tot, delta, dlt, rr,
                patience.getForm() == Patience.Form.EXPONENTIAL, patience.getTheta(),
                waitPoints, cdfWaitServed, cdfWaitAbandon, cdfWait);
    }

    /**
     * One step of eqs. (3.3)-(3.4) (hazard form) or (3.5)-(3.6) (ccdf form).
     *
     * <p>DIVERGENCE from the printed eqs. (3.5)-(3.6): they read
     * delta_j = int_{(j-1)/lambda}^{j/lambda} h(t) dt and
     * Delta_k = -log F^c(k/lambda), which are cumulative hazards, i.e.
     * dimensionless, while delta and Delta are rates everywhere else in the
     * paper. They are the AVERAGE hazard over an interval of length 1/lambda, so
     * the factor lambda is missing. Restoring it makes the ccdf form reduce to
     * the exact Erlang A rates under exponential patience, which the paper
     * states this approximation does (eq. 7.12).
     */
    private static double[] rateStep(int j, double lambda, double deltaPrev, Patience patience) {
        double deltaJ;
        double deltaTot;
        if (patience.getForm() == Patience.Form.CCDF) {
            double g = patience.ccdfAt(j / lambda);
            if (g <= 0) {
                throw new RuntimeException("qsys_mgisrgi_whitt: the patience ccdf vanishes at t = "
                        + (j / lambda) + ", so every customer has abandoned by then; supply a "
                        + "hazard instead");
            }
            deltaTot = -lambda * Math.log(g);
            deltaJ = deltaTot - deltaPrev;
        } else {
            deltaJ = patience.hazardAt(j / lambda);
            deltaTot = deltaPrev + deltaJ;
        }
        if (deltaJ < 0) {
            throw new RuntimeException(
                    "qsys_mgisrgi_whitt: the patience law produced a negative abandonment rate");
        }
        return new double[]{deltaJ, deltaTot};
    }

    /**
     * Eqs. (7.10)-(7.11): with k waiting, the total departure rate before the jth
     * departure epoch is s*mu + Delta_k - Delta_{j-1}, of which delta_j is the
     * share belonging to the customer of interest. Returns {rates, phi, survival},
     * the survival being prod_{l&lt;j}(1-phi_k(l)).
     */
    private static double[][] kernel(int k, double smu, double[] dlt, double[] delta) {
        double[] rate = new double[k];
        double[] phi = new double[k];
        double[] surv = new double[k];
        double running = 1.0;
        for (int j = 1; j <= k; j++) {
            rate[j - 1] = smu + dlt[k] - dlt[j - 1];
            phi[j - 1] = delta[j - 1] / rate[j - 1];
            surv[j - 1] = running;
            running *= (1.0 - phi[j - 1]);
        }
        return new double[][]{rate, phi, surv};
    }

    /** A conditional moment is 0/0 when the conditioning event cannot happen. */
    private static double ratio(double num, double den) {
        return den <= 0 ? 0.0 : num / den;
    }

    /**
     * Eqs. (7.22)-(7.23) when served, eqs. (7.32)-(7.33) otherwise. Both fold the
     * same per-position kernel: the wait is a sum of exponentials with rates
     * 1/m_k(j), truncated at the departure epoch that serves or loses the customer.
     */
    private static Complex transform(Complex z, double[] wArr, double[] sigma,
                                     List<double[][]> kernels, boolean served) {
        Complex val = Complex.ZERO;
        for (int k = 1; k <= wArr.length; k++) {
            double[][] kern = kernels.get(k - 1);
            double[] rate = kern[0];
            double[] phi = kern[1];
            double[] surv = kern[2];
            Complex chain = Complex.ONE;
            if (served) {
                for (int j = 0; j < rate.length; j++) {
                    chain = chain.multiply(z.add(rate[j]).reciprocal().multiply(rate[j]));
                }
                val = val.add(chain.multiply(wArr[k - 1] * sigma[k - 1]));
            } else {
                for (int j = 0; j < rate.length; j++) {
                    chain = chain.multiply(z.add(rate[j]).reciprocal().multiply(rate[j]));
                    val = val.add(chain.multiply(wArr[k - 1] * surv[j] * phi[j]));
                }
            }
        }
        return val;
    }
}
