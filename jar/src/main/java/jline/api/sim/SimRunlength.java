/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.sim;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

import jline.util.matrix.Matrix;

/**
 * Run-length planning for steady-state simulation.
 *
 * <p>THE QUANTITY THAT MATTERS is not the variance of the process but its
 * ASYMPTOTIC VARIANCE sigma^2 = lim t Var(time-average over [0,t]), twice the
 * integral of the autocovariance: a time average of a positively correlated
 * process converges at rate sigma^2/t, not Var(X)/t. Then
 * t* = (z/eps)^2 sigma^2/mean^2 is the run needed for relative precision eps.
 *
 * <p>For M/M/1, sigma^2 = 2 rho(1+rho)/(mu (1-rho)^4) in closed form. The FOURTH
 * power is the whole story; divided by the squared mean it leaves a run length
 * growing like (1-rho)^-2, so a queue at rho = 0.9 needs about 100 times the run
 * of one at rho = 0.
 *
 * <p>Port of MATLAB sim_runlength.m, sim_asymvar_mm1.m and sim_asymvar_ctmc.m.
 *
 * <p>Reference: W. Whitt (1989). Planning queueing simulations. Management
 * Science 35(11), 1341-1366.
 *
 * @since LINE 3.1.0
 */
public final class SimRunlength {

    private SimRunlength() {
    }

    /**
     * Asymptotic variance of the M/M/1 number-in-system process.
     *
     * @param lambda arrival rate
     * @param mu     service rate
     * @return map with mean, variance, asymptoticVariance and relaxationTime
     */
    public static Map<String, Double> sim_asymvar_mm1(double lambda, double mu) {
        if (lambda <= 0 || mu <= 0) {
            throw new RuntimeException("sim_asymvar_mm1: the rates must be positive");
        }
        double rho = lambda / mu;
        if (rho >= 1) {
            throw new RuntimeException("sim_asymvar_mm1: the queue must be stable, rho < 1");
        }
        double mean = rho / (1.0 - rho);
        double var = rho / ((1.0 - rho) * (1.0 - rho));
        double d = Math.pow(1.0 - rho, 4);
        double asym = 2.0 * rho * (1.0 + rho) / (mu * d);
        Map<String, Double> res = new HashMap<String, Double>();
        res.put("mean", mean);
        res.put("variance", var);
        res.put("asymptoticVariance", asym);
        res.put("relaxationTime", var > 0 ? asym / var : 0.0);
        return res;
    }

    /**
     * Asymptotic variance of a reward on a CTMC: 2 sum_x pi(x)g(x)d(x) with
     * g = f - E_pi[f] and A d = -g.
     *
     * <p>A d = -g pins d only up to a constant, so one equation of A is
     * redundant and one normalization replaces it. WHICH equation is dropped
     * matters: the rows of A are related by pi A = 0, so a row whose pi is tiny
     * is only nominally redundant, and dropping it loses real information -- on
     * a queue truncated where pi has underflowed, that alone puts sigma^2 out by
     * orders of magnitude. Dropping the row with the LARGEST pi is the
     * well-conditioned choice.
     *
     * @param A the generator, rows summing to zero
     * @param f the reward attached to each state
     * @return map with mean, variance and asymptoticVariance
     */
    public static Map<String, Double> sim_asymvar_ctmc(Matrix A, double[] f) {
        int n = A.getNumRows();
        if (A.getNumCols() != n) {
            throw new RuntimeException("sim_asymvar_ctmc: the generator must be square");
        }
        if (f.length != n) {
            throw new RuntimeException("sim_asymvar_ctmc: one reward per state is required");
        }
        for (int i = 0; i < n; i++) {
            double row = 0.0;
            for (int j = 0; j < n; j++) {
                row += A.get(i, j);
            }
            if (Math.abs(row) > 1e-8) {
                throw new RuntimeException("sim_asymvar_ctmc: the generator rows must sum to zero");
            }
        }
        // pi A = 0 with sum pi = 1.
        Matrix M = new Matrix(n, n);
        Matrix b = new Matrix(n, 1);
        for (int j = 0; j < n - 1; j++) {
            for (int i = 0; i < n; i++) {
                M.set(j, i, A.get(i, j));
            }
        }
        for (int i = 0; i < n; i++) {
            M.set(n - 1, i, 1.0);
        }
        b.set(n - 1, 0, 1.0);
        Matrix piM = new Matrix(n, 1);
        Matrix.solve(M, b, piM);
        double[] pi = new double[n];
        int drop = 0;
        for (int i = 0; i < n; i++) {
            pi[i] = piM.get(i, 0);
            if (pi[i] > pi[drop]) {
                drop = i;
            }
        }
        double mean = 0.0;
        for (int i = 0; i < n; i++) {
            mean += pi[i] * f[i];
        }
        double[] g = new double[n];
        for (int i = 0; i < n; i++) {
            g[i] = f[i] - mean;
        }
        Matrix M2 = new Matrix(n, n);
        Matrix b2 = new Matrix(n, 1);
        int r0 = 0;
        for (int i = 0; i < n; i++) {
            if (i == drop) {
                continue;
            }
            for (int j = 0; j < n; j++) {
                M2.set(r0, j, A.get(i, j));
            }
            b2.set(r0, 0, -g[i]);
            r0++;
        }
        for (int j = 0; j < n; j++) {
            M2.set(n - 1, j, pi[j]);
        }
        Matrix dM = new Matrix(n, 1);
        Matrix.solve(M2, b2, dM);
        double var = 0.0;
        double asym = 0.0;
        for (int i = 0; i < n; i++) {
            var += pi[i] * g[i] * g[i];
            asym += pi[i] * g[i] * dM.get(i, 0);
        }
        Map<String, Double> res = new HashMap<String, Double>();
        res.put("mean", mean);
        res.put("variance", var);
        res.put("asymptoticVariance", 2.0 * asym);
        res.put("relaxationTime", var > 0 ? 2.0 * asym / var : 0.0);
        return res;
    }

    /**
     * Run length for a steady-state estimate of a given relative precision.
     *
     * @param mean         the steady-state mean being estimated
     * @param asymVar      sigma^2 of that estimator
     * @param relPrecision the target half-width as a fraction of the mean
     * @param confidence   the confidence level
     * @param runLength    an actual run length, or 0 to skip the report
     * @return map with requiredRunLength, z and, when runLength is positive,
     *         halfWidth and achievedRelPrecision
     */
    public static Map<String, Double> sim_runlength(double mean, double asymVar,
                                                     double relPrecision, double confidence,
                                                     double runLength) {
        if (mean == 0) {
            throw new RuntimeException(
                    "sim_runlength: a relative precision is meaningless for a zero mean");
        }
        if (asymVar < 0) {
            throw new RuntimeException("sim_runlength: the asymptotic variance cannot be negative");
        }
        if (relPrecision <= 0) {
            throw new RuntimeException("sim_runlength: the relative precision must be positive");
        }
        if (confidence <= 0 || confidence >= 1) {
            throw new RuntimeException("sim_runlength: the confidence must lie in (0,1)");
        }
        double lo = 0.0;
        double hi = 40.0;
        double target = 1.0 - confidence;
        for (int i = 0; i < 200; i++) {
            double mid = 0.5 * (lo + hi);
            if (org.apache.commons.math3.special.Erf.erfc(mid / Math.sqrt(2.0)) > target) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        double z = 0.5 * (lo + hi);
        Map<String, Double> res = new HashMap<String, Double>();
        res.put("z", z);
        res.put("requiredRunLength", (z / relPrecision) * (z / relPrecision) * asymVar / (mean * mean));
        if (runLength > 0) {
            double hw = z * Math.sqrt(asymVar / runLength);
            res.put("halfWidth", hw);
            res.put("achievedRelPrecision", hw / Math.abs(mean));
        }
        return res;
    }

    /**
     * How long a simulation run should have been, from the one it already did.
     *
     * <p>A batch-means half-width H at confidence 1-alpha over a run of N
     * samples pins the ASYMPTOTIC variance of the estimator,
     * sigma^2 = (H/z)^2 N, and that is the quantity a run length is planned from
     * -- NOT the stationary variance, which on M/M/1 differs from it by a factor
     * blowing up like (1-rho)^-2. {@link #sim_runlength} then turns it into the
     * sample count that reaches a requested RELATIVE precision.
     *
     * <p>An entry with a non-positive mean or half-width is left NaN, since
     * there is nothing to plan from there.
     *
     * @param means        one mean per (station, class)
     * @param ciHalfWidth  the confidence-interval half-width of the same entries
     * @param samplesUsed  the run length those half-widths came from
     * @param relPrecision the relative precision to plan for
     * @param confidence   the level the half-widths were computed at
     * @return map with asymptoticVariance and requiredSamples (matrices), and
     *         relprecision, confidence and samplesUsed (scalars)
     */
    public static Map<String, Object> sim_runlength_plan(Matrix means, Matrix ciHalfWidth,
                                                          double samplesUsed, double relPrecision,
                                                          double confidence) {
        if (samplesUsed <= 0) {
            throw new RuntimeException(
                    "sim_runlength_plan: the number of samples already used must be positive");
        }
        int M = means.getNumRows();
        int K = means.getNumCols();
        Matrix asym = new Matrix(M, K);
        Matrix req = new Matrix(M, K);
        // The same z the interval itself was built with; sim_runlength above
        // computes it inline from erfc, so it is read back from there rather
        // than recomputed by a second rule.
        double z = sim_runlength(1.0, 0.0, 1.0, confidence, 0.0).get("z");
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < K; r++) {
                asym.set(i, r, Double.NaN);
                req.set(i, r, Double.NaN);
                double h = (ciHalfWidth.getNumRows() > i && ciHalfWidth.getNumCols() > r)
                        ? ciHalfWidth.get(i, r) : Double.NaN;
                double m = means.get(i, r);
                if (!Double.isFinite(h) || h <= 0 || !Double.isFinite(m) || m <= 0) {
                    continue;
                }
                double av = (h / z) * (h / z) * samplesUsed;
                asym.set(i, r, av);
                req.set(i, r, sim_runlength(m, av, relPrecision, confidence, 0.0)
                        .get("requiredRunLength"));
            }
        }
        Map<String, Object> out = new LinkedHashMap<String, Object>();
        out.put("relprecision", relPrecision);
        out.put("confidence", confidence);
        out.put("samplesUsed", samplesUsed);
        out.put("asymptoticVariance", asym);
        out.put("requiredSamples", req);
        return out;
    }
}
