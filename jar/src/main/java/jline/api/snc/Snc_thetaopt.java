/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.snc;

import java.util.function.DoubleUnaryOperator;

/**
 * Minimizes a Chernoff bound over the free parameter theta.
 *
 * <p>Every bound in the snc domain holds for each theta &gt; 0 for which the
 * arrival MGF is finite and the station is stable, so the reported bound is the
 * infimum over theta. The objective is evaluated on a logarithmic grid,
 * non-finite values (a diverging MGF, an unstable leftover rate) are discarded,
 * and the best grid point is refined by golden-section search in
 * log10(theta).</p>
 *
 * <p>THE TWO-STAGE SEARCH IS NOT AN OPTIMIZATION OF CONVENIENCE: the feasible
 * set is an interval whose endpoints are not known in closed form once
 * envelopes are composed, and an unguarded local search steps into the
 * infeasible region and terminates there.</p>
 *
 * <p>Port of matlab/src/api/snc/snc_thetaopt.m, whose refinement is
 * {@code fminbnd} (golden section plus parabolic interpolation); the plain
 * golden section here reaches the same optimum on these smooth objectives.</p>
 */
public final class Snc_thetaopt {
    private Snc_thetaopt() {}

    /** Substituted for a non-finite objective, so the search can still compare it. */
    private static final double INFEASIBLE = 1e300;
    private static final int GRID = 600;
    private static final double INV_PHI = (Math.sqrt(5.0) - 1.0) / 2.0;

    /**
     * @param fun the objective, a function of theta
     * @return the minimum and its theta; (Infinity, NaN) if nothing is feasible
     */
    public static SncResult snc_thetaopt(DoubleUnaryOperator fun) {
        return snc_thetaopt(fun, 1e3);
    }

    /**
     * @param fun      the objective, a function of theta
     * @param thetamax upper end of the search range
     * @return the minimum and its theta; (Infinity, NaN) if nothing is feasible
     */
    public static SncResult snc_thetaopt(DoubleUnaryOperator fun, double thetamax) {
        if (thetamax <= 0) {
            throw new IllegalArgumentException("snc_thetaopt: thetamax must be positive, got " + thetamax);
        }
        double loExp = -6.0;
        double hiExp = Math.log10(thetamax);
        double[] grid = new double[GRID];
        double[] fval = new double[GRID];
        int imin = 0;
        for (int i = 0; i < GRID; i++) {
            grid[i] = Math.pow(10.0, loExp + (hiExp - loExp) * i / (GRID - 1.0));
            fval[i] = safeval(fun, grid[i]);
            if (fval[i] < fval[imin]) {
                imin = i;
            }
        }
        double val = fval[imin];
        if (val >= 1e299) {
            return new SncResult(Double.POSITIVE_INFINITY, Double.NaN);
        }
        double theta = grid[imin];

        double lo = Math.log10(grid[Math.max(imin - 1, 0)]);
        double hi = Math.log10(grid[Math.min(imin + 1, GRID - 1)]);
        if (hi > lo) {
            double x1 = hi - INV_PHI * (hi - lo);
            double x2 = lo + INV_PHI * (hi - lo);
            double f1 = safeval(fun, Math.pow(10.0, x1));
            double f2 = safeval(fun, Math.pow(10.0, x2));
            for (int it = 0; it < 200 && (hi - lo) > 1e-12; it++) {
                if (f1 < f2) {
                    hi = x2;
                    x2 = x1;
                    f2 = f1;
                    x1 = hi - INV_PHI * (hi - lo);
                    f1 = safeval(fun, Math.pow(10.0, x1));
                } else {
                    lo = x1;
                    x1 = x2;
                    f1 = f2;
                    x2 = lo + INV_PHI * (hi - lo);
                    f2 = safeval(fun, Math.pow(10.0, x2));
                }
            }
            double xopt = 0.5 * (lo + hi);
            double vopt = safeval(fun, Math.pow(10.0, xopt));
            if (vopt < val) {
                val = vopt;
                theta = Math.pow(10.0, xopt);
            }
        }
        return new SncResult(val, theta);
    }

    private static double safeval(DoubleUnaryOperator fun, double theta) {
        double v;
        try {
            v = fun.applyAsDouble(theta);
        } catch (RuntimeException e) {
            return INFEASIBLE;
        }
        if (Double.isNaN(v) || Double.isInfinite(v)) {
            return INFEASIBLE;
        }
        return v;
    }
}
