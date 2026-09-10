/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.fluid;

/**
 * Clamped piecewise-linear interpolation of a vector-valued trajectory.
 *
 * <p>This is the shared time-varying-input evaluator of the fluid solver, used
 * both by the trajectory-based iteration (TBI) transient, which reads the
 * frozen complement trajectory of a cell, and by the time-varying rate
 * multiplier of the closing ODE ({@link jline.solvers.fluid.handlers.FluidRateMultiplier}).
 * It mirrors MATLAB {@code fluid_interpcols.m}, transposed to the row-major
 * (one row per time point) layout used throughout the Java fluid code: MATLAB
 * samples a quantity in the columns of a matrix, here each row {@code y[j]} is
 * the sample at time {@code t[j]}.</p>
 *
 * <p>Times outside {@code [t[0], t[n-1]]} are clamped to the boundary rows
 * (zero-order hold outside the grid).</p>
 */
public final class FluidInterp {

    private FluidInterp() {
    }

    /**
     * Interpolates the trajectory {@code (t, y)} at time {@code tq}.
     *
     * @param t   strictly increasing time grid, length n
     * @param y   n rows of dimension at least {@code dim}, row j sampled at t[j]
     * @param tq  query time
     * @param dim number of components to return
     * @return the interpolated (clamped) row, of length {@code dim}
     */
    public static double[] interp(double[] t, double[][] y, double tq, int dim) {
        int n = t.length;
        double[] out = new double[dim];
        if (n == 1 || tq <= t[0]) {
            for (int j = 0; j < dim; j++) {
                out[j] = y[0][j];
            }
            return out;
        }
        if (tq >= t[n - 1]) {
            for (int j = 0; j < dim; j++) {
                out[j] = y[n - 1][j];
            }
            return out;
        }
        int lo = 0;
        int hi = n - 1;
        while (hi - lo > 1) {
            int mid = (lo + hi) >>> 1;
            if (t[mid] <= tq) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        double denom = t[hi] - t[lo];
        double frac = (denom > 0.0) ? (tq - t[lo]) / denom : 0.0;
        for (int j = 0; j < dim; j++) {
            out[j] = y[lo][j] + frac * (y[hi][j] - y[lo][j]);
        }
        return out;
    }
}
