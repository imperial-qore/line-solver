/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.fluid.analyzers;

import java.util.List;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.fluid.handlers.FluidRateMultiplier;

import static jline.io.InputOutput.line_error;

/**
 * Time-varying per-state rate multiplier of the matrix-method fluid ODE.
 *
 * <p>The matrix method integrates {@code dx/dt = W' theta(x) + ALambda}, where
 * the drift out of an ODE state is linear in the service rate of the
 * (station, class) that owns it. A time-varying rate {@code rate(t)} for a
 * (station, class) is therefore exactly a multiplicative factor
 * {@code m(t) = rate(t)/nominal} applied to that state's entry of
 * {@code theta}, and to its nominal-rate throughput contribution. This is the
 * state-space counterpart of the per-event factor the closing ODE applies
 * through {@link FluidRateMultiplier} (MATLAB {@code solver_fluid_ratemult}),
 * and it lets the matrix method honour {@code options.config.rate_sched}
 * without switching solution method.</p>
 *
 * <p>Only the explicit {@code rate_sched} channel is expanded here. A schedule
 * on an EXT (Source) station modulates an open arrival stream, which enters the
 * matrix-method ODE through the constant {@code ALambda} inflow rather than
 * through {@code theta}; such an entry is rejected with a clear error pointing
 * at the closing ODE, which carries the arrival channel (and the NHPP path).</p>
 */
public class FluidStateRateMultiplier {

    private final double[] tgrid;
    private final double[][] mmat; // [ngrid][nStates]
    private final int nStates;

    private FluidStateRateMultiplier(double[] tgrid, double[][] mmat, int nStates) {
        this.tgrid = tgrid;
        this.mmat = mmat;
        this.nStates = nStates;
    }

    /**
     * Builds the multiplier for the kept ODE states, or returns null when no
     * schedule is configured, so the caller keeps the legacy autonomous ODE
     * numerically unchanged.
     *
     * @param options              solver options carrying {@code config.rate_sched}
     * @param keep                 indices of the retained states in the full state order
     * @param stationOfStateFull   station owning each full-order state
     * @param classOfStateFull     class owning each full-order state
     * @param sn                   network struct (for the nominal rates and the EXT check)
     */
    public static FluidStateRateMultiplier build(SolverOptions options, List<Integer> keep,
                                                 int[] stationOfStateFull, int[] classOfStateFull,
                                                 NetworkStruct sn) {
        if (options == null || options.config == null || options.config.rate_sched == null
                || options.config.rate_sched.isEmpty()) {
            return null;
        }
        List<FluidRateMultiplier.RateEntry> sched = options.config.rate_sched;

        // Union time grid of all entries, so every schedule is represented
        // without resampling loss.
        java.util.TreeSet<Double> ts = new java.util.TreeSet<Double>();
        for (int s = 0; s < sched.size(); s++) {
            FluidRateMultiplier.RateEntry e = sched.get(s);
            if (e.tgrid == null || e.rates == null || e.tgrid.length != e.rates.length
                    || e.tgrid.length == 0) {
                line_error(FluidStateRateMultiplier.class.getName(),
                        "rate_sched entry must carry tgrid and rates of the same non-zero length.");
            }
            for (int j = 0; j < e.tgrid.length; j++) {
                ts.add(e.tgrid[j]);
            }
        }
        double[] tgrid = new double[ts.size()];
        int p = 0;
        for (Double v : ts) {
            tgrid[p++] = v.doubleValue();
        }

        int nStates = keep.size();
        double[][] mmat = new double[tgrid.length][nStates];
        for (int j = 0; j < tgrid.length; j++) {
            for (int i = 0; i < nStates; i++) {
                mmat[j][i] = 1.0;
            }
        }

        boolean any = false;
        for (int s = 0; s < sched.size(); s++) {
            FluidRateMultiplier.RateEntry entry = sched.get(s);
            int ist = entry.station;
            int r = entry.jobclass;
            if (ist < 0 || ist >= sn.nstations || r < 0 || r >= sn.nclasses) {
                continue;
            }
            if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.EXT) {
                line_error(FluidStateRateMultiplier.class.getName(),
                        "rate_sched on a Source (EXT) station modulates an open arrival stream, "
                                + "which the matrix method carries as a constant inflow; "
                                + "use the 'closing' fluid method for that channel.");
            }
            double nominal = (entry.nominal != null) ? entry.nominal.doubleValue()
                    : sn.rates.get(ist, r);
            if (!(nominal > 0.0) && !(nominal < 0.0)) {
                continue; // no usable nominal rate: leave this channel at identity
            }
            for (int i = 0; i < nStates; i++) {
                int ki = keep.get(i).intValue();
                if (stationOfStateFull[ki] != ist || classOfStateFull[ki] != r) {
                    continue;
                }
                any = true;
                for (int j = 0; j < tgrid.length; j++) {
                    mmat[j][i] = interpClamped(entry.tgrid, entry.rates, tgrid[j]) / nominal;
                }
            }
        }
        if (!any) {
            return null;
        }
        return new FluidStateRateMultiplier(tgrid, mmat, nStates);
    }

    /** Per-state multiplier at time {@code t}, by clamped linear interpolation. */
    public double[] multAt(double t) {
        double[] out = new double[nStates];
        if (t <= tgrid[0]) {
            System.arraycopy(mmat[0], 0, out, 0, nStates);
            return out;
        }
        if (t >= tgrid[tgrid.length - 1]) {
            System.arraycopy(mmat[tgrid.length - 1], 0, out, 0, nStates);
            return out;
        }
        int lo = 0;
        int hi = tgrid.length - 1;
        while (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            if (tgrid[mid] <= t) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        double w = (t - tgrid[lo]) / (tgrid[lo + 1] - tgrid[lo]);
        for (int i = 0; i < nStates; i++) {
            out[i] = (1.0 - w) * mmat[lo][i] + w * mmat[lo + 1][i];
        }
        return out;
    }

    /** Clamped piecewise-linear interpolation of {@code (tg, val)} at {@code tt}. */
    private static double interpClamped(double[] tg, double[] val, double tt) {
        int n = tg.length;
        if (n == 1 || tt <= tg[0]) {
            return val[0];
        }
        if (tt >= tg[n - 1]) {
            return val[n - 1];
        }
        int lo = 0;
        int hi = n - 1;
        while (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            if (tg[mid] <= tt) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        double w = (tt - tg[lo]) / (tg[lo + 1] - tg[lo]);
        return (1.0 - w) * val[lo] + w * val[lo + 1];
    }
}
