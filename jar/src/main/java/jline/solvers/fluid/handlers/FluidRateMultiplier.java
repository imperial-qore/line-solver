/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.fluid.handlers;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import jline.lang.JobClass;
import jline.lang.nodes.Station;
import jline.lang.processes.NHPP;
import jline.solvers.SolverOptions;
import jline.solvers.fluid.FluidInterp;
import jline.util.matrix.Matrix;

import static jline.io.InputOutput.line_error;

/**
 * Time-varying per-event rate multiplier of the closing fluid ODE.
 *
 * <p>Mirrors MATLAB {@code solver_fluid_ratemult.m}. The multiplier is a
 * trajectory {@code m(t)} of length {@code numEvents}: event {@code e} of the
 * closing machinery is scaled at time {@code t} by {@code m(t)[e]}, evaluated
 * by clamped piecewise-linear interpolation ({@link FluidInterp}). The closing
 * rate is {@code rate = rateBase .* theta(x)} with {@code rateBase} linear in
 * the station-class service/arrival rate, so every time-varying source below
 * reduces to such a per-event multiplicative factor.</p>
 *
 * <p>Three independent, composable channels are honoured, all read from
 * {@link SolverOptions.Config} under the same names as the MATLAB options:</p>
 * <ol>
 *   <li>{@code rate_traj_tgrid} / {@code rate_traj_mmat}: a caller-supplied
 *       event multiplier matrix ({@code numEvents x ngrid}), used by the
 *       coupled LN layer transient.</li>
 *   <li>{@code nhpp_sched}: non-homogeneous Poisson source intensities. The
 *       nominal (time-average) rate baked into {@code rateBase} for the
 *       station-class is {@code mu(i,c)(0)}; the multiplier is
 *       {@code getRateAt(t)/nominal}, applied to every event sourced at that
 *       station and class.</li>
 *   <li>{@code rate_sched}: explicit per-(station,class) rate trajectories,
 *       reusing the same station-class to event expansion as the NHPP path.</li>
 * </ol>
 *
 * <p>{@link #build} returns null when no channel is configured, so the caller
 * keeps the legacy autonomous closure numerically unchanged.</p>
 */
public class FluidRateMultiplier {

    private final double[] tgrid;
    private final double[][] mmat; // [ngrid][numEvents]
    private final int numEvents;

    private FluidRateMultiplier(double[] tgrid, double[][] mmat, int numEvents) {
        this.tgrid = tgrid;
        this.mmat = mmat;
        this.numEvents = numEvents;
    }

    /**
     * Builds the multiplier for the closing machinery described by the given
     * event mapping, or returns null when no time-varying channel is set.
     *
     * @param numEvents  number of events of the (possibly immediate-eliminated) closing ODE
     * @param enabled    per-(station,class) enabled flags
     * @param qIndices   per-(station,class) starting state index
     * @param kic        per-(station,class) phase count
     * @param mu         per-(station,class) phase rate vectors
     * @param stations   station list of the sn station space
     * @param jobclasses job class list of the sn class space
     * @param eventIdx   per-event state index (numEvents x 1)
     * @param options    solver options carrying the config channels and timespan
     */
    public static FluidRateMultiplier build(int numEvents,
                                            boolean[][] enabled,
                                            Matrix qIndices,
                                            Matrix kic,
                                            Map<Station, Map<JobClass, Matrix>> mu,
                                            List<Station> stations,
                                            List<JobClass> jobclasses,
                                            Matrix eventIdx,
                                            SolverOptions options) {
        if (options == null || options.config == null || numEvents <= 0) {
            return null;
        }
        SolverOptions.Config cfg = options.config;

        // -- channel (1): caller-supplied rate_traj -------------------------------
        Traj user = null;
        if (cfg.rate_traj_tgrid != null && cfg.rate_traj_tgrid.length > 0 && cfg.rate_traj_mmat != null) {
            Matrix mmatIn = cfg.rate_traj_mmat;
            if (mmatIn.getNumRows() != numEvents) {
                line_error(FluidRateMultiplier.class.getName(),
                        String.format("rate_traj multiplier matrix has %d rows but the closing ODE has %d events.",
                                mmatIn.getNumRows(), numEvents));
            }
            int ngrid = cfg.rate_traj_tgrid.length;
            if (mmatIn.getNumCols() != ngrid) {
                line_error(FluidRateMultiplier.class.getName(),
                        String.format("rate_traj multiplier matrix has %d columns but the time grid has %d points.",
                                mmatIn.getNumCols(), ngrid));
            }
            double[][] mm = new double[ngrid][numEvents];
            for (int j = 0; j < ngrid; j++) {
                for (int e = 0; e < numEvents; e++) {
                    mm[j][e] = mmatIn.get(e, j);
                }
            }
            user = new Traj(cfg.rate_traj_tgrid.clone(), mm);
        }

        // horizon over which to expand a (possibly cyclic) schedule
        double t0 = 0.0;
        double tend = Double.POSITIVE_INFINITY;
        if (options.timespan != null && options.timespan.length >= 2) {
            if (Double.isFinite(options.timespan[0])) {
                t0 = options.timespan[0];
            }
            tend = options.timespan[1];
        }

        // -- channel (2): NHPP source intensities ---------------------------------
        Traj nhppTraj = null;
        if (cfg.nhpp_sched != null && !cfg.nhpp_sched.isEmpty()) {
            for (int s = 0; s < cfg.nhpp_sched.size(); s++) {
                NhppEntry entry = cfg.nhpp_sched.get(s);
                int i = entry.station;
                int c = entry.jobclass;
                if (i < 0 || c < 0 || i >= enabled.length || c >= enabled[i].length || !enabled[i][c]) {
                    continue; // class not served/active at this station in the fluid ODE
                }
                double nominal = nominalRate(mu, stations, jobclasses, i, c);
                if (!(nominal > 0.0)) {
                    continue;
                }
                NHPP nh = entry.nhpp;
                // horizon: for a non-finite timespan use a few periods so a cyclic
                // schedule is represented rather than clamped after one segment
                double period = nh.getPeriod();
                double thi;
                if (!Double.isFinite(tend)) {
                    if (Double.isFinite(period) && period > 0.0) {
                        thi = t0 + 3.0 * period;
                    } else {
                        thi = t0 + 1.0;
                    }
                } else {
                    thi = tend;
                }
                Traj steps = nhppSteps(nh, t0, thi);
                double[] rowmult = new double[steps.t.length];
                for (int j = 0; j < rowmult.length; j++) {
                    rowmult[j] = steps.m[j][0] / nominal;
                }
                Traj thisTraj = expandRows(numEvents, rowmult, steps.t, eventIdx, qIndices, kic, i, c);
                nhppTraj = merge(nhppTraj, thisTraj, numEvents);
            }
        }

        // -- channel (3): explicit per-(station,class) rate trajectories -----------
        Traj schedTraj = null;
        if (cfg.rate_sched != null && !cfg.rate_sched.isEmpty()) {
            for (int s = 0; s < cfg.rate_sched.size(); s++) {
                RateEntry entry = cfg.rate_sched.get(s);
                int i = entry.station;
                int c = entry.jobclass;
                if (i < 0 || c < 0 || i >= enabled.length || c >= enabled[i].length || !enabled[i][c]) {
                    continue;
                }
                double nominal;
                if (entry.nominal != null) {
                    nominal = entry.nominal.doubleValue();
                } else {
                    nominal = nominalRate(mu, stations, jobclasses, i, c);
                }
                if (!(nominal > 0.0)) {
                    continue;
                }
                if (entry.tgrid == null || entry.rates == null || entry.tgrid.length != entry.rates.length
                        || entry.tgrid.length == 0) {
                    line_error(FluidRateMultiplier.class.getName(),
                            "rate_sched entry must carry tgrid and rates of the same non-zero length.");
                }
                double[] rowmult = new double[entry.rates.length];
                for (int j = 0; j < rowmult.length; j++) {
                    rowmult[j] = entry.rates[j] / nominal;
                }
                Traj thisTraj = expandRows(numEvents, rowmult, entry.tgrid.clone(), eventIdx, qIndices, kic, i, c);
                schedTraj = merge(schedTraj, thisTraj, numEvents);
            }
        }

        // -- compose all channels -------------------------------------------------
        Traj composed = merge(merge(user, nhppTraj, numEvents), schedTraj, numEvents);
        if (composed == null) {
            return null;
        }
        return new FluidRateMultiplier(composed.t, composed.m, numEvents);
    }

    /** The per-event multiplier in force at time {@code t}. */
    public double[] evalAt(double t) {
        return FluidInterp.interp(tgrid, mmat, t, numEvents);
    }

    /** The time grid of the multiplier trajectory. */
    public double[] getTimeGrid() {
        return tgrid.clone();
    }

    private static double nominalRate(Map<Station, Map<JobClass, Matrix>> mu,
                                      List<Station> stations,
                                      List<JobClass> jobclasses,
                                      int i, int c) {
        Map<JobClass, Matrix> muI = mu.get(stations.get(i));
        if (muI == null) {
            return Double.NaN;
        }
        Matrix muIC = muI.get(jobclasses.get(c));
        if (muIC == null || muIC.isEmpty()) {
            return Double.NaN;
        }
        return muIC.get(0, 0);
    }

    /**
     * Expands a per-(station,class) multiplier trajectory to the event space:
     * every event sourced at (i,c) across its service phases carries the
     * multiplier, every other event carries 1.
     */
    private static Traj expandRows(int numEvents,
                                   double[] rowmult,
                                   double[] segT,
                                   Matrix eventIdx,
                                   Matrix qIndices,
                                   Matrix kic,
                                   int i, int c) {
        boolean[] rows = new boolean[numEvents];
        int base = (int) qIndices.get(i, c);
        int nphases = (int) kic.get(i, c);
        for (int e = 0; e < numEvents; e++) {
            int idx = (int) eventIdx.get(e, 0);
            if (idx >= base && idx < base + nphases) {
                rows[e] = true;
            }
        }
        double[][] m = new double[segT.length][numEvents];
        for (int j = 0; j < segT.length; j++) {
            for (int e = 0; e < numEvents; e++) {
                m[j][e] = rows[e] ? rowmult[j] : 1.0;
            }
        }
        return new Traj(segT, m);
    }

    /**
     * Step-faithful (time, rate) sampling of a piecewise-constant NHPP
     * intensity over {@code [t0, thi]}. Each segment contributes two samples,
     * at its start and just before its end, so clamped-linear interpolation
     * reproduces the step with a negligible transition ramp. The returned
     * trajectory carries the rate in column 0 of each row.
     */
    private static Traj nhppSteps(NHPP nh, double t0, double thi) {
        double[] bp = nh.getBreakpoints();
        double period = nh.getPeriod();
        TreeSet<Double> bounds = new TreeSet<Double>();
        if (nh.isCyclic() && Double.isFinite(period) && period > 0.0) {
            int kmax = (int) Math.ceil((thi - t0) / period) + 2;
            for (int k = -1; k <= kmax; k++) {
                for (int b = 0; b < bp.length; b++) {
                    double v = bp[b] + k * period;
                    if (v > t0 && v < thi) {
                        bounds.add(v);
                    }
                }
            }
        } else {
            for (int b = 0; b < bp.length; b++) {
                double v = bp[b];
                if (v > t0 && v < thi) {
                    bounds.add(v);
                }
            }
        }
        bounds.add(t0);
        bounds.add(thi);
        List<Double> blist = new ArrayList<Double>(bounds);

        double neps = Math.max(1e-9, 1e-6 * (thi - t0));
        int nb = blist.size() - 1;
        double[] segT = new double[2 * nb];
        double[][] segR = new double[2 * nb][1];
        for (int k = 0; k < nb; k++) {
            double a = blist.get(k).doubleValue();
            double b = blist.get(k + 1).doubleValue();
            double r = nh.getRateAt((a + b) / 2.0);
            segT[2 * k] = a;
            segR[2 * k][0] = r;
            segT[2 * k + 1] = Math.max(a + neps, b - neps);
            segR[2 * k + 1][0] = r;
        }
        return new Traj(segT, segR);
    }

    /**
     * Merges two event-multiplier trajectories onto the union time grid by
     * elementwise product (identity where a channel is silent).
     */
    private static Traj merge(Traj a, Traj b, int numEvents) {
        if (a == null) {
            return b;
        }
        if (b == null) {
            return a;
        }
        TreeSet<Double> union = new TreeSet<Double>();
        for (int j = 0; j < a.t.length; j++) {
            union.add(a.t[j]);
        }
        for (int j = 0; j < b.t.length; j++) {
            union.add(b.t[j]);
        }
        double[] tg = new double[union.size()];
        int p = 0;
        for (Double v : union) {
            tg[p++] = v.doubleValue();
        }
        double[][] mg = new double[tg.length][numEvents];
        for (int j = 0; j < tg.length; j++) {
            double[] ma = FluidInterp.interp(a.t, a.m, tg[j], numEvents);
            double[] mb = FluidInterp.interp(b.t, b.m, tg[j], numEvents);
            for (int e = 0; e < numEvents; e++) {
                mg[j][e] = ma[e] * mb[e];
            }
        }
        return new Traj(tg, mg);
    }

    /** A sampled trajectory: ascending time grid and one value row per time. */
    private static class Traj {
        final double[] t;
        final double[][] m;

        Traj(double[] t, double[][] m) {
            this.t = t;
            this.m = m;
        }
    }

    /**
     * An NHPP source intensity attached to a (station, class) of the sn index
     * space. The nominal rate baked into the closing rate base is taken from
     * {@code mu(station, jobclass)}.
     */
    public static class NhppEntry {
        public int station;
        public int jobclass;
        public NHPP nhpp;

        public NhppEntry(int station, int jobclass, NHPP nhpp) {
            this.station = station;
            this.jobclass = jobclass;
            this.nhpp = nhpp;
        }
    }

    /**
     * An explicit rate trajectory attached to a (station, class) of the sn
     * index space. When {@code nominal} is null the nominal rate baked into
     * the closing rate base is taken from {@code mu(station, jobclass)}.
     */
    public static class RateEntry {
        public int station;
        public int jobclass;
        public double[] tgrid;
        public double[] rates;
        public Double nominal;

        public RateEntry(int station, int jobclass, double[] tgrid, double[] rates, Double nominal) {
            this.station = station;
            this.jobclass = jobclass;
            this.tgrid = tgrid;
            this.rates = rates;
            this.nominal = nominal;
        }
    }
}
