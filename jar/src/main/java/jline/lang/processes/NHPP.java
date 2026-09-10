/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.GlobalConstants;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.io.Serializable;
import java.util.Random;

/**
 * A non-homogeneous Poisson process (NHPP) with a piecewise-constant intensity.
 *
 * <p>The intensity is a step function of the wall clock: segment {@code i}
 * covers {@code [breakpoints[i], breakpoints[i+1])} and carries rate
 * {@code rates[i]}, so {@code breakpoints} has one more entry than
 * {@code rates}. With {@code cyclic = true} the schedule repeats, giving a
 * cyclic Poisson process.
 *
 * <p>Two horizon conventions:
 * <ul>
 *   <li><b>cyclic</b>: the schedule repeats with period
 *       {@code T = breakpoints[n] - breakpoints[0]}; the active segment at time
 *       {@code t} is found from {@code (t - breakpoints[0]) mod T}.</li>
 *   <li><b>non-cyclic</b>: the intensity is zero outside
 *       {@code [breakpoints[0], breakpoints[n])}, so the process emits nothing
 *       once the schedule is exhausted. A non-cyclic NHPP is therefore a
 *       transient construct: run to steady state it converges to the empty
 *       system, and callers should use a time span within the horizon.</li>
 * </ul>
 *
 * <p><b>This is not a renewal process.</b> Successive intervals are dependent,
 * because the position within the schedule carries over from one event to the
 * next. Accordingly the scalar summaries that presuppose an i.i.d. interval
 * distribution -- {@link #getSCV()}, {@link #getSkewness()}, {@link #evalCDF},
 * {@link #evalLST} -- are undefined here and return NaN rather than a
 * representative exponential value, which would silently misreport the process
 * as Poisson. The schedule itself is the parameterisation: read it with
 * {@link #getRateSchedule()}. {@link #getMean()} is well defined and returns the
 * arrival-stationary (Palm) mean interval {@code 1/timeAverageRate}.
 *
 * <p>The process representation stores:
 * <ul>
 *   <li>{@code getProcess().get(0)} - 1x(n+1) matrix of breakpoints</li>
 *   <li>{@code getProcess().get(1)} - 1xn matrix of rates</li>
 *   <li>{@code getProcess().get(2)} - 1x1 matrix, 1 if cyclic else 0</li>
 * </ul>
 *
 * <p><b>Solver support.</b> Only the LDES simulation engine honours the exact
 * schedule; every analytical solver rejects a model using it via the standard
 * unsupported-feature check.
 */
public class NHPP extends ContinuousDistribution implements Serializable {

    private final double[] breakpoints;
    private final double[] rates;
    private final boolean cyclic;
    private final MatrixCell process;

    /**
     * Clock position of the next event to be drawn by {@link #sample}, in wall
     * time. The intensity depends on absolute time, so successive samples must
     * advance a clock; drawing each interval independently from the marginal
     * would discard the schedule entirely.
     */
    private double sampleClock;

    /**
     * Creates an NHPP with a piecewise-constant intensity.
     *
     * @param breakpoints strictly increasing segment boundaries, length n+1
     * @param rates       non-negative rate on each segment, length n
     * @param cyclic      whether the schedule repeats with the horizon as period
     */
    public NHPP(double[] breakpoints, double[] rates, boolean cyclic) {
        super("NHPP", 0, new Pair<Double, Double>(0.0, GlobalConstants.Inf));
        if (breakpoints == null || rates == null || rates.length == 0
                || breakpoints.length != rates.length + 1) {
            throw new IllegalArgumentException(
                    "NHPP: breakpoints must be non-empty with one more entry than rates");
        }
        for (int i = 0; i < rates.length; i++) {
            if (!(breakpoints[i + 1] > breakpoints[i])) {
                throw new IllegalArgumentException(
                        "NHPP: breakpoints must be strictly increasing (at index " + i + ")");
            }
            if (!(rates[i] >= 0.0) || Double.isInfinite(rates[i])) {
                throw new IllegalArgumentException(
                        "NHPP: rate at segment " + i + " must be finite and non-negative");
            }
        }
        double mass = 0.0;
        for (int i = 0; i < rates.length; i++) {
            mass += rates[i] * (breakpoints[i + 1] - breakpoints[i]);
        }
        if (!(mass > 0.0)) {
            throw new IllegalArgumentException(
                    "NHPP: the schedule has zero total intensity, so no event can ever occur");
        }
        this.breakpoints = breakpoints.clone();
        this.rates = rates.clone();
        this.cyclic = cyclic;

        int n = rates.length;
        Matrix bpMatrix = new Matrix(1, n + 1, n + 1);
        Matrix rateMatrix = new Matrix(1, n, n);
        for (int i = 0; i < n; i++) {
            bpMatrix.set(0, i, breakpoints[i]);
            rateMatrix.set(0, i, rates[i]);
        }
        bpMatrix.set(0, n, breakpoints[n]);
        Matrix cyclicMatrix = new Matrix(1, 1, 1);
        cyclicMatrix.set(0, 0, cyclic ? 1.0 : 0.0);
        this.process = new MatrixCell();
        this.process.set(0, bpMatrix);
        this.process.set(1, rateMatrix);
        this.process.set(2, cyclicMatrix);

        this.sampleClock = breakpoints[0];
        this.mean = 1.0 / getTimeAverageRate();
        this.immediate = false;
    }

    /** Creates a cyclic NHPP. */
    public NHPP(double[] breakpoints, double[] rates) {
        this(breakpoints, rates, true);
    }

    /** Returns a copy of the segment breakpoints, length n+1. */
    public double[] getBreakpoints() {
        return breakpoints.clone();
    }

    /** Returns a copy of the per-segment rates, length n. */
    public double[] getRates() {
        return rates.clone();
    }

    /** Whether the schedule repeats with period {@link #getPeriod()}. */
    public boolean isCyclic() {
        return cyclic;
    }

    /** Returns the number of segments in the schedule. */
    public int getNumSegments() {
        return rates.length;
    }

    /** Returns the horizon length, which is the period when cyclic. */
    public double getPeriod() {
        return breakpoints[breakpoints.length - 1] - breakpoints[0];
    }

    /**
     * Returns the time-average rate over the horizon,
     * {@code sum(rates[i]*width[i]) / sum(width[i])}.
     *
     * <p>For a cyclic schedule this is the long-run arrival rate. For a
     * non-cyclic one it is the average over the active horizon only, since the
     * intensity is zero afterwards.
     */
    public double getTimeAverageRate() {
        double mass = 0.0;
        for (int i = 0; i < rates.length; i++) {
            mass += rates[i] * (breakpoints[i + 1] - breakpoints[i]);
        }
        return mass / getPeriod();
    }

    /** The rate in force at wall-clock time {@code t}; zero past a non-cyclic horizon. */
    public double getRateAt(double t) {
        double period = getPeriod();
        double offset = t - breakpoints[0];
        if (cyclic) {
            offset = offset % period;
            if (offset < 0.0) {
                offset += period;
            }
        } else if (offset < 0.0 || offset >= period) {
            return 0.0;
        }
        double pos = breakpoints[0] + offset;
        for (int i = 0; i < rates.length; i++) {
            if (pos < breakpoints[i + 1]) {
                return rates[i];
            }
        }
        return rates[rates.length - 1];
    }

    /**
     * Returns the rate schedule. This is the parameterisation of the process;
     * the scalar interval summaries are not.
     */
    public RateSchedule getRateSchedule() {
        return new RateSchedule(breakpoints, rates, cyclic);
    }

    @Override
    public double getMean() {
        return 1.0 / getTimeAverageRate();
    }

    /**
     * Returns NaN: an NHPP is not a renewal process, so there is no i.i.d.
     * interval distribution for an SCV to summarise. Returning a representative
     * value here would report a time-varying process as an exponential one to
     * every consumer of {@code sn.scv}.
     */
    @Override
    public double getSCV() {
        return Double.NaN;
    }

    /** Returns NaN; see {@link #getSCV()}. */
    @Override
    public double getSkewness() {
        return Double.NaN;
    }

    /** Returns NaN; see {@link #getSCV()}. */
    @Override
    public double evalCDF(double t) {
        return Double.NaN;
    }

    /** Returns NaN; see {@link #getSCV()}. */
    @Override
    public double evalLST(double s) {
        return Double.NaN;
    }

    @Override
    public MatrixCell getProcess() {
        return process;
    }

    /** Restarts the {@link #sample} clock at the start of the schedule. */
    public void resetSampleClock() {
        this.sampleClock = breakpoints[0];
    }

    /**
     * Draws {@code n} successive interarrival times along one sample path.
     *
     * <p>The intensity depends on absolute time, so this advances an internal
     * clock across calls: consecutive samples form a realisation of the process
     * starting at {@code breakpoints[0]}, not independent draws from a marginal.
     * Use {@link #resetSampleClock()} to restart the path. A non-cyclic schedule
     * that runs out returns 0 for every remaining sample, the intensity there
     * being zero.
     */
    @Override
    public double[] sample(int n, Random random) {
        double[] samples = new double[n];
        for (int i = 0; i < n; i++) {
            double interval = nextInterval(random);
            samples[i] = interval;
            if (interval <= 0.0) {
                // Horizon exhausted: no further event can occur.
                break;
            }
            sampleClock += interval;
        }
        return samples;
    }

    /**
     * Time from {@link #sampleClock} to the next event, by inverse transform on
     * the cumulative intensity; 0 if no further event can occur.
     */
    private double nextInterval(Random random) {
        return nextInterval(sampleClock, -Math.log(1.0 - random.nextDouble()));
    }

    /**
     * Solves {@code integral_{from}^{from+x} lambda(u) du = residual} for
     * {@code x}, walking the schedule forward and consuming the budget segment
     * by segment. Returns 0 when a non-cyclic horizon is exhausted first, which
     * callers read as "no further event".
     *
     * <p>This is exact for an NHPP: conditional on no event since the last one,
     * the residual is governed by the intensity from the current instant
     * onward, so a holding time drawn under a rate that has since changed is not
     * a sample from this process.
     *
     * @param from     wall-clock origin
     * @param residual an Exp(1) budget of cumulative intensity
     */
    public double nextInterval(double from, double residual) {
        double period = getPeriod();
        double offset = from - breakpoints[0];
        if (cyclic) {
            offset = offset % period;
            if (offset < 0.0) {
                offset += period;
            }
        } else if (offset >= period) {
            return 0.0;
        } else if (offset < 0.0) {
            offset = 0.0;
        }
        double pos = breakpoints[0] + offset;
        int idx = 0;
        while (idx < rates.length - 1 && pos >= breakpoints[idx + 1]) {
            idx++;
        }
        double elapsed = 0.0;
        while (true) {
            double remainingInSegment = breakpoints[idx + 1] - pos;
            double massInSegment = rates[idx] * remainingInSegment;
            // The rate guard also keeps a zero-rate segment from dividing 0/0 on
            // the measure-zero draw residual == 0.
            if (rates[idx] > 0.0 && massInSegment >= residual) {
                return elapsed + residual / rates[idx];
            }
            residual -= massInSegment;
            elapsed += remainingInSegment;
            idx++;
            if (idx >= rates.length) {
                if (!cyclic) {
                    return 0.0;
                }
                idx = 0;
            }
            pos = breakpoints[idx];
        }
    }

    @Override
    public String toString() {
        return String.format("jline.NHPP(%d segments, %s, avgRate=%f)",
                rates.length, cyclic ? "cyclic" : "non-cyclic", getTimeAverageRate());
    }

    /**
     * A piecewise-constant rate schedule: the parameterisation shared by every
     * schedule-bearing process. Exposed so that model compilation can recognise
     * such a process by capability rather than by class name.
     */
    public static final class RateSchedule implements Serializable {
        private final double[] breakpoints;
        private final double[] rates;
        private final boolean cyclic;

        public RateSchedule(double[] breakpoints, double[] rates, boolean cyclic) {
            this.breakpoints = breakpoints.clone();
            this.rates = rates.clone();
            this.cyclic = cyclic;
        }

        public double[] getBreakpoints() {
            return breakpoints.clone();
        }

        public double[] getRates() {
            return rates.clone();
        }

        public boolean isCyclic() {
            return cyclic;
        }
    }
}
