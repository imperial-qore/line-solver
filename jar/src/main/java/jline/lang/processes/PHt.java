/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import jline.GlobalConstants;
import jline.api.mam.Map_lambda;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Time-inhomogeneous phase-type distribution (Ph_t).
 *
 * <p>Following Ko and Pender, "Diffusion limits for the (MAP_t/Ph_t/inf)^N queueing network",
 * Oper. Res. Lett. 45 (2017) 248-253, a Ph_t is an ordinary phase-type distribution whose
 * initial vector and sub-generator are functions of the wall clock, alpha(t) and S(t), required
 * only to be locally integrable. This class realises that definition with a piecewise-constant
 * schedule: segment k covers [breakpoints[k], breakpoints[k+1]) and carries the pair
 * (alpha[k], S[k]). The exit vector is s(t) = -S(t)e.
 *
 * <p>Because both the phase and the elapsed service depend on absolute time, a Ph_t service
 * time is a function of the epoch at which service starts: {@link #sampleFrom(double, Random)}
 * is the operative sampler, and {@link #sample(int, Random)} walks one path.
 *
 * <p>Like MAPt this does NOT extend Markovian, and the scalar summaries getSCV, getSkewness,
 * evalCDF and evalLST return NaN, the distribution of a service time being different at every
 * start epoch.
 *
 * <p>The MatrixCell layout is flat, [breakpoints, alpha_1..alpha_n, S_1..S_n, cyclic], each
 * alpha stored as a 1-by-h row so that every element of the cell is a matrix.
 */
public class PHt extends ContinuousDistribution implements Serializable {

    private static final long serialVersionUID = 1L;

    private final double[] breakpoints;
    private final List<Matrix> alpha;
    private final List<Matrix> subgen;
    private final boolean cyclic;
    private final MatrixCell process;

    /** Wall-clock position of the next sample; see {@link #sample(int, Random)}. */
    private double sampleClock;

    /**
     * Creates a Ph_t with a piecewise-constant schedule.
     *
     * @param breakpoints strictly increasing segment boundaries, length n+1
     * @param alpha       per-segment initial probability rows (1-by-h), length n
     * @param subgen      per-segment sub-generators (h-by-h), length n
     * @param cyclic      whether the schedule repeats with the horizon as period
     */
    public PHt(double[] breakpoints, List<Matrix> alpha, List<Matrix> subgen, boolean cyclic) {
        super("PHt", 0, new Pair<Double, Double>(0.0, GlobalConstants.Inf));
        if (breakpoints == null || alpha == null || subgen == null || subgen.isEmpty()
                || alpha.size() != subgen.size() || breakpoints.length != subgen.size() + 1) {
            throw new IllegalArgumentException(
                    "PHt: breakpoints must have one more entry than the number of segments, "
                            + "and alpha and S must be non-empty lists of equal length");
        }
        for (int i = 0; i + 1 < breakpoints.length; i++) {
            if (!(breakpoints[i + 1] > breakpoints[i])) {
                throw new IllegalArgumentException(
                        "PHt: breakpoints must be strictly increasing (at index " + i + ")");
            }
        }
        int h = subgen.get(0).getNumRows();
        for (int k = 0; k < subgen.size(); k++) {
            Matrix a = alpha.get(k);
            Matrix s = subgen.get(k);
            if (s.getNumRows() != h || s.getNumCols() != h || a.length() != h) {
                throw new IllegalArgumentException(
                        "PHt: every S must be square of order " + h
                                + " with a matching alpha; segment " + (k + 1) + " differs");
            }
            double total = 0.0;
            for (int i = 0; i < h; i++) {
                double ai = a.get(0, i);
                if (ai < 0.0) {
                    throw new IllegalArgumentException(
                            "PHt: alpha must be a probability vector in segment " + (k + 1));
                }
                total += ai;
            }
            if (Math.abs(total - 1.0) > 1e-10) {
                throw new IllegalArgumentException(
                        "PHt: alpha must be a probability vector in segment " + (k + 1));
            }
            for (int i = 0; i < h; i++) {
                double rowsum = 0.0;
                for (int j = 0; j < h; j++) {
                    if (i != j && s.get(i, j) < 0.0) {
                        throw new IllegalArgumentException(
                                "PHt: off-diagonal S entries must be non-negative in segment "
                                        + (k + 1));
                    }
                    rowsum += s.get(i, j);
                }
                if (-rowsum < -1e-10) {
                    throw new IllegalArgumentException(
                            "PHt: S must have non-positive row sums in segment " + (k + 1));
                }
            }
        }
        MAPt.checkCommonSupport(subgen, true, "off-diagonal S");
        MAPt.checkCommonSupport(alpha, false, "alpha");
        List<Matrix> exits = new ArrayList<Matrix>();
        for (int k = 0; k < subgen.size(); k++) {
            exits.add(exitVector(subgen.get(k)));
        }
        MAPt.checkCommonSupport(exits, false, "exit vector");
        boolean anyExit = false;
        for (int k = 0; k < exits.size() && !anyExit; k++) {
            for (int i = 0; i < h; i++) {
                if (exits.get(k).get(i, 0) > 0.0) {
                    anyExit = true;
                    break;
                }
            }
        }
        if (!anyExit) {
            throw new IllegalArgumentException(
                    "PHt: every segment has zero exit rate, so service can never complete");
        }

        this.breakpoints = breakpoints.clone();
        this.alpha = new ArrayList<Matrix>(alpha);
        this.subgen = new ArrayList<Matrix>(subgen);
        this.cyclic = cyclic;

        int n = subgen.size();
        Matrix bpMatrix = new Matrix(1, n + 1, n + 1);
        for (int i = 0; i <= n; i++) {
            bpMatrix.set(0, i, breakpoints[i]);
        }
        Matrix cyclicMatrix = new Matrix(1, 1, 1);
        cyclicMatrix.set(0, 0, cyclic ? 1.0 : 0.0);
        this.process = new MatrixCell();
        this.process.set(0, bpMatrix);
        for (int k = 0; k < n; k++) {
            this.process.set(1 + k, this.alpha.get(k));
        }
        for (int k = 0; k < n; k++) {
            this.process.set(1 + n + k, this.subgen.get(k));
        }
        this.process.set(1 + 2 * n, cyclicMatrix);

        this.sampleClock = breakpoints[0];
        this.mean = 1.0 / getTimeAverageRate();
        this.immediate = false;
    }

    /** Creates a cyclic Ph_t. */
    public PHt(double[] breakpoints, List<Matrix> alpha, List<Matrix> subgen) {
        this(breakpoints, alpha, subgen, true);
    }

    private static Matrix exitVector(Matrix s) {
        int h = s.getNumRows();
        Matrix out = new Matrix(h, 1, h);
        for (int i = 0; i < h; i++) {
            double rowsum = 0.0;
            for (int j = 0; j < h; j++) {
                rowsum += s.get(i, j);
            }
            out.set(i, 0, -rowsum);
        }
        return out;
    }

    public double[] getBreakpoints() {
        return breakpoints.clone();
    }

    public List<Matrix> getAlphaSegments() {
        return new ArrayList<Matrix>(alpha);
    }

    public List<Matrix> getSSegments() {
        return new ArrayList<Matrix>(subgen);
    }

    public boolean isCyclic() {
        return cyclic;
    }

    public int getNumSegments() {
        return subgen.size();
    }

    public int getNumberOfPhases() {
        return subgen.get(0).getNumRows();
    }

    /** Horizon length, which is the period when cyclic. */
    public double getPeriod() {
        return breakpoints[breakpoints.length - 1] - breakpoints[0];
    }

    /** Index of the segment in force at t, or -1 past a non-cyclic horizon. */
    public int getSegmentIndexAt(double t) {
        double period = getPeriod();
        double offset = t - breakpoints[0];
        if (cyclic) {
            offset = offset % period;
            if (offset < 0) {
                offset += period;
            }
        } else if (offset < 0.0 || offset >= period) {
            return -1;
        }
        double pos = breakpoints[0] + offset;
        for (int k = 0; k < subgen.size(); k++) {
            if (pos < breakpoints[k + 1]) {
                return k;
            }
        }
        return subgen.size() - 1;
    }

    /** alpha in force at t; the last segment's row past a non-cyclic horizon. */
    public Matrix getAlphaAt(double t) {
        int idx = getSegmentIndexAt(t);
        return idx < 0 ? alpha.get(alpha.size() - 1) : alpha.get(idx);
    }

    /** S in force at t; the zero matrix past a non-cyclic horizon. */
    public Matrix getSAt(double t) {
        int idx = getSegmentIndexAt(t);
        if (idx < 0) {
            return new Matrix(getNumberOfPhases(), getNumberOfPhases());
        }
        return subgen.get(idx);
    }

    /**
     * Width-weighted average nominal as a (D0, D1) MAP pair, D1 = s*alpha.
     *
     * <p>A convex combination of sub-generators is a sub-generator and of probability vectors a
     * probability vector, so the nominal is a valid phase-type; the fluid carrier reads it as
     * the equivalent MAP.
     */
    public MatrixCell getTimeAverageProcessMAP() {
        int h = getNumberOfPhases();
        Matrix abar = new Matrix(1, h, h);
        Matrix sbar = new Matrix(h, h);
        double total = getPeriod();
        for (int k = 0; k < subgen.size(); k++) {
            double w = (breakpoints[k + 1] - breakpoints[k]) / total;
            for (int i = 0; i < h; i++) {
                abar.set(0, i, abar.get(0, i) + w * alpha.get(k).get(0, i));
                for (int j = 0; j < h; j++) {
                    sbar.set(i, j, sbar.get(i, j) + w * subgen.get(k).get(i, j));
                }
            }
        }
        Matrix exit = exitVector(sbar);
        Matrix d1 = new Matrix(h, h);
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < h; j++) {
                d1.set(i, j, exit.get(i, 0) * abar.get(0, j));
            }
        }
        MatrixCell out = new MatrixCell();
        out.set(0, sbar);
        out.set(1, d1);
        return out;
    }

    /**
     * Completion rate of the time-averaged phase-type.
     *
     * <p>A phase-type read as a MAP with D1 = s*alpha has arrival rate 1/mean, so the shared
     * MAP routine answers this without a separate linear solve.
     */
    public double getTimeAverageRate() {
        MatrixCell nominal = getTimeAverageProcessMAP();
        return Map_lambda.map_lambda(nominal.get(0), nominal.get(1));
    }

    /** Completion rate of the phase-type in force at t; zero past a non-cyclic horizon. */
    public double getRateAt(double t) {
        int idx = getSegmentIndexAt(t);
        if (idx < 0) {
            return 0.0;
        }
        Matrix s = subgen.get(idx);
        Matrix a = alpha.get(idx);
        Matrix exit = exitVector(s);
        int h = s.getNumRows();
        Matrix d1 = new Matrix(h, h);
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < h; j++) {
                d1.set(i, j, exit.get(i, 0) * a.get(0, j));
            }
        }
        return Map_lambda.map_lambda(s, d1);
    }

    /**
     * The parameterisation of the process; the scalar summaries are not. Model compilation
     * recognises a schedule-bearing process by this capability rather than by class name.
     */
    public RateSchedule getRateSchedule() {
        return new RateSchedule(breakpoints, alpha, subgen, cyclic);
    }

    /** Mean of the time-averaged phase-type. */
    public double getMean() {
        return 1.0 / getTimeAverageRate();
    }

    public double getRate() {
        return getTimeAverageRate();
    }

    /**
     * NaN: the service-time distribution differs at every start epoch, so there is no single
     * i.i.d. law for an SCV to summarise.
     */
    public double getSCV() {
        return Double.NaN;
    }

    /** NaN; see {@link #getSCV()}. */
    public double getSkewness() {
        return Double.NaN;
    }

    /** NaN; see {@link #getSCV()}. */
    public double evalCDF(double t) {
        return Double.NaN;
    }

    /** NaN; see {@link #getSCV()}. */
    public double evalLST(double s) {
        return Double.NaN;
    }

    public MatrixCell getProcess() {
        return process;
    }

    /** Restart the sample path at the schedule start. */
    public void resetSampleClock() {
        this.sampleClock = breakpoints[0];
    }

    /**
     * Draws n successive service times along ONE sample path: sample i starts where sample i-1
     * completed, not at a fixed epoch.
     */
    public double[] sample(int n, Random random) {
        double[] out = new double[n];
        for (int i = 0; i < n; i++) {
            double interval = sampleFrom(sampleClock, random);
            out[i] = interval;
            if (interval <= 0.0) {
                break; // horizon exhausted: service can never complete
            }
            sampleClock += interval;
        }
        return out;
    }

    /**
     * Service time for a job whose service starts at wall clock t0.
     *
     * <p>Exact: within a segment the phase process is a homogeneous absorbing CTMC, and by the
     * memoryless property the residual holding time may be redrawn at a breakpoint.
     */
    public double sampleFrom(double t0, Random random) {
        int idx = getSegmentIndexAt(t0);
        if (idx < 0) {
            return 0.0;
        }
        int h = getNumberOfPhases();
        double u0 = random.nextDouble();
        double cum0 = 0.0;
        int phase = h - 1;
        for (int i = 0; i < h; i++) {
            cum0 += alpha.get(idx).get(0, i);
            if (u0 < cum0) {
                phase = i;
                break;
            }
        }
        double elapsed = 0.0;
        double pos = t0;
        while (true) {
            idx = getSegmentIndexAt(pos);
            if (idx < 0) {
                return 0.0;
            }
            double period = getPeriod();
            double offset = pos - breakpoints[0];
            if (cyclic) {
                offset = offset % period;
                if (offset < 0) {
                    offset += period;
                }
            }
            double toBoundary = (breakpoints[idx + 1] - breakpoints[0]) - offset;
            Matrix s = subgen.get(idx);
            double total = -s.get(phase, phase);
            if (total <= 0.0) {
                if (!cyclic && idx == subgen.size() - 1) {
                    return 0.0;
                }
                elapsed += toBoundary;
                pos += toBoundary;
                continue;
            }
            double holding = -Math.log(1.0 - random.nextDouble()) / total;
            if (holding >= toBoundary) {
                if (!cyclic && idx == subgen.size() - 1) {
                    return 0.0;
                }
                elapsed += toBoundary;
                pos += toBoundary;
                continue;
            }
            elapsed += holding;
            pos += holding;
            // Competing transitions: absorption first, then phase changes.
            double exit = 0.0;
            for (int j = 0; j < h; j++) {
                exit += s.get(phase, j);
            }
            exit = -exit;
            double u = random.nextDouble() * total;
            double cum = exit;
            if (u < cum) {
                return elapsed;
            }
            int chosen = -1;
            for (int j = 0; j < h; j++) {
                if (j == phase) {
                    continue;
                }
                cum += s.get(phase, j);
                if (u < cum) {
                    chosen = j;
                    break;
                }
            }
            if (chosen < 0) {
                return elapsed;
            }
            phase = chosen;
        }
    }

    @Override
    public String toString() {
        return "jline.PHt(" + getNumSegments() + " segments, " + getNumberOfPhases()
                + " phases, " + (cyclic ? "cyclic" : "non-cyclic")
                + ", avgRate=" + getTimeAverageRate() + ")";
    }

    /**
     * Capability type exposed so that model compilation can recognise a schedule-bearing
     * process by capability rather than by class name.
     */
    public static final class RateSchedule implements Serializable {
        private static final long serialVersionUID = 1L;
        private final double[] breakpoints;
        private final List<Matrix> alpha;
        private final List<Matrix> subgen;
        private final boolean cyclic;

        RateSchedule(double[] breakpoints, List<Matrix> alpha, List<Matrix> subgen,
                     boolean cyclic) {
            this.breakpoints = breakpoints.clone();
            this.alpha = new ArrayList<Matrix>(alpha);
            this.subgen = new ArrayList<Matrix>(subgen);
            this.cyclic = cyclic;
        }

        public double[] getBreakpoints() {
            return breakpoints.clone();
        }

        public List<Matrix> getAlphaSegments() {
            return new ArrayList<Matrix>(alpha);
        }

        public List<Matrix> getSSegments() {
            return new ArrayList<Matrix>(subgen);
        }

        public boolean isCyclic() {
            return cyclic;
        }
    }
}
