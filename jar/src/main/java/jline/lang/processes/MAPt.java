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
 * Time-inhomogeneous Markovian arrival process (MAP_t).
 *
 * <p>Following Ko and Pender, "Diffusion limits for the (MAP_t/Ph_t/inf)^N queueing network",
 * Oper. Res. Lett. 45 (2017) 248-253, a MAP_t is an ordinary MAP whose two matrices are
 * functions of the wall clock, D0(t) and D1(t), required only to be locally integrable. This
 * class realises that definition with a piecewise-constant schedule, which is dense in L1_loc
 * and is the form that serialises: segment k covers [breakpoints[k], breakpoints[k+1]) and
 * carries the pair (D0[k], D1[k]), so breakpoints has one more entry than the matrix lists.
 * D0 holds transition rates without an arrival, D1 the rates that generate one, and D0+D1 is a
 * generator in every segment.
 *
 * <p>Two horizon conventions, as for NHPP: a cyclic schedule repeats with period
 * breakpoints[n]-breakpoints[0]; a non-cyclic one is silent outside its horizon and is
 * therefore a transient construct.
 *
 * <p>Setting h = 1 with D0 = -lambda_k and D1 = lambda_k recovers exactly the NHPP with the
 * same breakpoints and rates.
 *
 * <p>This is neither a renewal process nor a time-homogeneous one, so the scalar summaries that
 * presuppose an i.i.d. interval distribution -- getSCV, getSkewness, evalCDF, evalLST -- return
 * NaN rather than a representative value that would misreport the process as stationary. The
 * schedule is the parameterisation: read it with {@link #getRateSchedule()}.
 *
 * <p>The class deliberately does NOT extend Markovian. Code gated on isMarkovian reads
 * getProcess as a single stationary (D0, D1) pair and would silently drop the schedule.
 *
 * <p>Constant support. The fluid solver expresses a segment as a per-entry multiplier on the
 * time-averaged nominal, and that multiplier is undefined where the nominal entry is zero, so
 * the constructor requires one sparsity pattern across segments.
 *
 * <p>The MatrixCell layout is flat, [breakpoints, D0_1..D0_n, D1_1..D1_n, cyclic], because a
 * MatrixCell cannot nest; n follows from the length.
 */
public class MAPt extends ContinuousDistribution implements Serializable {

    private static final long serialVersionUID = 1L;

    private final double[] breakpoints;
    private final List<Matrix> d0;
    private final List<Matrix> d1;
    private final boolean cyclic;
    private final MatrixCell process;

    /** Wall-clock position and phase of the next sample; see {@link #sample(int, Random)}. */
    private double sampleClock;
    private int samplePhase;

    /**
     * Creates a MAP_t with a piecewise-constant matrix schedule.
     *
     * @param breakpoints strictly increasing segment boundaries, length n+1
     * @param d0          per-segment no-arrival rate matrices, length n
     * @param d1          per-segment arrival-generating rate matrices, length n
     * @param cyclic      whether the schedule repeats with the horizon as period
     */
    public MAPt(double[] breakpoints, List<Matrix> d0, List<Matrix> d1, boolean cyclic) {
        super("MAPt", 0, new Pair<Double, Double>(0.0, GlobalConstants.Inf));
        if (breakpoints == null || d0 == null || d1 == null || d0.isEmpty()
                || d0.size() != d1.size() || breakpoints.length != d0.size() + 1) {
            throw new IllegalArgumentException(
                    "MAPt: breakpoints must have one more entry than the number of segments, "
                            + "and D0 and D1 must be non-empty lists of equal length");
        }
        for (int i = 0; i + 1 < breakpoints.length; i++) {
            if (!(breakpoints[i + 1] > breakpoints[i])) {
                throw new IllegalArgumentException(
                        "MAPt: breakpoints must be strictly increasing (at index " + i + ")");
            }
        }
        int h = d0.get(0).getNumRows();
        for (int k = 0; k < d0.size(); k++) {
            Matrix a = d0.get(k);
            Matrix b = d1.get(k);
            if (a.getNumRows() != h || a.getNumCols() != h
                    || b.getNumRows() != h || b.getNumCols() != h) {
                throw new IllegalArgumentException(
                        "MAPt: every D0 and D1 must be square of order " + h
                                + "; segment " + (k + 1) + " differs");
            }
            for (int i = 0; i < h; i++) {
                double rowsum = 0.0;
                for (int j = 0; j < h; j++) {
                    if (b.get(i, j) < 0.0) {
                        throw new IllegalArgumentException(
                                "MAPt: D1 must be non-negative in segment " + (k + 1));
                    }
                    if (i != j && a.get(i, j) < 0.0) {
                        throw new IllegalArgumentException(
                                "MAPt: off-diagonal D0 entries must be non-negative in segment "
                                        + (k + 1));
                    }
                    rowsum += a.get(i, j) + b.get(i, j);
                }
                if (Math.abs(rowsum) > 1e-10) {
                    throw new IllegalArgumentException(
                            "MAPt: D0+D1 must have zero row sums (generator) in segment " + (k + 1));
                }
            }
        }
        checkCommonSupport(d0, true, "off-diagonal D0");
        checkCommonSupport(d1, false, "D1");
        boolean anyArrival = false;
        for (int k = 0; k < d1.size() && !anyArrival; k++) {
            for (int i = 0; i < h && !anyArrival; i++) {
                for (int j = 0; j < h; j++) {
                    if (d1.get(k).get(i, j) > 0.0) {
                        anyArrival = true;
                        break;
                    }
                }
            }
        }
        if (!anyArrival) {
            throw new IllegalArgumentException(
                    "MAPt: every segment has zero arrival intensity, so no event can ever occur");
        }

        this.breakpoints = breakpoints.clone();
        this.d0 = new ArrayList<Matrix>(d0);
        this.d1 = new ArrayList<Matrix>(d1);
        this.cyclic = cyclic;

        int n = d0.size();
        Matrix bpMatrix = new Matrix(1, n + 1, n + 1);
        for (int i = 0; i <= n; i++) {
            bpMatrix.set(0, i, breakpoints[i]);
        }
        Matrix cyclicMatrix = new Matrix(1, 1, 1);
        cyclicMatrix.set(0, 0, cyclic ? 1.0 : 0.0);
        this.process = new MatrixCell();
        this.process.set(0, bpMatrix);
        for (int k = 0; k < n; k++) {
            this.process.set(1 + k, this.d0.get(k));
        }
        for (int k = 0; k < n; k++) {
            this.process.set(1 + n + k, this.d1.get(k));
        }
        this.process.set(1 + 2 * n, cyclicMatrix);

        this.sampleClock = breakpoints[0];
        this.samplePhase = 0;
        this.mean = 1.0 / getTimeAverageRate();
        this.immediate = false;
    }

    /** Creates a cyclic MAP_t. */
    public MAPt(double[] breakpoints, List<Matrix> d0, List<Matrix> d1) {
        this(breakpoints, d0, d1, true);
    }

    /**
     * Rejects a schedule whose matrices do not share one sparsity pattern.
     *
     * <p>The fluid solver expresses a segment as a per-entry multiplier on a nominal matrix,
     * and that multiplier is undefined where the nominal entry is zero. Requiring one pattern
     * is what makes the nominal nonzero wherever any segment is.
     */
    public static void checkCommonSupport(List<Matrix> mats, boolean ignoreDiagonal, String label) {
        Matrix ref = mats.get(0);
        for (int k = 1; k < mats.size(); k++) {
            Matrix cur = mats.get(k);
            for (int i = 0; i < ref.getNumRows(); i++) {
                for (int j = 0; j < ref.getNumCols(); j++) {
                    if (ignoreDiagonal && i == j) {
                        continue;
                    }
                    boolean a = ref.get(i, j) != 0.0;
                    boolean b = cur.get(i, j) != 0.0;
                    if (a != b) {
                        throw new IllegalArgumentException(
                                "MAPt/PHt: the " + label + " sparsity pattern must be identical "
                                        + "across segments; segment " + (k + 1) + " differs from "
                                        + "segment 1. A schedule that switches a transition on or "
                                        + "off cannot be expressed as a per-entry multiplier on "
                                        + "the time-averaged process. Keep the entry present with "
                                        + "a small positive rate instead.");
                    }
                }
            }
        }
    }

    public double[] getBreakpoints() {
        return breakpoints.clone();
    }

    public List<Matrix> getD0Segments() {
        return new ArrayList<Matrix>(d0);
    }

    public List<Matrix> getD1Segments() {
        return new ArrayList<Matrix>(d1);
    }

    public boolean isCyclic() {
        return cyclic;
    }

    public int getNumSegments() {
        return d0.size();
    }

    public int getNumberOfPhases() {
        return d0.get(0).getNumRows();
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
        for (int k = 0; k < d0.size(); k++) {
            if (pos < breakpoints[k + 1]) {
                return k;
            }
        }
        return d0.size() - 1;
    }

    /** D0 in force at t; the zero matrix past a non-cyclic horizon. */
    public Matrix getD0At(double t) {
        int idx = getSegmentIndexAt(t);
        if (idx < 0) {
            return new Matrix(getNumberOfPhases(), getNumberOfPhases());
        }
        return d0.get(idx);
    }

    /** D1 in force at t; the zero matrix past a non-cyclic horizon. */
    public Matrix getD1At(double t) {
        int idx = getSegmentIndexAt(t);
        if (idx < 0) {
            return new Matrix(getNumberOfPhases(), getNumberOfPhases());
        }
        return d1.get(idx);
    }

    /**
     * Width-weighted average (D0bar, D1bar) over the horizon.
     *
     * <p>This is the nominal stationary MAP that carries the phase structure where a solver
     * needs a time-homogeneous carrier. A convex combination of generators is a generator, so
     * it is itself a valid MAP.
     */
    public MatrixCell getTimeAverageProcess() {
        int h = getNumberOfPhases();
        Matrix d0bar = new Matrix(h, h);
        Matrix d1bar = new Matrix(h, h);
        double total = getPeriod();
        for (int k = 0; k < d0.size(); k++) {
            double w = (breakpoints[k + 1] - breakpoints[k]) / total;
            for (int i = 0; i < h; i++) {
                for (int j = 0; j < h; j++) {
                    d0bar.set(i, j, d0bar.get(i, j) + w * d0.get(k).get(i, j));
                    d1bar.set(i, j, d1bar.get(i, j) + w * d1.get(k).get(i, j));
                }
            }
        }
        MatrixCell out = new MatrixCell();
        out.set(0, d0bar);
        out.set(1, d1bar);
        return out;
    }

    /**
     * Arrival rate of the time-averaged MAP. For h = 1 this is exactly the NHPP width-weighted
     * average intensity.
     */
    public double getTimeAverageRate() {
        MatrixCell nominal = getTimeAverageProcess();
        return Map_lambda.map_lambda(nominal.get(0), nominal.get(1));
    }

    /**
     * Arrival rate of the MAP in force at t; zero past a non-cyclic horizon. This is the
     * stationary rate of that segment, not the instantaneous conditional intensity, which
     * depends on the current phase.
     */
    public double getRateAt(double t) {
        int idx = getSegmentIndexAt(t);
        if (idx < 0) {
            return 0.0;
        }
        return Map_lambda.map_lambda(d0.get(idx), d1.get(idx));
    }

    /**
     * The parameterisation of the process; the scalar summaries are not. Model compilation
     * recognises a schedule-bearing process by this capability rather than by class name.
     */
    public RateSchedule getRateSchedule() {
        return new RateSchedule(breakpoints, d0, d1, cyclic);
    }

    /** Palm mean interval of the time-averaged MAP. */
    public double getMean() {
        return 1.0 / getTimeAverageRate();
    }

    public double getRate() {
        return getTimeAverageRate();
    }

    /**
     * NaN: a MAP_t is neither renewal nor time-homogeneous, so there is no i.i.d. interval
     * distribution for an SCV to summarise. Returning the SCV of the time-averaged MAP would
     * report a time-varying process as a stationary one to every consumer of sn.scv.
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

    /** Restart the sample path at the schedule start, in phase 1. */
    public void resetSampleClock() {
        this.sampleClock = breakpoints[0];
        this.samplePhase = 0;
    }

    /**
     * Draws n successive interarrival times along ONE sample path. Both the intensity and the
     * phase depend on absolute time, so this advances an internal clock and phase across calls.
     */
    public double[] sample(int n, Random random) {
        double[] out = new double[n];
        for (int i = 0; i < n; i++) {
            double[] step = nextArrival(sampleClock, samplePhase, random);
            out[i] = step[0];
            if (step[0] <= 0.0) {
                break; // horizon exhausted: no further arrival can occur
            }
            sampleClock += step[0];
            samplePhase = (int) step[1];
        }
        return out;
    }

    /**
     * Time to the next arrival from wall clock from in the given phase.
     *
     * <p>Exact: within a segment the phase process is a homogeneous CTMC, and by the memoryless
     * property the residual holding time may be redrawn at a breakpoint, so the boundary is
     * crossed by advancing the clock and resampling under the new matrices.
     *
     * @return {interval, phase after the arrival}; interval 0 when a non-cyclic horizon is
     *         exhausted, which callers read as "no further arrival"
     */
    public double[] nextArrival(double from, int phase, Random random) {
        int h = getNumberOfPhases();
        double elapsed = 0.0;
        double pos = from;
        while (true) {
            int idx = getSegmentIndexAt(pos);
            if (idx < 0) {
                return new double[]{0.0, phase};
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
            Matrix a = d0.get(idx);
            Matrix b = d1.get(idx);
            double total = -a.get(phase, phase);
            if (total <= 0.0) {
                if (!cyclic && idx == d0.size() - 1) {
                    return new double[]{0.0, phase};
                }
                elapsed += toBoundary;
                pos += toBoundary;
                continue;
            }
            double holding = -Math.log(1.0 - random.nextDouble()) / total;
            if (holding >= toBoundary) {
                if (!cyclic && idx == d0.size() - 1) {
                    return new double[]{0.0, phase};
                }
                elapsed += toBoundary;
                pos += toBoundary;
                continue;
            }
            elapsed += holding;
            pos += holding;
            // Competing transitions out of the current phase, arrivals first.
            double u = random.nextDouble() * total;
            double cum = 0.0;
            int chosen = -1;
            for (int j = 0; j < 2 * h; j++) {
                double w;
                if (j < h) {
                    w = b.get(phase, j);
                } else {
                    w = (j - h == phase) ? 0.0 : a.get(phase, j - h);
                }
                cum += w;
                if (u < cum) {
                    chosen = j;
                    break;
                }
            }
            if (chosen < 0) {
                // Rounding shortfall: attribute the draw to the last positive entry.
                for (int j = 2 * h - 1; j >= 0; j--) {
                    double w = (j < h) ? b.get(phase, j)
                            : ((j - h == phase) ? 0.0 : a.get(phase, j - h));
                    if (w > 0.0) {
                        chosen = j;
                        break;
                    }
                }
            }
            if (chosen < h) {
                return new double[]{elapsed, chosen};
            }
            phase = chosen - h;
        }
    }

    @Override
    public String toString() {
        return "jline.MAPt(" + getNumSegments() + " segments, " + getNumberOfPhases()
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
        private final List<Matrix> d0;
        private final List<Matrix> d1;
        private final boolean cyclic;

        RateSchedule(double[] breakpoints, List<Matrix> d0, List<Matrix> d1, boolean cyclic) {
            this.breakpoints = breakpoints.clone();
            this.d0 = new ArrayList<Matrix>(d0);
            this.d1 = new ArrayList<Matrix>(d1);
            this.cyclic = cyclic;
        }

        public double[] getBreakpoints() {
            return breakpoints.clone();
        }

        public List<Matrix> getD0Segments() {
            return new ArrayList<Matrix>(d0);
        }

        public List<Matrix> getD1Segments() {
            return new ArrayList<Matrix>(d1);
        }

        public boolean isCyclic() {
            return cyclic;
        }
    }
}
