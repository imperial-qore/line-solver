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
 * Time-inhomogeneous MARKED Markovian arrival process (MMAP_t).
 *
 * <p>The two axes of {@link MAPt} and {@link MarkedMAP} crossed: arrivals are labelled with
 * one of K marks, AND the matrices that generate them are functions of the wall clock. A
 * segment k covers [breakpoints[k], breakpoints[k+1]) and carries D0[k] together with the K
 * blocks D1^(1)[k], ..., D1^(K)[k]; D0 holds transition rates without an arrival, D1^(c) the
 * rates that generate an arrival of mark c, and
 *
 * <pre>  D0[k] + sum_c D1^(c)[k]</pre>
 *
 * is a generator in every segment. The aggregate sum_c D1^(c)[k] is the D1 of the underlying
 * {@link MAPt}, so hiding the marks recovers exactly that process.
 *
 * <p>Two horizon conventions, as for MAPt and NHPP: a cyclic schedule repeats with period
 * breakpoints[n]-breakpoints[0]; a non-cyclic one is silent outside its horizon and is
 * therefore a transient construct.
 *
 * <p><b>Reductions.</b> With K = 1 this is exactly the MAPt with the same matrices, and it is
 * held to the same constructor rules so the reduction is exact rather than merely close. With
 * one segment, or with every segment identical, it is exactly the stationary {@link MMAP}.
 *
 * <p>Like MAPt this is neither a renewal process nor a time-homogeneous one, so getSCV,
 * getSkewness, evalCDF and evalLST return NaN rather than a value that would misreport the
 * process as stationary. It deliberately does NOT extend Markovian or MarkedMAP: code gated on
 * isMarkovian reads getProcess as a single stationary pair and would silently drop the schedule.
 *
 * <p><b>Constant support.</b> The constructor requires one sparsity pattern across segments,
 * per mark block and for the off-diagonal of D0, reusing {@link MAPt#checkCommonSupport}. The
 * rule is what keeps the K = 1 reduction to MAPt exact, and it leaves the fluid path open: a
 * segment is expressed there as a per-entry multiplier on a nominal, which is undefined where
 * the nominal entry is zero.
 *
 * <p>At a Source the mark selects the class of the arriving job
 * ({@code Source.setMarkedArrival}, {@code sn.markidx}). As a SERVICE process it is sampled
 * for its duration and the mark is discarded, exactly as an MMAP is.
 *
 * <p>The MatrixCell layout is flat, because a MatrixCell cannot nest:
 *
 * <pre>  [breakpoints, K, D0_1..D0_n, D1^(1)_1..D1^(1)_n, ..., D1^(K)_1..D1^(K)_n, cyclic]</pre>
 *
 * with K carried explicitly in a 1-by-1 matrix, since the length alone gives only n(K+1) and
 * cannot separate the two. The total length is 3 + n(K+1).
 *
 * <p>References: Q.-M. He, "The versatility of MMAP[K] and the MMAP[K]/G[K]/1 queue", Queueing
 * Systems 38(4), 2001, for the marked structure; Y. M. Ko and J. Pender, "Diffusion limits for
 * the (MAP_t/Ph_t/inf)^N queueing network", Oper. Res. Lett. 45(3), 2017, for the
 * time-inhomogeneous one.
 */
public class MMAPt extends ContinuousDistribution implements Serializable {

    private static final long serialVersionUID = 1L;

    private final double[] breakpoints;
    private final List<Matrix> d0;
    /** d1k.get(c).get(j) is the block of mark c+1 in segment j. */
    private final List<List<Matrix>> d1k;
    /** Per-segment aggregate sum_c D1^(c), the D1 of the underlying MAPt. */
    private final List<Matrix> d1;
    private final boolean cyclic;
    private final MatrixCell process;

    /** Wall-clock position and phase of the next sample; see {@link #sample(int, Random)}. */
    private double sampleClock;
    private int samplePhase;
    /** The 1-based mark of the interval last returned by {@link #sample(int, Random)}. */
    private int lastMark;

    /**
     * Creates an MMAP_t with a piecewise-constant marked matrix schedule.
     *
     * @param breakpoints strictly increasing segment boundaries, length n+1
     * @param d0          per-segment no-arrival rate matrices, length n
     * @param d1k         the K mark blocks, each a list of n per-segment matrices
     * @param cyclic      whether the schedule repeats with the horizon as period
     */
    public MMAPt(double[] breakpoints, List<Matrix> d0, List<List<Matrix>> d1k, boolean cyclic) {
        super("MMAPt", 0, new Pair<Double, Double>(0.0, GlobalConstants.Inf));
        if (breakpoints == null || d0 == null || d1k == null || d0.isEmpty() || d1k.isEmpty()
                || breakpoints.length != d0.size() + 1) {
            throw new IllegalArgumentException(
                    "MMAPt: breakpoints must have one more entry than the number of segments, "
                            + "and D0 and at least one mark block must be non-empty");
        }
        int n = d0.size();
        int K = d1k.size();
        for (int c = 0; c < K; c++) {
            if (d1k.get(c) == null || d1k.get(c).size() != n) {
                throw new IllegalArgumentException(
                        "MMAPt: mark block " + (c + 1) + " has "
                                + (d1k.get(c) == null ? 0 : d1k.get(c).size())
                                + " segments against " + n + " in D0; every mark must be defined "
                                + "on the whole schedule");
            }
        }
        for (int i = 0; i + 1 < breakpoints.length; i++) {
            if (!(breakpoints[i + 1] > breakpoints[i])) {
                throw new IllegalArgumentException(
                        "MMAPt: breakpoints must be strictly increasing (at index " + i + ")");
            }
        }
        int h = d0.get(0).getNumRows();
        List<Matrix> agg = new ArrayList<Matrix>();
        for (int j = 0; j < n; j++) {
            Matrix a = d0.get(j);
            if (a.getNumRows() != h || a.getNumCols() != h) {
                throw new IllegalArgumentException(
                        "MMAPt: every D0 must be square of order " + h + "; segment " + (j + 1)
                                + " differs");
            }
            Matrix sum = new Matrix(h, h, h * h);
            for (int c = 0; c < K; c++) {
                Matrix b = d1k.get(c).get(j);
                if (b.getNumRows() != h || b.getNumCols() != h) {
                    throw new IllegalArgumentException(
                            "MMAPt: every mark block must be square of order " + h + "; mark "
                                    + (c + 1) + " of segment " + (j + 1) + " differs");
                }
                for (int p = 0; p < h; p++) {
                    for (int q = 0; q < h; q++) {
                        double v = b.get(p, q);
                        if (v < 0.0) {
                            throw new IllegalArgumentException(
                                    "MMAPt: mark block " + (c + 1) + " must be non-negative in "
                                            + "segment " + (j + 1));
                        }
                        if (v != 0.0) {
                            sum.set(p, q, sum.get(p, q) + v);
                        }
                    }
                }
            }
            for (int p = 0; p < h; p++) {
                double rowsum = 0.0;
                for (int q = 0; q < h; q++) {
                    if (p != q && a.get(p, q) < 0.0) {
                        throw new IllegalArgumentException(
                                "MMAPt: off-diagonal D0 entries must be non-negative in segment "
                                        + (j + 1));
                    }
                    rowsum += a.get(p, q) + sum.get(p, q);
                }
                if (Math.abs(rowsum) > 1e-10) {
                    throw new IllegalArgumentException(
                            "MMAPt: D0 plus the mark blocks must have zero row sums (generator) "
                                    + "in segment " + (j + 1) + "; row " + p + " sums to "
                                    + rowsum);
                }
            }
            agg.add(sum);
        }
        MAPt.checkCommonSupport(d0, true, "off-diagonal D0");
        for (int c = 0; c < K; c++) {
            MAPt.checkCommonSupport(d1k.get(c), false, "D1 of mark " + (c + 1));
        }
        boolean anyArrival = false;
        for (int j = 0; j < n && !anyArrival; j++) {
            for (int p = 0; p < h && !anyArrival; p++) {
                for (int q = 0; q < h; q++) {
                    if (agg.get(j).get(p, q) > 0.0) {
                        anyArrival = true;
                        break;
                    }
                }
            }
        }
        if (!anyArrival) {
            throw new IllegalArgumentException(
                    "MMAPt: every segment has zero arrival intensity, so no event can ever occur");
        }

        this.breakpoints = breakpoints.clone();
        this.d0 = new ArrayList<Matrix>(d0);
        this.d1k = new ArrayList<List<Matrix>>();
        for (int c = 0; c < K; c++) {
            this.d1k.add(new ArrayList<Matrix>(d1k.get(c)));
        }
        this.d1 = agg;
        this.cyclic = cyclic;

        Matrix bpMatrix = new Matrix(1, n + 1, n + 1);
        for (int i = 0; i <= n; i++) {
            bpMatrix.set(0, i, breakpoints[i]);
        }
        Matrix kMatrix = new Matrix(1, 1, 1);
        kMatrix.set(0, 0, K);
        Matrix cyclicMatrix = new Matrix(1, 1, 1);
        cyclicMatrix.set(0, 0, cyclic ? 1.0 : 0.0);
        this.process = new MatrixCell();
        this.process.set(0, bpMatrix);
        this.process.set(1, kMatrix);
        for (int j = 0; j < n; j++) {
            this.process.set(2 + j, this.d0.get(j));
        }
        for (int c = 0; c < K; c++) {
            for (int j = 0; j < n; j++) {
                this.process.set(2 + (c + 1) * n + j, this.d1k.get(c).get(j));
            }
        }
        this.process.set(2 + (K + 1) * n, cyclicMatrix);

        this.sampleClock = breakpoints[0];
        this.samplePhase = 0;
        this.lastMark = 0;
        this.mean = 1.0 / getTimeAverageRate();
        this.immediate = false;
    }

    /** Creates a cyclic MMAP_t. */
    public MMAPt(double[] breakpoints, List<Matrix> d0, List<List<Matrix>> d1k) {
        this(breakpoints, d0, d1k, true);
    }

    public double[] getBreakpoints() {
        return breakpoints.clone();
    }

    public List<Matrix> getD0Segments() {
        return new ArrayList<Matrix>(d0);
    }

    /** Per-segment aggregate sum_c D1^(c), i.e. the D1 of the underlying MAPt. */
    public List<Matrix> getD1Segments() {
        return new ArrayList<Matrix>(d1);
    }

    /**
     * Per-segment blocks of one mark.
     *
     * @param k the 1-based mark index
     * @return the n per-segment matrices of that mark
     */
    public List<Matrix> getD1Segments(int k) {
        if (k < 1 || k > d1k.size()) {
            throw new IllegalArgumentException("MMAPt: mark index out of range: " + k);
        }
        return new ArrayList<Matrix>(d1k.get(k - 1));
    }

    /** All mark blocks, mark-major: get(c).get(j) is mark c+1 in segment j. */
    public List<List<Matrix>> getMarkSegments() {
        List<List<Matrix>> out = new ArrayList<List<Matrix>>();
        for (int c = 0; c < d1k.size(); c++) {
            out.add(new ArrayList<Matrix>(d1k.get(c)));
        }
        return out;
    }

    public int getNumberOfTypes() {
        return d1k.size();
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

    /**
     * The UNMARKED schedule, i.e. the MAPt whose D1 is the per-segment aggregate.
     *
     * <p>Hiding the marks is the exact operation: an arrival of the MMAPt is an arrival of
     * this process regardless of its label.
     *
     * @return the aggregate MAPt
     */
    public MAPt toMAPt() {
        return new MAPt(breakpoints.clone(), new ArrayList<Matrix>(d0),
                new ArrayList<Matrix>(d1), cyclic);
    }

    /**
     * The MARGINAL schedule of one mark: arrivals fire only on that mark's blocks, while the
     * other marks' transitions become hidden phase changes. Segment by segment this is
     * MAPt(D0 + D1_agg - D1^(k), D1^(k)), the time-varying analogue of
     * {@link MarkedMAP#toMAPs(int)}.
     *
     * @param k the 1-based mark index
     * @return the marginal MAPt of that mark
     */
    public MAPt toMAPts(int k) {
        if (k < 1 || k > d1k.size()) {
            throw new IllegalArgumentException("MMAPt: mark index out of range: " + k);
        }
        int n = d0.size();
        List<Matrix> hidden = new ArrayList<Matrix>();
        for (int j = 0; j < n; j++) {
            Matrix m = d0.get(j).copy();
            m.addEq(1.0, d1.get(j));
            m.addEq(-1.0, d1k.get(k - 1).get(j));
            hidden.add(m);
        }
        return new MAPt(breakpoints.clone(), hidden,
                new ArrayList<Matrix>(d1k.get(k - 1)), cyclic);
    }

    /**
     * Width-weighted average pair over the horizon, {D0bar, D1bar}.
     *
     * <p>This is the stationary carrier of the phase structure where a solver needs a
     * time-homogeneous one; a convex combination of generators is a generator, so it is
     * itself a valid MAP.
     *
     * @return the nominal (D0, D1) pair
     */
    public MatrixCell getTimeAverageProcess() {
        return toMAPt().getTimeAverageProcess();
    }

    /**
     * Width-weighted average of one mark's blocks, which is the D1 of that mark in the
     * nominal marked process.
     *
     * @param k the 1-based mark index
     * @return the nominal block of that mark
     */
    public Matrix getTimeAverageMark(int k) {
        if (k < 1 || k > d1k.size()) {
            throw new IllegalArgumentException("MMAPt: mark index out of range: " + k);
        }
        int h = getNumberOfPhases();
        Matrix out = new Matrix(h, h, h * h);
        double total = getPeriod();
        for (int j = 0; j < d0.size(); j++) {
            double w = (breakpoints[j + 1] - breakpoints[j]) / total;
            out.addEq(w, d1k.get(k - 1).get(j));
        }
        return out;
    }

    /** Arrival rate of the time-averaged aggregate MAP. */
    public double getTimeAverageRate() {
        MatrixCell nominal = getTimeAverageProcess();
        return Map_lambda.map_lambda(nominal.get(0), nominal.get(1));
    }

    /**
     * Per-mark arrival rates of the time-averaged process.
     *
     * <p>These sum to {@link #getTimeAverageRate()}, which is the identity a marked stream
     * has to satisfy: labelling the arrivals cannot change how many there are.
     *
     * @return the K nominal per-mark rates
     */
    public double[] getTimeAverageMarkRates() {
        // Through the canonical mmap_count_lambda rather than a second implementation of
        // theta*D1k*e: the rate of a marked counting process is weighted by the STATIONARY
        // phase distribution map_piq, not by the embedded departure law, and the two differ
        // whenever the phases do not all arrive at the same rate.
        MatrixCell nominal = getTimeAverageProcess();
        MatrixCell cell = new MatrixCell();
        cell.set(0, nominal.get(0));
        cell.set(1, nominal.get(1));
        for (int c = 0; c < d1k.size(); c++) {
            cell.set(2 + c, getTimeAverageMark(c + 1));
        }
        Matrix lk = jline.api.mam.Mmap_count_lambda.mmap_count_lambda(cell);
        double[] out = new double[d1k.size()];
        for (int c = 0; c < d1k.size(); c++) {
            out[c] = lk.get(0, c);
        }
        return out;
    }

    @Override
    public double getMean() {
        return 1.0 / getTimeAverageRate();
    }

    @Override
    public double getRate() {
        return getTimeAverageRate();
    }

    /**
     * NaN: an MMAP_t is neither renewal nor time-homogeneous, so there is no i.i.d. interval
     * distribution for an SCV to summarise. Returning the SCV of the time-averaged process
     * would report a time-varying one as stationary to every consumer of sn.scv.
     */
    @Override
    public double getSCV() {
        return Double.NaN;
    }

    /** NaN; see {@link #getSCV()}. */
    @Override
    public double getSkewness() {
        return Double.NaN;
    }

    /** NaN; see {@link #getSCV()}. */
    @Override
    public double evalCDF(double t) {
        return Double.NaN;
    }

    /** NaN; see {@link #getSCV()}. */
    @Override
    public double evalLST(double s) {
        return Double.NaN;
    }

    @Override
    public MatrixCell getProcess() {
        return process;
    }

    /** Restarts the sample path at the schedule start. */
    public void resetSampleClock() {
        this.sampleClock = breakpoints[0];
        this.samplePhase = 0;
        this.lastMark = 0;
    }

    /** The 1-based mark of the interval last returned by {@link #sample(int, Random)}. */
    public int getLastMark() {
        return lastMark;
    }

    /**
     * Draws n successive interarrival times along ONE sample path.
     *
     * <p>Both the intensity and the phase depend on absolute time, so this advances an
     * internal clock and phase across calls; use {@link #resetSampleClock()} to restart. A
     * non-cyclic schedule that runs out returns 0 for every remaining sample.
     *
     * @param n      the number of samples
     * @param random the generator
     * @return the n interarrival times
     */
    @Override
    public double[] sample(int n, Random random) {
        double[] out = new double[n];
        for (int i = 0; i < n; i++) {
            double[] step = nextArrival(sampleClock, samplePhase, random);
            out[i] = step[0];
            if (step[0] <= 0.0) {
                break;
            }
            sampleClock += step[0];
            samplePhase = (int) step[1];
            lastMark = (int) step[2];
        }
        return out;
    }

    @Override
    public double[] sample(int n) {
        return sample(n, new Random());
    }

    /**
     * Time to the next arrival from wall clock {@code from} in {@code phase}, the phase after
     * it and the 1-based mark it carries.
     *
     * <p>Exact: within a segment the phase process is a homogeneous CTMC, and by the
     * memoryless property the residual holding time may be redrawn at a breakpoint, so the
     * boundary is crossed by advancing the clock and resampling under the new matrices. The
     * mark is decided by WHICH block's transition fired, in the same competing-transitions
     * draw that ends the interval, so it is not an extra layer on top of an unmarked walk.
     *
     * @param from   the wall-clock start
     * @param phase  the current phase
     * @param random the generator
     * @return {interval, next phase, mark}; the interval is 0 once a non-cyclic horizon is
     *         exhausted
     */
    public double[] nextArrival(double from, int phase, Random random) {
        double elapsed = 0.0;
        double pos = from;
        int h = getNumberOfPhases();
        int K = d1k.size();
        int guard = 0;
        while (guard++ < 1000000) {
            int idx = getSegmentIndexAt(pos);
            if (idx < 0) {
                return new double[] {0.0, phase, 0};
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
            Matrix dz = d0.get(idx);
            double total = -dz.get(phase, phase);
            if (total <= 0.0) {
                // An absorbing phase in this segment: only a boundary can free it.
                if (!cyclic && idx == d0.size() - 1) {
                    return new double[] {0.0, phase, 0};
                }
                elapsed += toBoundary;
                pos += toBoundary;
                continue;
            }
            double holding = -Math.log(1.0 - random.nextDouble()) / total;
            if (holding >= toBoundary) {
                if (!cyclic && idx == d0.size() - 1) {
                    return new double[] {0.0, phase, 0};
                }
                elapsed += toBoundary;
                pos += toBoundary;
                continue;
            }
            elapsed += holding;
            pos += holding;
            // Competing transitions out of the current phase: every mark's arrivals first, in
            // mark order, then the hidden phase changes. The winning column says both where
            // the chain goes and, when it is an arrival, which mark the event carries.
            double u = random.nextDouble() * total;
            double cum = 0.0;
            // DESTINATION-MAJOR, then MARK-MINOR. The order is load-bearing, not cosmetic:
            // only this way does the running total after all marks of a destination equal
            // the aggregate total after that destination, which is what makes a K = 1 MMAPt
            // pick the same destination the MAPt would for the same uniform. It also matches
            // the MATLAB, python and C++ samplers and this codebase's own engine walk,
            // Solver_ssj.sampleMMAPtInterval; before 2026-09-12 this method alone scanned
            // mark-major and so disagreed with all four on one seed.
            for (int j = 0; j < h; j++) {
                for (int c = 0; c < K; c++) {
                    cum += d1k.get(c).get(idx).get(phase, j);
                    if (u < cum) {
                        return new double[] {elapsed, j, c + 1};
                    }
                }
            }
            boolean moved = false;
            for (int j = 0; j < h; j++) {
                if (j == phase) {
                    continue;
                }
                cum += dz.get(phase, j);
                if (u < cum) {
                    phase = j;
                    moved = true;
                    break;
                }
            }
            if (!moved) {
                // Rounding left u at or past the total: fall back to the last mark that has
                // any mass out of this phase, as MAPt does, rather than redrawing a holding
                // time and biasing the interval upwards.
                for (int j = h - 1; j >= 0; j--) {
                    for (int c = K - 1; c >= 0; c--) {
                        if (d1k.get(c).get(idx).get(phase, j) > 0.0) {
                            return new double[] {elapsed, j, c + 1};
                        }
                    }
                }
            }
        }
        return new double[] {elapsed, phase, 0};
    }
}
