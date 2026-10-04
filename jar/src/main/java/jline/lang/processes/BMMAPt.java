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
 * BATCH, MARKED, time-inhomogeneous Markovian arrival process (BMMAP_t).
 *
 * <p>The three axes LINE models for an arrival stream, crossed. A block is indexed by segment,
 * mark and batch size: D^(c,b)[j] holds the rates that, in segment j, release a batch of b jobs
 * ALL carrying mark c, and
 *
 * <pre>  D0[j] + sum_c sum_b D^(c,b)[j]</pre>
 *
 * is a generator in every segment. A BATCH IS HOMOGENEOUS IN ITS MARK: one epoch releases b
 * jobs that all carry mark c. This is the BMMAP[K] of the queueing literature and it is what
 * the simulation engines can release, since a batch is dispatched under a single class.
 *
 * <p>Two derived levels are stored beside the blocks and they are what keeps every existing
 * consumer working:
 *
 * <pre>  D1^(c)[j] = sum_b D^(c,b)[j]   the per-mark schedule, batches hidden
 *  D1[j]     = sum_c D1^(c)[j]     the aggregate schedule, a MAPt</pre>
 *
 * So a consumer that ignores batches reads exactly the {@link MMAPt} of D1^(c), and one that
 * ignores marks as well reads exactly the {@link MAPt} of D1.
 *
 * <p><b>Reductions.</b> With every batch size 1 this is the MMAPt with the same blocks; with
 * K = 1 it is the unmarked batch schedule; with both it is the MAPt with the same matrices;
 * with every segment identical it is the stationary {@link BMAP}. The sample-path reductions
 * are exact because {@link #nextArrival} accumulates competing transitions destination-major,
 * then mark-minor, then batch-minor, so neither label costs an extra draw and dropping the
 * batch axis leaves the scan order untouched.
 *
 * <p>Like MAPt and MMAPt this is neither renewal nor time-homogeneous, so getSCV, getSkewness,
 * evalCDF and evalLST return NaN rather than a value that would misreport the process as
 * stationary. It deliberately does NOT extend Markovian, and it deliberately does NOT extend
 * MMAPt: a consumer testing {@code instanceof MMAPt} would then accept it and silently DROP
 * the batches, which is a wrong answer rather than a refusal.
 *
 * <p><b>Rates.</b> {@link #getTimeAverageRate()} is the EVENT (batch epoch) rate and
 * {@link #getMean()} its reciprocal, the Palm mean inter-batch interval, matching BMAP whose
 * inter-batch process is the aggregate pair. The JOB rate, which is what a station's
 * throughput must balance against, is {@link #getTimeAverageJobRate()} and equals
 * sum_c sum_b b*rate(c,b).
 *
 * <p><b>Constant support.</b> The constructor requires one sparsity pattern across segments for
 * the off-diagonal of D0, for each batch-aggregated mark block D1^(c) and for the aggregate D1,
 * reusing {@link MAPt#checkCommonSupport}. Those three are precisely what {@link #toMAPt()} and
 * {@link #toMMAPt()} pass to constructors that enforce the rule themselves, so the reductions
 * stay constructible. The individual (mark, batch) blocks are DELIBERATELY EXEMPT: requiring
 * one pattern there too would forbid the composition changing with the segment, pairs in the
 * morning and singles at night, which is the one thing this family exists to express and which
 * neither reduction needs. The batch axis is DENSE in b = 1..B, as BMAP's {D0, D1, ..., Dk} is,
 * so an unused batch size is declared as a zero block rather than omitted.
 *
 * <p>The MatrixCell layout is flat, because a MatrixCell cannot nest:
 *
 * <pre>  [breakpoints, K, B, D0_1..D0_n, D^(1,1)_1..n, D^(1,2)_1..n, ..., D^(K,B)_1..n, cyclic]</pre>
 *
 * with K and B carried explicitly in 1-by-1 matrices, since the length alone gives only
 * n(1 + K*B) and cannot separate the three. The total length is 4 + n(1 + K*B). Blocks are
 * MARK-MAJOR then BATCH, matching the wire.
 *
 * <p>References: Q.-M. He, "The versatility of MMAP[K] and the MMAP[K]/G[K]/1 queue", Queueing
 * Systems 38(4), 2001, for the marked structure; D. M. Lucantoni, "New results on the single
 * server queue with a batch Markovian arrival process", Stochastic Models 7(1), 1991, for the
 * batch structure; Y. M. Ko and J. Pender, "Diffusion limits for the (MAP_t/Ph_t/inf)^N
 * queueing network", Oper. Res. Lett. 45(3), 2017, for the time-inhomogeneous one.
 */
public class BMMAPt extends ContinuousDistribution implements Serializable {

    private static final long serialVersionUID = 1L;

    private final double[] breakpoints;
    private final List<Matrix> d0;
    /** d1kb.get(c).get(b).get(j) is the block of mark c+1, batch size b+1, in segment j. */
    private final List<List<List<Matrix>>> d1kb;
    /** d1k.get(c).get(j) is the batch-aggregated block of mark c+1 in segment j. */
    private final List<List<Matrix>> d1k;
    /** Per-segment aggregate sum_c sum_b D^(c,b), the D1 of the underlying MAPt. */
    private final List<Matrix> d1;
    private final boolean cyclic;
    private final MatrixCell process;

    /** Wall-clock position and phase of the next sample; see {@link #sample(int, Random)}. */
    private double sampleClock;
    private int samplePhase;
    /** The 1-based mark and the batch size of the interval last returned by sample. */
    private int lastMark;
    private int lastBatch;

    /**
     * Creates a BMMAP_t with a piecewise-constant batch marked matrix schedule.
     *
     * @param breakpoints strictly increasing segment boundaries, length n+1
     * @param d0          per-segment no-arrival rate matrices, length n
     * @param d1kb        the blocks, mark-major then batch then segment: entry (c)(b)(j)
     *                    releases b+1 jobs of mark c+1 in segment j
     * @param cyclic      whether the schedule repeats with the horizon as period
     */
    public BMMAPt(double[] breakpoints, List<Matrix> d0, List<List<List<Matrix>>> d1kb,
                  boolean cyclic) {
        super("BMMAPt", 0, new Pair<Double, Double>(0.0, GlobalConstants.Inf));
        if (breakpoints == null || d0 == null || d1kb == null || d0.isEmpty() || d1kb.isEmpty()
                || breakpoints.length != d0.size() + 1) {
            throw new IllegalArgumentException(
                    "BMMAPt: breakpoints must have one more entry than the number of segments, "
                            + "and D0 and at least one mark block must be non-empty");
        }
        int n = d0.size();
        int K = d1kb.size();
        if (d1kb.get(0) == null || d1kb.get(0).isEmpty()) {
            throw new IllegalArgumentException("BMMAPt: mark 1 declares no batch size");
        }
        int B = d1kb.get(0).size();
        for (int c = 0; c < K; c++) {
            if (d1kb.get(c) == null || d1kb.get(c).size() != B) {
                throw new IllegalArgumentException(
                        "BMMAPt: mark " + (c + 1) + " declares "
                                + (d1kb.get(c) == null ? 0 : d1kb.get(c).size())
                                + " batch sizes against " + B + " in mark 1; the batch axis is "
                                + "dense and every mark must span it (use a zero block for an "
                                + "unused batch size)");
            }
            for (int b = 0; b < B; b++) {
                if (d1kb.get(c).get(b) == null || d1kb.get(c).get(b).size() != n) {
                    throw new IllegalArgumentException(
                            "BMMAPt: block (mark " + (c + 1) + ", batch " + (b + 1) + ") has "
                                    + (d1kb.get(c).get(b) == null ? 0 : d1kb.get(c).get(b).size())
                                    + " segments against " + n + " in D0; every block must be "
                                    + "defined on the whole schedule");
                }
            }
        }
        for (int i = 0; i + 1 < breakpoints.length; i++) {
            if (!(breakpoints[i + 1] > breakpoints[i])) {
                throw new IllegalArgumentException(
                        "BMMAPt: breakpoints must be strictly increasing (at index " + i + ")");
            }
        }
        int h = d0.get(0).getNumRows();
        List<Matrix> agg = new ArrayList<Matrix>();
        List<List<Matrix>> perMark = new ArrayList<List<Matrix>>();
        for (int c = 0; c < K; c++) {
            perMark.add(new ArrayList<Matrix>());
        }
        for (int j = 0; j < n; j++) {
            Matrix a = d0.get(j);
            if (a.getNumRows() != h || a.getNumCols() != h) {
                throw new IllegalArgumentException(
                        "BMMAPt: every D0 must be square of order " + h + "; segment " + (j + 1)
                                + " differs");
            }
            Matrix sum = new Matrix(h, h, h * h);
            for (int c = 0; c < K; c++) {
                Matrix markSum = new Matrix(h, h, h * h);
                for (int b = 0; b < B; b++) {
                    Matrix block = d1kb.get(c).get(b).get(j);
                    if (block.getNumRows() != h || block.getNumCols() != h) {
                        throw new IllegalArgumentException(
                                "BMMAPt: every block must be square of order " + h + "; block "
                                        + "(mark " + (c + 1) + ", batch " + (b + 1) + ") of "
                                        + "segment " + (j + 1) + " differs");
                    }
                    for (int p = 0; p < h; p++) {
                        for (int q = 0; q < h; q++) {
                            double v = block.get(p, q);
                            if (v < 0.0) {
                                throw new IllegalArgumentException(
                                        "BMMAPt: block (mark " + (c + 1) + ", batch " + (b + 1)
                                                + ") must be non-negative in segment " + (j + 1));
                            }
                            if (v != 0.0) {
                                markSum.set(p, q, markSum.get(p, q) + v);
                                sum.set(p, q, sum.get(p, q) + v);
                            }
                        }
                    }
                }
                perMark.get(c).add(markSum);
            }
            for (int p = 0; p < h; p++) {
                double rowsum = 0.0;
                for (int q = 0; q < h; q++) {
                    if (p != q && a.get(p, q) < 0.0) {
                        throw new IllegalArgumentException(
                                "BMMAPt: off-diagonal D0 entries must be non-negative in segment "
                                        + (j + 1));
                    }
                    rowsum += a.get(p, q) + sum.get(p, q);
                }
                if (Math.abs(rowsum) > 1e-10) {
                    throw new IllegalArgumentException(
                            "BMMAPt: D0 plus every batch block must have zero row sums "
                                    + "(generator) in segment " + (j + 1) + "; row " + p
                                    + " sums to " + rowsum);
                }
            }
            agg.add(sum);
        }
        MAPt.checkCommonSupport(d0, true, "off-diagonal D0");
        for (int c = 0; c < K; c++) {
            MAPt.checkCommonSupport(perMark.get(c), false,
                    "batch-aggregated D1 of mark " + (c + 1));
        }
        MAPt.checkCommonSupport(agg, false, "aggregate D1");
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
                    "BMMAPt: every segment has zero arrival intensity, so no event can ever "
                            + "occur");
        }

        this.breakpoints = breakpoints.clone();
        this.d0 = new ArrayList<Matrix>(d0);
        this.d1kb = new ArrayList<List<List<Matrix>>>();
        for (int c = 0; c < K; c++) {
            List<List<Matrix>> byBatch = new ArrayList<List<Matrix>>();
            for (int b = 0; b < B; b++) {
                byBatch.add(new ArrayList<Matrix>(d1kb.get(c).get(b)));
            }
            this.d1kb.add(byBatch);
        }
        this.d1k = perMark;
        this.d1 = agg;
        this.cyclic = cyclic;

        Matrix bpMatrix = new Matrix(1, n + 1, n + 1);
        for (int i = 0; i <= n; i++) {
            bpMatrix.set(0, i, breakpoints[i]);
        }
        Matrix kMatrix = new Matrix(1, 1, 1);
        kMatrix.set(0, 0, K);
        Matrix bMatrix = new Matrix(1, 1, 1);
        bMatrix.set(0, 0, B);
        Matrix cyclicMatrix = new Matrix(1, 1, 1);
        cyclicMatrix.set(0, 0, cyclic ? 1.0 : 0.0);
        this.process = new MatrixCell();
        this.process.set(0, bpMatrix);
        this.process.set(1, kMatrix);
        this.process.set(2, bMatrix);
        for (int j = 0; j < n; j++) {
            this.process.set(3 + j, this.d0.get(j));
        }
        for (int c = 0; c < K; c++) {
            for (int b = 0; b < B; b++) {
                int base = 3 + n + (c * B + b) * n;
                for (int j = 0; j < n; j++) {
                    this.process.set(base + j, this.d1kb.get(c).get(b).get(j));
                }
            }
        }
        this.process.set(3 + n + K * B * n, cyclicMatrix);

        this.sampleClock = breakpoints[0];
        this.samplePhase = 0;
        this.lastMark = 0;
        this.lastBatch = 0;
        this.mean = 1.0 / getTimeAverageRate();
        this.immediate = false;
    }

    /** Creates a cyclic BMMAP_t. */
    public BMMAPt(double[] breakpoints, List<Matrix> d0, List<List<List<Matrix>>> d1kb) {
        this(breakpoints, d0, d1kb, true);
    }

    public double[] getBreakpoints() {
        return breakpoints.clone();
    }

    public List<Matrix> getD0Segments() {
        return new ArrayList<Matrix>(d0);
    }

    /** Per-segment aggregate sum_c sum_b D^(c,b), i.e. the D1 of the underlying MAPt. */
    public List<Matrix> getD1Segments() {
        return new ArrayList<Matrix>(d1);
    }

    /**
     * Per-segment BATCH-AGGREGATED blocks of one mark, i.e. that mark's MMAPt block.
     *
     * @param k the 1-based mark index
     * @return the n per-segment matrices of that mark
     */
    public List<Matrix> getD1Segments(int k) {
        checkMark(k);
        return new ArrayList<Matrix>(d1k.get(k - 1));
    }

    /**
     * Per-segment blocks of one (mark, batch size) pair.
     *
     * @param k the 1-based mark index
     * @param b the batch size, 1-based and at most {@link #getMaxBatchSize()}
     * @return the n per-segment matrices of that block
     */
    public List<Matrix> getBatchSegments(int k, int b) {
        checkMark(k);
        checkBatch(b);
        return new ArrayList<Matrix>(d1kb.get(k - 1).get(b - 1));
    }

    /** All batch-aggregated mark blocks, mark-major: get(c).get(j) is mark c+1 in segment j. */
    public List<List<Matrix>> getMarkSegments() {
        List<List<Matrix>> out = new ArrayList<List<Matrix>>();
        for (int c = 0; c < d1k.size(); c++) {
            out.add(new ArrayList<Matrix>(d1k.get(c)));
        }
        return out;
    }

    /** All blocks, mark-major then batch then segment. */
    public List<List<List<Matrix>>> getBatchBlocks() {
        List<List<List<Matrix>>> out = new ArrayList<List<List<Matrix>>>();
        for (int c = 0; c < d1kb.size(); c++) {
            List<List<Matrix>> byBatch = new ArrayList<List<Matrix>>();
            for (int b = 0; b < d1kb.get(c).size(); b++) {
                byBatch.add(new ArrayList<Matrix>(d1kb.get(c).get(b)));
            }
            out.add(byBatch);
        }
        return out;
    }

    public int getNumberOfTypes() {
        return d1kb.size();
    }

    /** The largest batch size the schedule declares. */
    public int getMaxBatchSize() {
        return d1kb.get(0).size();
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
     * The UNMARKED, UNBATCHED schedule, i.e. the MAPt whose D1 is the per-segment aggregate.
     *
     * <p>Hiding both labels is the exact operation: an epoch of the BMMAPt is an epoch of this
     * process whatever it released.
     *
     * @return the aggregate MAPt
     */
    public MAPt toMAPt() {
        return new MAPt(breakpoints.clone(), new ArrayList<Matrix>(d0),
                new ArrayList<Matrix>(d1), cyclic);
    }

    /**
     * The BATCH-BLIND marked schedule, i.e. the MMAPt whose mark blocks are the per-mark
     * aggregates over batch size. Every epoch keeps its mark and releases one job.
     *
     * @return the batch-aggregated MMAPt
     */
    public MMAPt toMMAPt() {
        return new MMAPt(breakpoints.clone(), new ArrayList<Matrix>(d0), getMarkSegments(),
                cyclic);
    }

    /**
     * The MARGINAL schedule of one mark: epochs fire only on that mark's blocks, while the
     * other marks' transitions become hidden phase changes. Segment by segment this is
     * MAPt(D0 + D1_agg - D1^(k), D1^(k)).
     *
     * @param k the 1-based mark index
     * @return the marginal MAPt of that mark
     */
    public MAPt toMAPts(int k) {
        checkMark(k);
        int n = d0.size();
        List<Matrix> hidden = new ArrayList<Matrix>();
        for (int j = 0; j < n; j++) {
            Matrix m = d0.get(j).copy();
            m.addEq(1.0, d1.get(j));
            m.addEq(-1.0, d1k.get(k - 1).get(j));
            hidden.add(m);
        }
        return new MAPt(breakpoints.clone(), hidden, new ArrayList<Matrix>(d1k.get(k - 1)),
                cyclic);
    }

    /**
     * The width-weighted time average as a stationary BMAP, i.e. the batch structure with the
     * schedule averaged out and the marks hidden.
     *
     * @return the nominal BMAP, whose block b releases b jobs
     */
    public BMAP toBMAP() {
        MatrixCell nominal = getTimeAverageProcess();
        int h = getNumberOfPhases();
        int B = getMaxBatchSize();
        // Through the varargs constructor, NOT BMAP(MatrixCell): that overload is documented
        // as taking a raw {D0, D1, ..., Dk} but hands the cell straight to MarkedMAP, whose
        // layout puts the AGGREGATE at index 1 (getMaxBatchSize() is size()-2). The varargs
        // form builds the aggregate slot itself and is the only unambiguous one.
        Matrix[] blocks = new Matrix[B];
        for (int b = 1; b <= B; b++) {
            Matrix acc = new Matrix(h, h, h * h);
            for (int c = 1; c <= getNumberOfTypes(); c++) {
                acc.addEq(1.0, getTimeAverageBatch(c, b));
            }
            blocks[b - 1] = acc;
        }
        return new BMAP(nominal.get(0), blocks);
    }

    /**
     * Width-weighted average pair over the horizon, {D0bar, D1bar}.
     *
     * <p>This is the stationary carrier of the phase structure where a solver needs a
     * time-homogeneous one; a convex combination of generators is a generator, so it is itself
     * a valid MAP.
     *
     * @return the nominal (D0, D1) pair
     */
    public MatrixCell getTimeAverageProcess() {
        return toMAPt().getTimeAverageProcess();
    }

    /**
     * Width-weighted average of one mark's batch-aggregated blocks.
     *
     * @param k the 1-based mark index
     * @return the nominal block of that mark
     */
    public Matrix getTimeAverageMark(int k) {
        checkMark(k);
        return widthAverage(d1k.get(k - 1));
    }

    /**
     * Width-weighted average of one (mark, batch size) block.
     *
     * @param k the 1-based mark index
     * @param b the batch size
     * @return the nominal block of that pair
     */
    public Matrix getTimeAverageBatch(int k, int b) {
        checkMark(k);
        checkBatch(b);
        return widthAverage(d1kb.get(k - 1).get(b - 1));
    }

    /**
     * The EVENT rate of the time-averaged aggregate MAP, i.e. batch epochs per unit time.
     *
     * <p>This is NOT the job rate once any batch exceeds one; see
     * {@link #getTimeAverageJobRate()}.
     *
     * @return epochs per unit time
     */
    public double getTimeAverageRate() {
        MatrixCell nominal = getTimeAverageProcess();
        return Map_lambda.map_lambda(nominal.get(0), nominal.get(1));
    }

    /**
     * JOBS per unit time, sum_c sum_b b*rate(c,b).
     *
     * <p>This is the quantity a station's throughput balances against, and it exceeds
     * {@link #getTimeAverageRate()} whenever a batch larger than one has mass.
     *
     * @return jobs per unit time
     */
    public double getTimeAverageJobRate() {
        double[] byBatch = getTimeAverageBatchRates();
        double out = 0.0;
        for (int b = 0; b < byBatch.length; b++) {
            out += (b + 1) * byBatch[b];
        }
        return out;
    }

    /**
     * Per-mark EVENT rates of the time-averaged process.
     *
     * <p>These sum to {@link #getTimeAverageRate()}, which is the identity a marked stream has
     * to satisfy: labelling the epochs cannot change how many there are.
     *
     * @return the K nominal per-mark epoch rates
     */
    public double[] getTimeAverageMarkRates() {
        List<Matrix> blocks = new ArrayList<Matrix>();
        for (int c = 1; c <= getNumberOfTypes(); c++) {
            blocks.add(getTimeAverageMark(c));
        }
        return markedRates(blocks);
    }

    /**
     * Per-mark JOB rates, sum_b b*rate(c,b). These sum to
     * {@link #getTimeAverageJobRate()}.
     *
     * @return the K nominal per-mark job rates
     */
    public double[] getTimeAverageMarkJobRates() {
        int K = getNumberOfTypes();
        int B = getMaxBatchSize();
        List<Matrix> blocks = new ArrayList<Matrix>();
        for (int c = 1; c <= K; c++) {
            for (int b = 1; b <= B; b++) {
                blocks.add(getTimeAverageBatch(c, b));
            }
        }
        double[] flat = markedRates(blocks);
        double[] out = new double[K];
        for (int c = 0; c < K; c++) {
            for (int b = 0; b < B; b++) {
                out[c] += (b + 1) * flat[c * B + b];
            }
        }
        return out;
    }

    /**
     * Per-batch-size EVENT rates of the time-averaged process, marks hidden. Entry b-1 is the
     * rate of batches of size b, and they sum to {@link #getTimeAverageRate()}.
     *
     * @return the B nominal per-batch-size epoch rates
     */
    public double[] getTimeAverageBatchRates() {
        int K = getNumberOfTypes();
        int B = getMaxBatchSize();
        int h = getNumberOfPhases();
        List<Matrix> blocks = new ArrayList<Matrix>();
        for (int b = 1; b <= B; b++) {
            Matrix acc = new Matrix(h, h, h * h);
            for (int c = 1; c <= K; c++) {
                acc.addEq(1.0, getTimeAverageBatch(c, b));
            }
            blocks.add(acc);
        }
        return markedRates(blocks);
    }

    /**
     * Mean jobs per epoch of the time-averaged process.
     *
     * @return the mean batch size, 0 when the nominal has no arrival mass
     */
    public double getMeanBatchSize() {
        double[] rates = getTimeAverageBatchRates();
        double total = 0.0;
        double weighted = 0.0;
        for (int b = 0; b < rates.length; b++) {
            total += rates[b];
            weighted += (b + 1) * rates[b];
        }
        return total > 0.0 ? weighted / total : 0.0;
    }

    /**
     * The Palm mean interval BETWEEN EPOCHS of the time-averaged aggregate MAP, matching BMAP
     * whose inter-batch process is the aggregate.
     */
    @Override
    public double getMean() {
        return 1.0 / getTimeAverageRate();
    }

    @Override
    public double getRate() {
        return getTimeAverageRate();
    }

    /**
     * NaN: a BMMAP_t is neither renewal nor time-homogeneous, so there is no i.i.d. interval
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
        this.lastBatch = 0;
    }

    /** The 1-based mark of the interval last returned by {@link #sample(int, Random)}. */
    public int getLastMark() {
        return lastMark;
    }

    /** The batch size of the interval last returned by {@link #sample(int, Random)}. */
    public int getLastBatch() {
        return lastBatch;
    }

    /**
     * Draws n successive inter-epoch times along ONE sample path.
     *
     * <p>Both the intensity and the phase depend on absolute time, so this advances an internal
     * clock and phase across calls; use {@link #resetSampleClock()} to restart. A non-cyclic
     * schedule that runs out returns 0 for every remaining sample.
     *
     * @param n      the number of samples
     * @param random the generator
     * @return the n inter-epoch times
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
            lastBatch = (int) step[3];
        }
        return out;
    }

    @Override
    public double[] sample(int n) {
        return sample(n, new Random());
    }

    /**
     * Time to the next epoch from wall clock {@code from} in {@code phase}, the phase after it,
     * the 1-based mark it carries and how many jobs it releases.
     *
     * <p>Exact: within a segment the phase process is a homogeneous CTMC, and by the memoryless
     * property the residual holding time may be redrawn at a breakpoint, so the boundary is
     * crossed by advancing the clock and resampling under the new matrices.
     *
     * <p>NEITHER LABEL COSTS AN EXTRA DRAW: the competing transitions are accumulated
     * destination-major, then mark-minor, then batch-minor, so the running total after all
     * (mark, batch) pairs of a destination equals the aggregate total after that destination
     * and the winning destination is the one the unlabelled walk would choose for the same
     * uniform. Putting batch INSIDE mark is what makes a B = 1 BMMAPt reproduce the MMAPt
     * sample path for sample path, exactly as mark-inside-destination makes a K = 1 MMAPt
     * reproduce the MAPt.
     *
     * @param from   the wall-clock start
     * @param phase  the current phase
     * @param random the generator
     * @return {interval, next phase, mark, batch}; the interval is 0 once a non-cyclic horizon
     *         is exhausted
     */
    public double[] nextArrival(double from, int phase, Random random) {
        double elapsed = 0.0;
        double pos = from;
        int h = getNumberOfPhases();
        int K = d1kb.size();
        int B = getMaxBatchSize();
        int guard = 0;
        while (guard++ < 1000000) {
            int idx = getSegmentIndexAt(pos);
            if (idx < 0) {
                return new double[] {0.0, phase, 0, 0};
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
                    return new double[] {0.0, phase, 0, 0};
                }
                elapsed += toBoundary;
                pos += toBoundary;
                continue;
            }
            double holding = -Math.log(1.0 - random.nextDouble()) / total;
            if (holding >= toBoundary) {
                if (!cyclic && idx == d0.size() - 1) {
                    return new double[] {0.0, phase, 0, 0};
                }
                elapsed += toBoundary;
                pos += toBoundary;
                continue;
            }
            elapsed += holding;
            pos += holding;
            double u = random.nextDouble() * total;
            double cum = 0.0;
            for (int j = 0; j < h; j++) {
                for (int c = 0; c < K; c++) {
                    for (int b = 0; b < B; b++) {
                        cum += d1kb.get(c).get(b).get(idx).get(phase, j);
                        if (u < cum) {
                            return new double[] {elapsed, j, c + 1, b + 1};
                        }
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
                // Rounding left u at or past the total: fall back to the last block with any
                // mass out of this phase, as MMAPt does, rather than redrawing a holding time
                // and biasing the interval upwards.
                for (int j = h - 1; j >= 0; j--) {
                    for (int c = K - 1; c >= 0; c--) {
                        for (int b = B - 1; b >= 0; b--) {
                            if (d1kb.get(c).get(b).get(idx).get(phase, j) > 0.0) {
                                return new double[] {elapsed, j, c + 1, b + 1};
                            }
                        }
                    }
                }
            }
        }
        return new double[] {elapsed, phase, 0, 0};
    }

    private void checkMark(int k) {
        if (k < 1 || k > d1kb.size()) {
            throw new IllegalArgumentException("BMMAPt: mark index out of range: " + k);
        }
    }

    private void checkBatch(int b) {
        if (b < 1 || b > getMaxBatchSize()) {
            throw new IllegalArgumentException("BMMAPt: batch size out of range: " + b);
        }
    }

    /** Width-weighted mean of one block's n per-segment matrices. */
    private Matrix widthAverage(List<Matrix> segs) {
        int h = getNumberOfPhases();
        Matrix out = new Matrix(h, h, h * h);
        double total = getPeriod();
        for (int j = 0; j < segs.size(); j++) {
            double w = (breakpoints[j + 1] - breakpoints[j]) / total;
            out.addEq(w, segs.get(j));
        }
        return out;
    }

    /**
     * Per-block rates of the time-averaged process, through the canonical mmap_count_lambda so
     * that the stationary phase weighting is the one M3A uses rather than a second
     * implementation of it. The blocks must partition the aggregate D1bar.
     */
    private double[] markedRates(List<Matrix> blocks) {
        MatrixCell nominal = getTimeAverageProcess();
        MatrixCell cell = new MatrixCell();
        cell.set(0, nominal.get(0));
        cell.set(1, nominal.get(1));
        for (int c = 0; c < blocks.size(); c++) {
            cell.set(2 + c, blocks.get(c));
        }
        Matrix lk = jline.api.mam.Mmap_count_lambda.mmap_count_lambda(cell);
        double[] out = new double[blocks.size()];
        for (int c = 0; c < blocks.size(); c++) {
            out[c] = lk.get(0, c);
        }
        return out;
    }
}
