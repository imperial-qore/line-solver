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
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Time-inhomogeneous MARKED phase-type distribution (MPH_t).
 *
 * <p>An {@link MPH} whose three ingredients are functions of the wall clock: segment j covers
 * [breakpoints[j], breakpoints[j+1]) and carries an entry law alpha_j, a sub-generator S_j and
 * K exit vectors s_(1,j), ..., s_(K,j) satisfying the partition identity
 *
 * <pre>  sum_c s_(c,j) = -S_j e</pre>
 *
 * in every segment. This is to {@link PHt} what MPH is to {@link PH}.
 *
 * <p><b>It is STORED LOWERED to {@link MMAPt} form</b>, segment by segment, by
 *
 * <pre>  D0_j = S_j,   D1^(c)_j = s_(c,j) alpha_j</pre>
 *
 * exactly as a PHt is stored as its equivalent MAPt pair. One marked schedule walk therefore
 * serves both families and there is a single sn.proc shape behind
 * {@code ProcessType.isMarkedSchedule}. The original (alpha, S, exits) are kept on the object
 * for the getters and for serialization, so a round trip returns an MPHt and not the MMAPt it
 * lowers to.
 *
 * <p>The lowering makes each segment's arrivals RENEWAL within that segment: D1^(c)_j
 * factorises through alpha_j, so the phase after a completion does not depend on the phase
 * before it. Across a breakpoint the process is still time-varying, which is what separates an
 * MPH_t from an MPH.
 *
 * <p>At a Source the mark selects the class of the arriving job
 * ({@code Source.setMarkedArrival}, {@code sn.markidx}). As a SERVICE process it is sampled
 * for its duration and the mark is discarded, exactly as an MMAP is.
 *
 * <p>References: Q.-M. He, "The versatility of MMAP[K] and the MMAP[K]/G[K]/1 queue", Queueing
 * Systems 38(4), 2001; Y. M. Ko and J. Pender, "Diffusion limits for the (MAP_t/Ph_t/inf)^N
 * queueing network", Oper. Res. Lett. 45(3), 2017.
 */
public class MPHt extends ContinuousDistribution implements Serializable {

    private static final long serialVersionUID = 1L;

    /** Entry-law mass and the partition identity are checked to this absolute tolerance. */
    private static final double TOL = 1e-10;

    private final double[] breakpoints;
    private final List<Matrix> alpha;
    private final List<Matrix> subgen;
    /** exit.get(c).get(j) is the exit vector of mark c+1 in segment j. */
    private final List<List<Matrix>> exit;
    private final boolean cyclic;
    private final MMAPt lowered;

    /**
     * Creates an MPH_t with a piecewise-constant marked phase-type schedule.
     *
     * @param breakpoints strictly increasing segment boundaries, length n+1
     * @param alpha       per-segment 1-by-h entry laws, length n
     * @param S           per-segment h-by-h sub-generators, length n
     * @param exit        the K mark exits, each a list of n per-segment h-by-1 vectors
     * @param cyclic      whether the schedule repeats with the horizon as period
     */
    public MPHt(double[] breakpoints, List<Matrix> alpha, List<Matrix> S,
                List<List<Matrix>> exit, boolean cyclic) {
        super("MPHt", 0, new Pair<Double, Double>(0.0, GlobalConstants.Inf));
        if (breakpoints == null || alpha == null || S == null || exit == null
                || alpha.isEmpty() || exit.isEmpty() || S.size() != alpha.size()) {
            throw new IllegalArgumentException(
                    "MPHt: alpha, S and a non-empty list of mark exits are all required, with "
                            + "alpha and S of equal length");
        }
        int n = alpha.size();
        int K = exit.size();
        for (int c = 0; c < K; c++) {
            if (exit.get(c) == null || exit.get(c).size() != n) {
                throw new IllegalArgumentException(
                        "MPHt: mark exit " + (c + 1) + " has "
                                + (exit.get(c) == null ? 0 : exit.get(c).size())
                                + " segments against " + n + " in alpha; every mark must be "
                                + "defined on the whole schedule");
            }
        }
        int h = alpha.get(0).getNumCols();
        for (int j = 0; j < n; j++) {
            Matrix a = alpha.get(j);
            Matrix s = S.get(j);
            if (a.getNumRows() != 1 || a.getNumCols() != h) {
                throw new IllegalArgumentException(
                        "MPHt: every alpha must be a 1-by-" + h + " row vector; segment "
                                + (j + 1) + " differs");
            }
            if (s.getNumRows() != h || s.getNumCols() != h) {
                throw new IllegalArgumentException(
                        "MPHt: every S must be square of order " + h + "; segment " + (j + 1)
                                + " differs");
            }
            double amass = 0.0;
            for (int i = 0; i < h; i++) {
                if (a.get(0, i) < 0.0) {
                    throw new IllegalArgumentException(
                            "MPHt: alpha must be non-negative in segment " + (j + 1));
                }
                amass += a.get(0, i);
            }
            if (Math.abs(amass - 1.0) > TOL) {
                throw new IllegalArgumentException(
                        "MPHt: alpha must sum to one in every segment; segment " + (j + 1)
                                + " sums to " + amass);
            }
            for (int i = 0; i < h; i++) {
                if (s.get(i, i) >= 0.0) {
                    throw new IllegalArgumentException(
                            "MPHt: the diagonal of S must be negative; segment " + (j + 1)
                                    + " entry (" + i + "," + i + ") is " + s.get(i, i));
                }
                for (int q = 0; q < h; q++) {
                    if (i != q && s.get(i, q) < 0.0) {
                        throw new IllegalArgumentException(
                                "MPHt: off-diagonal entries of S must be non-negative in segment "
                                        + (j + 1));
                    }
                }
            }
            for (int c = 0; c < K; c++) {
                Matrix sk = exit.get(c).get(j);
                if (sk == null || sk.getNumRows() != h || sk.getNumCols() != 1) {
                    throw new IllegalArgumentException(
                            "MPHt: exit vector of mark " + (c + 1) + " in segment " + (j + 1)
                                    + " must be " + h + "-by-1");
                }
                for (int i = 0; i < h; i++) {
                    if (sk.get(i, 0) < 0.0) {
                        throw new IllegalArgumentException(
                                "MPHt: exit vectors must be non-negative; mark " + (c + 1)
                                        + " of segment " + (j + 1) + " is not");
                    }
                }
            }
            for (int i = 0; i < h; i++) {
                double rowsum = 0.0;
                for (int q = 0; q < h; q++) {
                    rowsum += s.get(i, q);
                }
                double marked = 0.0;
                for (int c = 0; c < K; c++) {
                    marked += exit.get(c).get(j).get(i, 0);
                }
                if (Math.abs(marked + rowsum) > TOL) {
                    throw new IllegalArgumentException(
                            "MPHt: the exit vectors must partition the absorption rate of S; in "
                                    + "segment " + (j + 1) + " phase " + i + " they sum to "
                                    + marked + " against an absorption rate of " + (-rowsum));
                }
            }
        }

        // The lowering. Its constructor re-checks the generator property and the common
        // support, so a schedule that passes here is a valid MMAPt by construction.
        List<Matrix> d0 = new ArrayList<Matrix>();
        for (int j = 0; j < n; j++) {
            d0.add(S.get(j).copy());
        }
        List<List<Matrix>> d1k = new ArrayList<List<Matrix>>();
        for (int c = 0; c < K; c++) {
            List<Matrix> blocks = new ArrayList<Matrix>();
            for (int j = 0; j < n; j++) {
                Matrix block = new Matrix(h, h, h * h);
                Matrix sk = exit.get(c).get(j);
                Matrix a = alpha.get(j);
                for (int i = 0; i < h; i++) {
                    double si = sk.get(i, 0);
                    if (si == 0.0) {
                        continue;
                    }
                    for (int q = 0; q < h; q++) {
                        double v = si * a.get(0, q);
                        if (v != 0.0) {
                            block.set(i, q, v);
                        }
                    }
                }
                blocks.add(block);
            }
            d1k.add(blocks);
        }
        this.lowered = new MMAPt(breakpoints, d0, d1k, cyclic);

        this.breakpoints = breakpoints.clone();
        this.alpha = new ArrayList<Matrix>();
        this.subgen = new ArrayList<Matrix>();
        for (int j = 0; j < n; j++) {
            this.alpha.add(alpha.get(j).copy());
            this.subgen.add(S.get(j).copy());
        }
        this.exit = new ArrayList<List<Matrix>>();
        for (int c = 0; c < K; c++) {
            List<Matrix> col = new ArrayList<Matrix>();
            for (int j = 0; j < n; j++) {
                col.add(exit.get(c).get(j).copy());
            }
            this.exit.add(col);
        }
        this.cyclic = cyclic;
        this.mean = lowered.getMean();
        this.immediate = false;
    }

    /** Creates a cyclic MPH_t. */
    public MPHt(double[] breakpoints, List<Matrix> alpha, List<Matrix> S,
                List<List<Matrix>> exit) {
        this(breakpoints, alpha, S, exit, true);
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

    /**
     * Per-segment exit vectors of one mark.
     *
     * @param k the 1-based mark index
     * @return the n per-segment exit vectors
     */
    public List<Matrix> getExitSegments(int k) {
        if (k < 1 || k > exit.size()) {
            throw new IllegalArgumentException("MPHt: mark index out of range: " + k);
        }
        return new ArrayList<Matrix>(exit.get(k - 1));
    }

    /** All exit vectors, mark-major: get(c).get(j) is mark c+1 in segment j. */
    public List<List<Matrix>> getExitSegments() {
        List<List<Matrix>> out = new ArrayList<List<Matrix>>();
        for (int c = 0; c < exit.size(); c++) {
            out.add(new ArrayList<Matrix>(exit.get(c)));
        }
        return out;
    }

    public int getNumberOfTypes() {
        return exit.size();
    }

    public boolean isCyclic() {
        return cyclic;
    }

    public int getNumSegments() {
        return alpha.size();
    }

    public int getNumberOfPhases() {
        return alpha.get(0).getNumCols();
    }

    /** Horizon length, which is the period when cyclic. */
    public double getPeriod() {
        return breakpoints[breakpoints.length - 1] - breakpoints[0];
    }

    /** Index of the segment in force at t, or -1 past a non-cyclic horizon. */
    public int getSegmentIndexAt(double t) {
        return lowered.getSegmentIndexAt(t);
    }

    /**
     * The equivalent MMAP_t, which is the form this process is stored and sampled in.
     *
     * @return the lowered marked schedule
     */
    public MMAPt toMMAPt() {
        return lowered;
    }

    /** The UNMARKED schedule of the duration, PHt(alpha, S). */
    public PHt toPHt() {
        return new PHt(breakpoints.clone(), new ArrayList<Matrix>(alpha),
                new ArrayList<Matrix>(subgen), cyclic);
    }

    /** The marginal schedule of one mark; see {@link MMAPt#toMAPts(int)}. */
    public MAPt toMAPts(int k) {
        return lowered.toMAPts(k);
    }

    /** Width-weighted average pair over the horizon, {D0bar, D1bar} of the lowered process. */
    public MatrixCell getTimeAverageProcess() {
        return lowered.getTimeAverageProcess();
    }

    /** Arrival rate of the time-averaged aggregate MAP. */
    public double getTimeAverageRate() {
        return lowered.getTimeAverageRate();
    }

    /** Per-mark arrival rates of the time-averaged process; see {@link MMAPt}. */
    public double[] getTimeAverageMarkRates() {
        return lowered.getTimeAverageMarkRates();
    }

    @Override
    public double getMean() {
        return lowered.getMean();
    }

    @Override
    public double getRate() {
        return lowered.getRate();
    }

    /** NaN: an MPH_t is time-varying, so no i.i.d. interval law summarises it. */
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

    /**
     * The MMAPt-shaped cell of the lowering, which is what reaches sn.proc.
     *
     * <p>The procid stays MPHT, so the featset gate and the wire type still tell the two
     * families apart; only the representation is shared.
     */
    @Override
    public MatrixCell getProcess() {
        return lowered.getProcess();
    }

    /** Restarts the sample path at the schedule start. */
    public void resetSampleClock() {
        lowered.resetSampleClock();
    }

    /** The 1-based mark of the interval last returned by {@link #sample(int, Random)}. */
    public int getLastMark() {
        return lowered.getLastMark();
    }

    @Override
    public double[] sample(int n, Random random) {
        return lowered.sample(n, random);
    }

    @Override
    public double[] sample(int n) {
        return lowered.sample(n);
    }

    /** Time to the next completion, the phase after it and its mark; see {@link MMAPt}. */
    public double[] nextArrival(double from, int phase, Random random) {
        return lowered.nextArrival(from, phase, random);
    }
}
