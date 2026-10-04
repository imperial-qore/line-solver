/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * A Marked Phase-Type distribution (MPH).
 *
 * <p>A PH whose absorption is LABELLED: the chain runs on h transient phases with entry
 * law alpha and sub-generator S, and absorption happens through one of K exit vectors
 * s_1, ..., s_K, so a completion carries both a duration and the mark of the exit that
 * produced it. The partition identity is
 *
 * <pre>  sum_k s_k = -S e</pre>
 *
 * which is what makes the K exits account for exactly the absorption the sub-generator
 * leaves, no more and no less.
 *
 * <p><b>An MPH is the RENEWAL special case of an MMAP.</b> Lowering it by
 *
 * <pre>  D0 = S,   D1k = s_k alpha,   D1 = sum_k D1k = (-S e) alpha</pre>
 *
 * gives a marked Markovian arrival process whose successive intervals are independent
 * and PH(alpha, S) distributed, because D1k(i, j) factorises through alpha and so the
 * phase after an event does not depend on the phase before it. That is exactly the
 * cell this class hands to {@link MarkedMAP}, so every consumer of the M3A layout
 * {D0, D1_agg, D11, ..., D1K} serves an MPH unchanged and no separate sampler,
 * serializer or state machinery is needed.
 *
 * <p>It nevertheless carries a {@link jline.lang.constant.ProcessType#MPH} of its own
 * rather than aliasing MMAP, so a solver that cannot honour a marked renewal process
 * refuses it by name instead of inheriting MMAP's support. Because this class extends
 * MarkedMAP, any dispatch that tests {@code instanceof MarkedMAP} must test MPH FIRST
 * or it will tag a marked PH as an ordinary MMAP; that is the same ordering trap
 * {@link BMAP} already carries.
 *
 * <p>At a Source the mark selects the class of the arriving job
 * ({@code Source.setMarkedArrival}, {@code sn.markidx}). As a SERVICE process it is
 * sampled for its duration and the mark is discarded, exactly as an MMAP is.
 *
 * <p>Reference: Q.-M. He, "The versatility of MMAP[K] and the MMAP[K]/G[K]/1 queue",
 * Queueing Systems 38(4), 2001.
 */
public class MPH extends MarkedMAP implements Serializable {

    private static final long serialVersionUID = 1L;

    /** Row sums, entry law and partition identity are checked to this absolute tolerance. */
    private static final double TOL = 1e-10;

    private final Matrix alpha;
    private final Matrix subgen;
    private final List<Matrix> exit;

    /**
     * Creates a marked phase-type distribution.
     *
     * @param alpha  1-by-h entry law, non-negative and summing to one
     * @param S      h-by-h sub-generator: non-negative off-diagonal, negative diagonal
     * @param exit   the K exit vectors, each h-by-1 and non-negative, with sum_k s_k = -S e
     */
    public MPH(Matrix alpha, Matrix S, List<Matrix> exit) {
        super(toM3ACell(alpha, S, exit), exit == null ? 0 : exit.size());
        this.name = "MPH";
        this.alpha = alpha.copy();
        this.subgen = S.copy();
        this.exit = new ArrayList<Matrix>();
        for (int k = 0; k < exit.size(); k++) {
            this.exit.add(exit.get(k).copy());
        }
        this.immediate = false;
    }

    /**
     * Lowers (alpha, S, {s_k}) to the M3A cell {D0, D1_agg, D11, ..., D1K}, validating as
     * it goes. Static because it runs before the superclass constructor.
     */
    private static MatrixCell toM3ACell(Matrix alpha, Matrix S, List<Matrix> exit) {
        if (alpha == null || S == null || exit == null || exit.isEmpty()) {
            throw new IllegalArgumentException(
                    "MPH: alpha, S and a non-empty list of exit vectors are all required");
        }
        if (alpha.getNumRows() != 1) {
            throw new IllegalArgumentException("MPH: alpha must be a row vector");
        }
        int h = alpha.getNumCols();
        if (S.getNumRows() != h || S.getNumCols() != h) {
            throw new IllegalArgumentException(
                    "MPH: S must be square of order " + h + " to match alpha");
        }
        double amass = 0.0;
        for (int i = 0; i < h; i++) {
            double a = alpha.get(0, i);
            if (a < 0.0) {
                throw new IllegalArgumentException("MPH: alpha must be non-negative");
            }
            amass += a;
        }
        if (Math.abs(amass - 1.0) > TOL) {
            throw new IllegalArgumentException(
                    "MPH: alpha must sum to one; it sums to " + amass
                            + ". A defective entry law would put mass on an instantaneous "
                            + "completion that carries no mark.");
        }
        for (int i = 0; i < h; i++) {
            if (S.get(i, i) >= 0.0) {
                throw new IllegalArgumentException(
                        "MPH: the diagonal of S must be negative; entry (" + i + "," + i
                                + ") is " + S.get(i, i));
            }
            for (int j = 0; j < h; j++) {
                if (i != j && S.get(i, j) < 0.0) {
                    throw new IllegalArgumentException(
                            "MPH: off-diagonal entries of S must be non-negative; entry ("
                                    + i + "," + j + ") is " + S.get(i, j));
                }
            }
        }
        int K = exit.size();
        for (int k = 0; k < K; k++) {
            Matrix sk = exit.get(k);
            if (sk == null || sk.getNumRows() != h || sk.getNumCols() != 1) {
                throw new IllegalArgumentException(
                        "MPH: exit vector " + (k + 1) + " must be " + h + "-by-1");
            }
            for (int i = 0; i < h; i++) {
                if (sk.get(i, 0) < 0.0) {
                    throw new IllegalArgumentException(
                            "MPH: exit vector " + (k + 1) + " must be non-negative");
                }
            }
        }
        // sum_k s_k = -S e, phase by phase. This is the identity that makes D0 + sum_k D1k
        // a generator, so a violation here is a process that either loses probability mass
        // or manufactures it, not a rounding question.
        for (int i = 0; i < h; i++) {
            double rowsum = 0.0;
            for (int j = 0; j < h; j++) {
                rowsum += S.get(i, j);
            }
            double marked = 0.0;
            for (int k = 0; k < K; k++) {
                marked += exit.get(k).get(i, 0);
            }
            if (Math.abs(marked + rowsum) > TOL) {
                throw new IllegalArgumentException(
                        "MPH: the exit vectors must partition the absorption rate of S; in phase "
                                + i + " they sum to " + marked + " against an absorption rate of "
                                + (-rowsum) + ". Every way of leaving the transient phases must "
                                + "carry exactly one mark.");
            }
        }
        boolean anyExit = false;
        for (int k = 0; k < K && !anyExit; k++) {
            for (int i = 0; i < h; i++) {
                if (exit.get(k).get(i, 0) > 0.0) {
                    anyExit = true;
                    break;
                }
            }
        }
        if (!anyExit) {
            throw new IllegalArgumentException(
                    "MPH: every exit vector is zero, so absorption can never occur");
        }

        MatrixCell cell = new MatrixCell();
        cell.set(0, S.copy());
        Matrix agg = new Matrix(h, h, h * h);
        List<Matrix> blocks = new ArrayList<Matrix>();
        for (int k = 0; k < K; k++) {
            Matrix d1k = new Matrix(h, h, h * h);
            for (int i = 0; i < h; i++) {
                double si = exit.get(k).get(i, 0);
                if (si == 0.0) {
                    continue;
                }
                for (int j = 0; j < h; j++) {
                    double v = si * alpha.get(0, j);
                    if (v != 0.0) {
                        d1k.set(i, j, v);
                        agg.set(i, j, agg.get(i, j) + v);
                    }
                }
            }
            blocks.add(d1k);
        }
        cell.set(1, agg);
        for (int k = 0; k < K; k++) {
            cell.set(2 + k, blocks.get(k));
        }
        return cell;
    }

    /** The 1-by-h entry law. */
    public Matrix getAlpha() {
        return alpha.copy();
    }

    /** The h-by-h sub-generator. */
    public Matrix getSubgenerator() {
        return subgen.copy();
    }

    /**
     * The h-by-1 exit vector of mark k.
     *
     * @param k the 1-based mark index
     * @return the exit vector
     */
    public Matrix getExitVector(int k) {
        if (k < 1 || k > exit.size()) {
            throw new IllegalArgumentException("MPH: mark index out of range: " + k);
        }
        return exit.get(k - 1).copy();
    }

    /** The exit vectors, ordered by mark. */
    public List<Matrix> getExitVectors() {
        List<Matrix> out = new ArrayList<Matrix>();
        for (int k = 0; k < exit.size(); k++) {
            out.add(exit.get(k).copy());
        }
        return out;
    }

    /** Number of transient phases. Widens Markovian.getNumberOfPhases, which returns long. */
    @Override
    public long getNumberOfPhases() {
        return alpha.getNumCols();
    }

    /**
     * The UNMARKED phase-type distribution of the duration, PH(alpha, S).
     *
     * <p>Hiding the marks leaves the interval law untouched, which is what separates an
     * MPH from a general MMAP: the aggregate of an MMAP is a MAP, and only here is it a
     * renewal process with a phase-type marginal.
     *
     * @return the duration's PH law
     */
    public PH toPH() {
        return new PH(alpha.copy(), subgen.copy());
    }

    /**
     * The probability that a completion carries mark k, alpha (-S)^-1 s_k.
     *
     * <p>Successive completions are independent, so this is both the long-run fraction of
     * mark-k completions and the probability of any one of them.
     *
     * @return the K mark probabilities, ordered by mark
     */
    public double[] getMarkProbabilities() {
        int h = alpha.getNumCols();
        Matrix inv = subgen.neg().inv();
        Matrix tau = alpha.mult(inv);
        double[] out = new double[exit.size()];
        for (int k = 0; k < exit.size(); k++) {
            double p = 0.0;
            for (int i = 0; i < h; i++) {
                p += tau.get(0, i) * exit.get(k).get(i, 0);
            }
            out[k] = p;
        }
        return out;
    }
}
