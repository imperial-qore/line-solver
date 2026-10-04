/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.RandomManager;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.io.Serializable;
import java.util.Random;

import static jline.api.mam.Mmap_timereverse.mmap_timereverse;

/**
 * A Marked Markovian Arrival Process
 */
public class MarkedMAP extends Marked implements Serializable {

    public MarkedMAP() {
        super("MarkedMAP", 0);
    }

    /**
     * Constructs a MarkedMAP from a MatrixCell, normalizing it to the M3A
     * layout {D0, D1_agg, D11, ..., D1K} used throughout the mam API
     * (mmap_sample, mmap_normalize, mmap_super) and by Marked.D(1,k).
     *
     * Accepted input layouts:
     * - {D0, D1_agg, D11, ..., D1K} (M3A): detected when cell 1 equals the sum
     *   of cells 2..end; stored as-is.
     * - {D0, D11, ..., D1K} (per-type, no aggregate): the aggregate is inserted
     *   at index 1.
     * - {D0} or {D0, D1}: unmarked (K=0); stored as-is.
     *
     * The detection is a sum check with relative tolerance; for ambiguous
     * inputs use {@link #MarkedMAP(MatrixCell, int)} with an explicit type
     * count, mirroring the MATLAB MarkedMAP(D,K) constructor.
     */
    public MarkedMAP(MatrixCell mmap) {
        this(toM3A(mmap), 0, true);
    }

    /**
     * Constructs a MarkedMAP with an explicit number of arrival types,
     * mirroring MATLAB MarkedMAP(D,K): K == size-2 means the cell is already
     * in M3A layout {D0, D1_agg, D11..D1K}; K == size-1 means the cell is
     * {D0, D11..D1K} and the aggregate is inserted at index 1.
     *
     * @param mmap the process matrices
     * @param K    the number of arrival types (marks)
     */
    public MarkedMAP(MatrixCell mmap, int K) {
        this(explicitToM3A(mmap, K), 0, true);
    }

    private MarkedMAP(MatrixCell m3a, int unusedDisambiguator, boolean fromM3A) {
        super("MMAP", m3a.size());
        setProcess(m3a);
        Matrix D0 = m3a.get(0);
        Matrix D1 = m3a.size() > 1 ? m3a.get(1) : null;
        this.setParam(1, "D0", D0);
        if (D1 != null) {
            this.setParam(2, "D1", D1);
        }
        for (int i = 2; i < m3a.size(); i++) {
            this.setParam(i + 1, "D1" + i, m3a.get(i));
        }
        this.nPhases = D0.getNumCols();
    }

    /**
     * Number of arrival types (marks) K under the M3A layout.
     *
     * @return K = size-2, or 0 for an unmarked {D0[,D1]} process
     */
    public int getNumberOfTypes() {
        return Math.max(0, getProcess().size() - 2);
    }

    /** Normalizes a constructor input cell to the M3A layout (see javadoc). */
    private static MatrixCell toM3A(MatrixCell mmap) {
        if (mmap == null || mmap.size() <= 2) {
            return mmap;
        }
        Matrix sumTail = new Matrix(mmap.get(0).getNumRows(), mmap.get(0).getNumCols());
        for (int k = 2; k < mmap.size(); k++) {
            sumTail.addEq(1.0, mmap.get(k));
        }
        if (matrixApproxEqual(mmap.get(1), sumTail)) {
            return mmap; // already M3A
        }
        return insertAggregate(mmap);
    }

    private static MatrixCell explicitToM3A(MatrixCell mmap, int K) {
        if (K == mmap.size() - 2) {
            return mmap;
        } else if (K == mmap.size() - 1) {
            return insertAggregate(mmap);
        }
        throw new IllegalArgumentException(
                "Inconsistency between the number of types K and the number of supplied D matrices.");
    }

    /** Builds {D0, sum(D1k), D11..D1K} from a {D0, D11..D1K} cell. */
    private static MatrixCell insertAggregate(MatrixCell mmap) {
        MatrixCell out = new MatrixCell();
        out.set(0, mmap.get(0));
        Matrix agg = new Matrix(mmap.get(0).getNumRows(), mmap.get(0).getNumCols());
        for (int k = 1; k < mmap.size(); k++) {
            agg.addEq(1.0, mmap.get(k));
        }
        out.set(1, agg);
        for (int k = 1; k < mmap.size(); k++) {
            out.set(1 + k, mmap.get(k));
        }
        return out;
    }

    private static boolean matrixApproxEqual(Matrix a, Matrix b) {
        if (a.getNumRows() != b.getNumRows() || a.getNumCols() != b.getNumCols()) {
            return false;
        }
        for (int i = 0; i < a.getNumRows(); i++) {
            for (int j = 0; j < a.getNumCols(); j++) {
                double av = a.get(i, j);
                double bv = b.get(i, j);
                if (Math.abs(av - bv) > 1e-10 * Math.max(1.0, Math.abs(av))) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Normalizes the MarkedMAP so that D0+D1_agg rows form a proper
     * infinitesimal generator, with the aggregate at index 1 recomputed from
     * the per-type matrices (M3A layout). Delegates to the mam API
     * mmap_normalize; falls back to a {D0,D1} MAP-style normalization for an
     * unmarked process.
     */
    public void normalize() {
        MatrixCell proc = getProcess();
        if (proc.size() > 2) {
            MatrixCell normalized = jline.api.mam.Mmap_normalize.mmap_normalize(proc);
            setProcess(normalized);
            this.setParam(1, "D0", normalized.get(0));
            this.setParam(2, "D1", normalized.get(1));
            for (int i = 2; i < normalized.size(); i++) {
                this.setParam(i + 1, "D1" + i, normalized.get(i));
            }
            return;
        }
        // Unmarked {D0[,D1]}: clip negatives and rebalance the D0 diagonal
        Matrix D0 = D(0).copy();
        Matrix D1 = proc.size() > 1 ? proc.get(1) : null;
        for (int i = 0; i < nPhases; i++) {
            for (int j = 0; j < nPhases; j++) {
                if (i != j && D0.get(i, j) < 0) {
                    D0.set(i, j, 0.0);
                }
                if (D1 != null && D1.get(i, j) < 0) {
                    D1.set(i, j, 0.0);
                }
            }
        }
        for (int i = 0; i < nPhases; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < nPhases; j++) {
                if (i != j) {
                    rowSum += D0.get(i, j);
                }
                if (D1 != null) {
                    rowSum += D1.get(i, j);
                }
            }
            D0.set(i, i, -rowSum);
        }
        proc.set(0, D0);
        setProcess(proc);
    }

    @Override
    public double[] sample(int n) {
        return this.sample(n, RandomManager.getThreadRandomAsRandom());
    }

    /**
     * Generates samples from the MarkedMAP using the specified random generator.
     * Returns inter-arrival times with corresponding marks.
     *
     * @param n      the number of samples to generate
     * @param random the random number generator to use
     * @return array of inter-arrival times (marks are not returned in this basic implementation)
     */
    @Override
    public double[] sample(int n, Random random) {
        try {
            // see _kb/01-model-classes.md (Java process-construction notes) for rationale
            Matrix D0 = D(0);
            Matrix D1_agg = getProcess().get(1);
            return jline.api.mam.Map_sample.map_sample(D0, D1_agg, n, random);
        } catch (Exception e) {
            // Fallback: generate exponential samples with the mean rate
            double[] samples = new double[n];
            double mean = getMean();
            for (int i = 0; i < n; i++) {
                samples[i] = -mean * Math.log(random.nextDouble());
            }
            return samples;
        }
    }

    public MarkedMAP toTimeReversed() {
        MatrixCell reversed = mmap_timereverse(this.process);
        return new MarkedMAP(reversed);
    }

    /**
     * Aggregate MAP over all marks: MAP(D0, D1_agg). Mirrors MATLAB
     * MarkedMAP.toMAP.
     *
     * @return the aggregate (unmarked) MAP
     */
    public MAP toMAP() {
        return new MAP(D(0).copy(), getProcess().get(1).copy());
    }

    /**
     * Marginal MAP of mark k: arrivals fire only on D1k, while the other
     * marks' transitions become hidden phase changes, i.e. MAP(D0 + D1_agg
     * - D1k, D1k). Mirrors MATLAB MarkedMAP.toMAPs for a single type.
     *
     * @param k the mark index (1-based)
     * @return the marginal MAP of mark k
     */
    public MAP toMAPs(int k) {
        if (k < 1 || k > getNumberOfTypes()) {
            throw new IllegalArgumentException("Mark index out of range: " + k);
        }
        Matrix D1k = getProcess().get(1 + k);
        Matrix Df = D(0).copy();
        Df.addEq(1.0, getProcess().get(1));
        Df.addEq(-1.0, D1k);
        return new MAP(Df, D1k.copy());
    }

}
