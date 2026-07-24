/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.qsys;

import java.util.ArrayList;
import java.util.Iterator;

import org.ejml.data.DMatrixSparseCSC;
import org.ejml.data.DMatrixSparseTriplet;
import org.ejml.ops.DConvertMatrixStruct;

import jline.GlobalConstants;
import jline.api.mc.Ctmc_solve;
import jline.api.mc.Dtmc_solve;
import jline.io.InputOutput;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixEntry;

/**
 * Analyzes a BMAP/M/1 queue by the matrix-analytic (M/G/1-type) method.
 *
 * <p>Beyond the usual performance measures the result exposes the intermediate
 * matrix-analytic quantities themselves (theta, lambda, rho, the uniformization
 * constant q, the randomized blocks A0/A1/B0/Bk, A, alpha, the matrix G, the
 * drift, the measured decay rate and the level probabilities), so that the
 * algorithm can be inspected and taught rather than only its output.</p>
 *
 * <p>Port of MATLAB qsys_bmapm1.m.</p>
 */
public final class Qsys_bmapm1 {

    /** Default maximum functional iterations for G. */
    public static final int DEFAULT_MAX_ITER = 10000;
    /** Default convergence tolerance for G. */
    public static final double DEFAULT_TOLERANCE = 1e-12;
    /** Default relative truncation target for the level distribution. */
    public static final double DEFAULT_TAIL_TOLERANCE = 1e-10;

    private Qsys_bmapm1() {
    }

    /**
     * Analyzes a BMAP/M/1 queue with the default uniformization constant and
     * adaptive level truncation.
     *
     * @param D  BMAP matrices {D0, D1, ..., DK}; D0 carries the hidden transitions,
     *           Dk (k &gt;= 1) the transitions that release a batch of k customers
     * @param mu exponential service rate
     * @return the matrix-analytic result
     */
    public static QsysBmapM1Result qsys_bmapm1(Matrix[] D, double mu) {
        return qsys_bmapm1(D, mu, Double.NaN, DEFAULT_MAX_ITER, DEFAULT_TOLERANCE, -1,
                DEFAULT_TAIL_TOLERANCE);
    }

    /**
     * Analyzes a BMAP/M/1 queue.
     *
     * @param D             BMAP matrices {D0, D1, ..., DK}
     * @param mu            exponential service rate
     * @param uniformization uniformization constant q used to randomize the generator
     *                      into a discrete-time M/G/1-type chain; it must dominate every
     *                      total outflow rate. Pass NaN to use max_i(-D0(i,i)) + mu
     * @param maxIter       maximum functional iterations for G
     * @param tolerance     convergence tolerance for G
     * @param maxLevel      fixed level truncation for the queue-length distribution;
     *                      non-positive selects the adaptive rule
     * @param tailTolerance relative truncation target for the level distribution
     * @return the matrix-analytic result
     */
    public static QsysBmapM1Result qsys_bmapm1(Matrix[] D, double mu, double uniformization,
                                               int maxIter, double tolerance, int maxLevel,
                                               double tailTolerance) {
        // Validate inputs
        if (D == null || D.length < 2) {
            InputOutput.line_error("qsys_bmapm1",
                    "The BMAP must be given as a cell array {D0,D1,...,DK} with at least D0 and D1.");
        }
        int V = D[0].getNumRows();
        for (int k = 0; k < D.length; k++) {
            Matrix Dk = D[k];
            if (Dk == null || Dk.getNumRows() != V || Dk.getNumCols() != V || !Dk.isFinite()) {
                InputOutput.line_error("qsys_bmapm1", String.format(
                        "BMAP matrix D{%d} must be a finite %dx%d matrix.", k + 1, V, V));
            }
            if (k >= 1 && minEntry(Dk) < -GlobalConstants.FineTol) {
                InputOutput.line_error("qsys_bmapm1", String.format(
                        "BMAP arrival matrix D{%d} must be non-negative.", k + 1));
            }
        }
        Matrix Dsum = new Matrix(V, V);
        Dsum.zero();
        for (int k = 0; k < D.length; k++) {
            Dsum = Dsum.add(1.0, D[k]);
        }
        for (int i = 0; i < V; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < V; j++) {
                rowSum += Dsum.get(i, j);
            }
            if (Math.abs(rowSum) > Math.sqrt(GlobalConstants.FineTol)) {
                InputOutput.line_error("qsys_bmapm1",
                        "BMAP matrices are inconsistent: sum_k D_k must have zero row sums.");
            }
        }
        if (!Double.isFinite(mu) || mu <= 0) {
            InputOutput.line_error("qsys_bmapm1", "The service rate mu must be a finite positive scalar.");
        }

        int K = D.length - 1;

        // Arrival characterization
        Matrix theta = Ctmc_solve.ctmc_solve(Dsum);
        theta = asRow(theta);
        Matrix sumKDk = new Matrix(V, V);
        sumKDk.zero();
        for (int k = 1; k <= K; k++) {
            sumKDk = sumKDk.add((double) k, D[k]);
        }
        double lambda = 0.0;
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                lambda += theta.get(0, i) * sumKDk.get(i, j);
            }
        }
        double rho = lambda / mu;

        // Randomization (uniformization) into a discrete-time M/G/1-type chain
        double maxOutflow = 0.0;
        for (int i = 0; i < V; i++) {
            maxOutflow = Math.max(maxOutflow, -D[0].get(i, i));
        }
        double q = Double.isNaN(uniformization) ? (maxOutflow + mu) : uniformization;
        if (q < maxOutflow + mu - GlobalConstants.FineTol) {
            InputOutput.line_error("qsys_bmapm1", String.format(
                    "The uniformization constant q = %g does not dominate the total outflow rate %g; "
                    + "the randomized chain would have negative entries.", q, maxOutflow + mu));
        }

        Matrix A0 = Matrix.scaleMult(Matrix.eye(V), mu / q);          // level down by one
        Matrix A1 = Matrix.scaleMult(D[0].add(-mu, Matrix.eye(V)), 1.0 / q).add(1.0, Matrix.eye(V));
        Matrix B0 = Matrix.scaleMult(D[0], 1.0 / q).add(1.0, Matrix.eye(V));
        Matrix[] Bk = new Matrix[K + 1];
        for (int k = 1; k <= K; k++) {
            Bk[k] = Matrix.scaleMult(D[k], 1.0 / q);                  // level up by k
        }

        Matrix A = A0.add(1.0, A1);
        for (int k = 1; k <= K; k++) {
            A = A.add(1.0, Bk[k]);
        }
        Matrix alpha = asRow(Dtmc_solve.dtmc_solve(A));

        // Matrix G: minimal non-negative solution of the M/G/1-type equation
        Matrix G = new Matrix(V, V);
        G.zero();
        boolean converged = false;
        double lastChange = Double.POSITIVE_INFINITY;
        for (int iter = 0; iter < maxIter; iter++) {
            Matrix Gpow = G.copy();
            Matrix Gnew = A0.add(1.0, A1.mult(G));
            for (int k = 1; k <= K; k++) {
                Gpow = Gpow.mult(G);                                  // G^(k+1)
                Gnew = Gnew.add(1.0, Bk[k].mult(Gpow));
            }
            lastChange = maxAbsDiff(Gnew, G);
            G = Gnew;
            if (lastChange < tolerance) {
                converged = true;
                break;
            }
        }
        if (!converged) {
            InputOutput.line_warning("qsys_bmapm1", String.format(
                    "The functional iteration for G did not converge to %g in %d iterations "
                    + "(last change %g). The queue may be unstable.", tolerance, maxIter, lastChange));
        }

        // Stability drift
        Matrix upDrift = new Matrix(V, V);
        upDrift.zero();
        for (int k = 1; k <= K; k++) {
            upDrift = upDrift.add((double) k, Bk[k]);
        }
        double drift = rowWeightedSum(alpha, upDrift) - rowWeightedSum(alpha, A0);

        // see _kb/03-api-layer.md for rationale
        int levelMax;
        Matrix levelProb;
        double truncError;
        if (maxLevel > 0) {
            levelMax = maxLevel;
            levelProb = solveLevels(D, mu, V, K, levelMax);
            truncError = levelTailError(levelProb, levelMax);
        } else {
            levelMax = Math.max(50, (int) Math.ceil(20.0 / Math.max(1.0 - Math.min(rho, 0.999), 2.2e-16)));
            while (true) {
                levelProb = solveLevels(D, mu, V, K, levelMax);
                truncError = levelTailError(levelProb, levelMax);
                if (truncError <= tailTolerance || ((long) 2 * levelMax + 1) * V > 2e5) {
                    break;
                }
                levelMax = 2 * levelMax;
            }
            if (truncError > tailTolerance) {
                InputOutput.line_warning("qsys_bmapm1", String.format(
                        "The level distribution did not reach the requested accuracy: residual %.3e > "
                        + "TailTolerance %.3e at level %d.", truncError, tailTolerance, levelMax));
            }
        }

        int nLevels = levelProb.getNumRows();
        double[] levelMass = new double[nLevels];
        for (int i = 0; i < nLevels; i++) {
            for (int j = 0; j < V; j++) {
                levelMass[i] += levelProb.get(i, j);
            }
        }
        // Measured decay rate: the ratio settles geometrically, so read it where the
        // mass is still numerically meaningful rather than at the truncation boundary.
        int usable = -1;
        for (int i = nLevels - 1; i >= 0; i--) {
            if (levelMass[i] > 1e-12) {
                usable = i;
                break;
            }
        }
        double decayRate;
        if (usable < 2) {
            decayRate = Double.NaN;
        } else {
            int ref = Math.max(1, usable / 2);
            decayRate = levelMass[ref + 1] / levelMass[ref];
        }

        double meanQueueLength = 0.0;
        for (int i = 0; i < nLevels; i++) {
            meanQueueLength += ((double) i) * levelMass[i];
        }

        return new QsysBmapM1Result(theta, lambda, rho, q, A0, A1, B0, Bk, A, alpha, G, drift,
                decayRate, levelProb, levelMass[0], meanQueueLength, rho, lambda,
                nLevels - 1, truncError, "LINE:qsys_bmapm1");
    }

    /** Relative contribution the truncated tail would add to the mean level. */
    private static double levelTailError(Matrix levelProb, int levelMax) {
        int n = levelProb.getNumRows();
        int V = levelProb.getNumCols();
        double meanLevel = 0.0;
        double topMass = 0.0;
        for (int i = 0; i < n; i++) {
            double m = 0.0;
            for (int j = 0; j < V; j++) {
                m += levelProb.get(i, j);
            }
            meanLevel += ((double) i) * m;
            if (i == n - 1) {
                topMass = m;
            }
        }
        return levelMax * topMass / Math.max(meanLevel, Double.MIN_NORMAL);
    }

    /**
     * Level-truncated CTMC generator of the BMAP/M/1 queue and its stationary
     * distribution. Level n holds n customers in the system; the phase is the BMAP
     * state. Service fires only above level 0.
     */
    private static Matrix solveLevels(Matrix[] D, double mu, int V, int K, int levelMax) {
        int totalDim = (levelMax + 1) * V;

        ArrayList<int[]> idx = new ArrayList<int[]>();
        ArrayList<Double> vals = new ArrayList<Double>();

        // Local blocks: D0 on every level.
        for (int lvl = 0; lvl <= levelMax; lvl++) {
            addBlock(idx, vals, D[0], lvl * V, lvl * V, 1.0);
        }
        // Service: level n -> n-1 for n >= 1.
        for (int lvl = 1; lvl <= levelMax; lvl++) {
            for (int i = 0; i < V; i++) {
                idx.add(new int[] { lvl * V + i, (lvl - 1) * V + i });
                vals.add(Double.valueOf(mu));
            }
        }
        // Batch arrivals: level n -> n+k.
        for (int k = 1; k <= K; k++) {
            for (int lvl = 0; lvl + k <= levelMax; lvl++) {
                addBlock(idx, vals, D[k], lvl * V, (lvl + k) * V, 1.0);
            }
        }

        double[] rowSum = new double[totalDim];
        for (int e = 0; e < idx.size(); e++) {
            rowSum[idx.get(e)[0]] += vals.get(e).doubleValue();
        }

        DMatrixSparseTriplet triplet = new DMatrixSparseTriplet(totalDim, totalDim, idx.size() + 2 * totalDim);
        int lastCol = totalDim - 1;
        boolean[] diagWritten = new boolean[totalDim];
        for (int e = 0; e < idx.size(); e++) {
            int r = idx.get(e)[0];
            int c = idx.get(e)[1];
            double v = vals.get(e).doubleValue();
            if (r == c) {
                v -= rowSum[r];
                diagWritten[r] = true;
            }
            if (v != 0.0 && c != lastCol) {
                triplet.addItem(r, c, v);
            }
        }
        for (int r = 0; r < totalDim; r++) {
            if (!diagWritten[r] && rowSum[r] != 0.0 && r != lastCol) {
                triplet.addItem(r, r, -rowSum[r]);
            }
        }
        // Replace the last column with ones for normalization.
        for (int r = 0; r < totalDim; r++) {
            triplet.addItem(r, lastCol, 1.0);
        }

        DMatrixSparseCSC Qcsc = DConvertMatrixStruct.convert(triplet, (DMatrixSparseCSC) null);
        Matrix Q = new Matrix((org.ejml.data.DMatrix) Qcsc);

        Matrix Qt = Q.transpose();
        Matrix b = new Matrix(totalDim, 1);
        b.zero();
        b.set(totalDim - 1, 0, 1.0);
        Matrix piCol = new Matrix(totalDim, 1);
        Matrix.solveSafe(Qt, b, piCol);

        Matrix levelProb = new Matrix(levelMax + 1, V);
        double total = 0.0;
        for (int lvl = 0; lvl <= levelMax; lvl++) {
            for (int j = 0; j < V; j++) {
                double v = piCol.get(lvl * V + j, 0);
                if (v < 0) {
                    v = 0;
                }
                levelProb.set(lvl, j, v);
                total += v;
            }
        }
        if (total > 0) {
            for (int lvl = 0; lvl <= levelMax; lvl++) {
                for (int j = 0; j < V; j++) {
                    levelProb.set(lvl, j, levelProb.get(lvl, j) / total);
                }
            }
        }
        return levelProb;
    }

    private static void addBlock(ArrayList<int[]> idx, ArrayList<Double> vals, Matrix A,
                                 int rowBase, int colBase, double scale) {
        Iterator<MatrixEntry> it = A.nonZeroIterator();
        while (it.hasNext()) {
            MatrixEntry e = it.next();
            if (e.value == 0.0) {
                continue;
            }
            idx.add(new int[] { rowBase + e.row, colBase + e.col });
            vals.add(Double.valueOf(scale * e.value));
        }
    }

    private static Matrix asRow(Matrix v) {
        if (v.getNumRows() == 1) {
            return v;
        }
        Matrix out = new Matrix(1, v.length());
        for (int i = 0; i < v.length(); i++) {
            out.set(0, i, v.get(i));
        }
        return out;
    }

    private static double minEntry(Matrix A) {
        double m = Double.POSITIVE_INFINITY;
        for (int i = 0; i < A.getNumRows(); i++) {
            for (int j = 0; j < A.getNumCols(); j++) {
                m = Math.min(m, A.get(i, j));
            }
        }
        return m;
    }

    private static double maxAbsDiff(Matrix A, Matrix B) {
        double m = 0.0;
        for (int i = 0; i < A.getNumRows(); i++) {
            for (int j = 0; j < A.getNumCols(); j++) {
                m = Math.max(m, Math.abs(A.get(i, j) - B.get(i, j)));
            }
        }
        return m;
    }

    /** Computes v * A * e for a row vector v and a square matrix A. */
    private static double rowWeightedSum(Matrix v, Matrix A) {
        double s = 0.0;
        for (int i = 0; i < A.getNumRows(); i++) {
            for (int j = 0; j < A.getNumCols(); j++) {
                s += v.get(0, i) * A.get(i, j);
            }
        }
        return s;
    }
}
