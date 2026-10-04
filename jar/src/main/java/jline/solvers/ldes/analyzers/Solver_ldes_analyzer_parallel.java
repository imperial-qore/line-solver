/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ldes.analyzers;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.distribution.TDistribution;

import jline.lang.NetworkStruct;
import jline.solvers.ldes.LDESOptions;
import jline.solvers.ldes.LDESResult;
import jline.solvers.ldes.SolverLDES;
import jline.solvers.ldes.handlers.Solver_ssj;
import jline.util.matrix.Matrix;

/**
 * Parallel replication analyzer for LDES solver.
 *
 * Runs multiple independent simulation replications in parallel and aggregates
 * results using cross-replication statistics for confidence intervals.
 */
public final class Solver_ldes_analyzer_parallel {
    private Solver_ldes_analyzer_parallel() {}

    public static LDESResult solver_ldes_analyzer_parallel(
            NetworkStruct sn,
            LDESOptions options,
            SolverLDES solverLDES) {
        long Tstart = System.nanoTime();
        final int numReplications = options.replications;

        final int M = sn.nstations;
        final int K = sn.nclasses;

        // see _kb/09-ldes-and-cache.md (Ensemble transient section)
        final boolean isTransient = options.timespan != null
                && options.timespan.length >= 2
                && Double.isFinite(options.timespan[1]);
        final ConcurrentHashMap<Integer, LDESResult> tranResults =
                new ConcurrentHashMap<Integer, LDESResult>();

        // Thread-safe storage for per-replication results
        final ConcurrentHashMap<Integer, Matrix> QNs = new ConcurrentHashMap<Integer, Matrix>();
        final ConcurrentHashMap<Integer, Matrix> UNs = new ConcurrentHashMap<Integer, Matrix>();
        final ConcurrentHashMap<Integer, Matrix> RNs = new ConcurrentHashMap<Integer, Matrix>();
        final ConcurrentHashMap<Integer, Matrix> TNs = new ConcurrentHashMap<Integer, Matrix>();
        final ConcurrentHashMap<Integer, Matrix> ANs = new ConcurrentHashMap<Integer, Matrix>();
        final ConcurrentHashMap<Integer, Matrix> CNs = new ConcurrentHashMap<Integer, Matrix>();
        final ConcurrentHashMap<Integer, Matrix> XNs = new ConcurrentHashMap<Integer, Matrix>();

        // Impatience metrics storage
        final ConcurrentHashMap<Integer, Matrix> renegedCounts = new ConcurrentHashMap<Integer, Matrix>();
        final ConcurrentHashMap<Integer, Matrix> avgRenegingWaits = new ConcurrentHashMap<Integer, Matrix>();
        final ConcurrentHashMap<Integer, Matrix> renegingRates = new ConcurrentHashMap<Integer, Matrix>();
        final ConcurrentHashMap<Integer, Matrix> balkedCounts = new ConcurrentHashMap<Integer, Matrix>();
        final ConcurrentHashMap<Integer, Matrix> balkingProbs = new ConcurrentHashMap<Integer, Matrix>();
        final ConcurrentHashMap<Integer, Matrix> retriedCounts = new ConcurrentHashMap<Integer, Matrix>();
        final ConcurrentHashMap<Integer, Matrix> retrialDroppedCounts = new ConcurrentHashMap<Integer, Matrix>();
        final ConcurrentHashMap<Integer, Matrix> avgOrbitSizes = new ConcurrentHashMap<Integer, Matrix>();
        // Fork-Join quorum sibling-drop rate (station-indexed)
        final ConcurrentHashMap<Integer, Matrix> dropRateJoins = new ConcurrentHashMap<Integer, Matrix>();

        // Divide samples across replications
        final int samplesPerReplication = Math.max(1000, (options.samples + numReplications - 1) / numReplications);

        final NetworkStruct snFinal = sn;
        final LDESOptions optionsFinal = options;
        final SolverLDES solverLDESFinal = solverLDES;

        final List<Future<?>> futures = new ArrayList<Future<?>>();
        for (int repIdx = 0; repIdx < numReplications; repIdx++) {
            final int replicationIndex = repIdx;
            futures.add(solverLDES.threadPool.submit(new Runnable() {
                @Override
                public void run() {
                    // Reproducibility comes from repOptions.seed below, which is
                    // derived from the master seed and the replication index.

                    // Create per-replication options copy
                    LDESOptions repOptions = optionsFinal.copy();
                    repOptions.samples = samplesPerReplication;
                    repOptions.seed = optionsFinal.seed + replicationIndex;

                    // see _kb/09-ldes-and-cache.md (Ensemble transient section)
                    LDESResult repResult;
                    if (isTransient) {
                        repResult = Solver_ssj.solver_ssj_transient(snFinal, repOptions, null);
                        tranResults.put(replicationIndex, repResult);
                    } else {
                        repResult = Solver_ssj.solver_ssj(snFinal, repOptions, null);
                    }

                    QNs.put(replicationIndex, repResult.QN);
                    UNs.put(replicationIndex, repResult.UN);
                    RNs.put(replicationIndex, repResult.RN);
                    TNs.put(replicationIndex, repResult.TN);
                    if (repResult.AN != null) {
                        ANs.put(replicationIndex, repResult.AN);
                    } else {
                        ANs.put(replicationIndex, new Matrix(M, K));
                    }
                    CNs.put(replicationIndex, repResult.CN);
                    XNs.put(replicationIndex, repResult.XN);

                    if (repResult.renegedCustomers != null) {
                        renegedCounts.put(replicationIndex, repResult.renegedCustomers);
                    }
                    if (repResult.avgRenegingWaitTime != null) {
                        avgRenegingWaits.put(replicationIndex, repResult.avgRenegingWaitTime);
                    }
                    if (repResult.renegingRate != null) {
                        renegingRates.put(replicationIndex, repResult.renegingRate);
                    }
                    if (repResult.balkedCustomers != null) {
                        balkedCounts.put(replicationIndex, repResult.balkedCustomers);
                    }
                    if (repResult.balkingProbability != null) {
                        balkingProbs.put(replicationIndex, repResult.balkingProbability);
                    }
                    if (repResult.retriedCustomers != null) {
                        retriedCounts.put(replicationIndex, repResult.retriedCustomers);
                    }
                    if (repResult.retrialDropped != null) {
                        retrialDroppedCounts.put(replicationIndex, repResult.retrialDropped);
                    }
                    if (repResult.avgOrbitSize != null) {
                        avgOrbitSizes.put(replicationIndex, repResult.avgOrbitSize);
                    }
                    if (repResult.DropRateJoin != null) {
                        dropRateJoins.put(replicationIndex, repResult.DropRateJoin);
                    }
                }
            }));
        }

        solverLDES.threadPool.shutdown();
        try {
            if (!solverLDES.threadPool.awaitTermination(600, TimeUnit.SECONDS)) {
                solverLDES.threadPool.shutdownNow();
                throw new RuntimeException(
                        "Parallel LDES execution timed out after 600s with "
                                + QNs.size() + " of " + numReplications + " replications complete");
            }
        } catch (InterruptedException e) {
            solverLDES.threadPool.shutdownNow();
            Thread.currentThread().interrupt();
            throw new RuntimeException("Parallel LDES execution interrupted", e);
        }

        // see _kb/09-ldes-and-cache.md (Ensemble transient section (replication failure surfacing))
        for (Future<?> f : futures) {
            try {
                f.get();
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                throw new RuntimeException("Parallel LDES execution interrupted", e);
            } catch (ExecutionException e) {
                throw new RuntimeException("LDES replication failed: " + e.getCause(), e.getCause());
            }
        }
        if (QNs.size() != numReplications) {
            throw new RuntimeException("LDES produced " + QNs.size() + " of "
                    + numReplications + " replications");
        }

        Matrix QN = averageMatrices(QNs.values(), numReplications);
        Matrix UN = averageMatrices(UNs.values(), numReplications);
        Matrix RN = averageMatrices(RNs.values(), numReplications);
        Matrix TN = averageMatrices(TNs.values(), numReplications);
        Matrix AN = averageMatrices(ANs.values(), numReplications);
        Matrix CN = averageMatrices(CNs.values(), numReplications);
        Matrix XN = averageMatrices(XNs.values(), numReplications);

        double alpha = 1.0 - options.confint;
        double tCrit;
        if (numReplications > 1) {
            TDistribution tDist = new TDistribution((double) (numReplications - 1));
            tCrit = tDist.inverseCumulativeProbability(1.0 - alpha / 2.0);
        } else {
            tCrit = 0.0;
        }

        Matrix QNCI = computeCrossReplicationCI(new ArrayList<Matrix>(QNs.values()), QN, tCrit, numReplications);
        Matrix UNCI = computeCrossReplicationCI(new ArrayList<Matrix>(UNs.values()), UN, tCrit, numReplications);
        Matrix RNCI = computeCrossReplicationCI(new ArrayList<Matrix>(RNs.values()), RN, tCrit, numReplications);
        Matrix TNCI = computeCrossReplicationCI(new ArrayList<Matrix>(TNs.values()), TN, tCrit, numReplications);
        Matrix ANCI = computeCrossReplicationCI(new ArrayList<Matrix>(ANs.values()), AN, tCrit, numReplications);
        Matrix WNCI = RNCI.copy();

        Matrix QNRelPrec = computeRelativePrecision(QNCI, QN);
        Matrix UNRelPrec = computeRelativePrecision(UNCI, UN);
        Matrix RNRelPrec = computeRelativePrecision(RNCI, RN);
        Matrix TNRelPrec = computeRelativePrecision(TNCI, TN);

        LDESResult result = new LDESResult();
        result.QN = QN;
        result.UN = UN;
        result.RN = RN;
        result.TN = TN;
        result.AN = AN;
        result.CN = CN;
        result.XN = XN;
        result.sn = sn;
        result.method = "parallel";
        result.runtime = ((double) (System.nanoTime() - Tstart)) / 1_000_000_000.0;

        result.QNCI = QNCI;
        result.UNCI = UNCI;
        result.RNCI = RNCI;
        result.TNCI = TNCI;
        result.ANCI = ANCI;
        result.WNCI = WNCI;
        result.QNRelPrec = QNRelPrec;
        result.UNRelPrec = UNRelPrec;
        result.RNRelPrec = RNRelPrec;
        result.TNRelPrec = TNRelPrec;

        result.converged = checkConvergence(QNRelPrec, UNRelPrec, RNRelPrec, TNRelPrec, options.cnvgtol);
        result.stoppingReason = result.converged ? "convergence" : "max_replications";
        result.convergenceBatches = numReplications;

        // Transient ensemble: average the per-bucket series across replications.
        if (isTransient && tranResults.size() == numReplications) {
            List<LDESResult> ordered = new ArrayList<LDESResult>();
            for (int r = 0; r < numReplications; r++) {
                ordered.add(tranResults.get(r));
            }
            LDESResult avg = Solver_ssj.averageTransientResults(ordered);
            result.t = avg.t;
            result.QNt = avg.QNt;
            result.UNt = avg.UNt;
            result.TNt = avg.TNt;
        }

        if (!renegedCounts.isEmpty()) {
            result.renegedCustomers = averageMatrices(renegedCounts.values(), numReplications);
        }
        if (!avgRenegingWaits.isEmpty()) {
            result.avgRenegingWaitTime = averageMatrices(avgRenegingWaits.values(), numReplications);
        }
        if (!renegingRates.isEmpty()) {
            result.renegingRate = averageMatrices(renegingRates.values(), numReplications);
        }
        if (!balkedCounts.isEmpty()) {
            result.balkedCustomers = averageMatrices(balkedCounts.values(), numReplications);
        }
        if (!balkingProbs.isEmpty()) {
            result.balkingProbability = averageMatrices(balkingProbs.values(), numReplications);
        }
        if (!retriedCounts.isEmpty()) {
            result.retriedCustomers = averageMatrices(retriedCounts.values(), numReplications);
        }
        if (!retrialDroppedCounts.isEmpty()) {
            result.retrialDropped = averageMatrices(retrialDroppedCounts.values(), numReplications);
        }
        if (!avgOrbitSizes.isEmpty()) {
            result.avgOrbitSize = averageMatrices(avgOrbitSizes.values(), numReplications);
        }
        if (!dropRateJoins.isEmpty()) {
            result.DropRateJoin = averageMatrices(dropRateJoins.values(), numReplications);
        }

        return result;
    }

    /**
     * Averages the per-replication matrices. The divisor is the number of matrices
     * actually supplied, not the number of replications requested: dividing by the
     * latter would scale the result by (supplied / requested) whenever a replication
     * failed to store its result, which is a silent error rather than a loud one.
     */
    private static Matrix averageMatrices(Collection<Matrix> matrices, int count) {
        if (matrices.isEmpty()) {
            throw new IllegalArgumentException("Cannot average empty collection");
        }
        Matrix first = matrices.iterator().next();
        Matrix result = new Matrix(first.getNumRows(), first.getNumCols());
        result.zero();
        for (Matrix m : matrices) {
            result.addEq(m);
        }
        return Matrix.scaleMult(result, 1.0 / (double) matrices.size());
    }

    private static Matrix computeCrossReplicationCI(
            List<Matrix> replicationResults,
            Matrix grandMean,
            double tCrit,
            int numReplications) {
        int rows = grandMean.getNumRows();
        int cols = grandMean.getNumCols();
        Matrix ciHalfWidth = new Matrix(rows, cols);

        if (numReplications < 2) {
            ciHalfWidth.fill(0.0);
            return ciHalfWidth;
        }

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                double mean = grandMean.get(i, j);
                double sumSquaredDiff = 0.0;
                for (Matrix rep : replicationResults) {
                    double diff = rep.get(i, j) - mean;
                    sumSquaredDiff += diff * diff;
                }
                double variance = sumSquaredDiff / (numReplications - 1);
                double stdError = Math.sqrt(variance / numReplications);
                ciHalfWidth.set(i, j, tCrit * stdError);
            }
        }

        return ciHalfWidth;
    }

    private static Matrix computeRelativePrecision(Matrix ciHalfWidth, Matrix mean) {
        Matrix result = new Matrix(mean.getNumRows(), mean.getNumCols());
        for (int i = 0; i < mean.getNumRows(); i++) {
            for (int j = 0; j < mean.getNumCols(); j++) {
                double m = mean.get(i, j);
                double ci = ciHalfWidth.get(i, j);
                if (m > 0) {
                    result.set(i, j, ci / m);
                } else {
                    result.set(i, j, 0.0);
                }
            }
        }
        return result;
    }

    private static boolean checkConvergence(
            Matrix qnRelPrec,
            Matrix unRelPrec,
            Matrix rnRelPrec,
            Matrix tnRelPrec,
            double tolerance) {
        Matrix[] allMatrices = new Matrix[] { qnRelPrec, unRelPrec, rnRelPrec, tnRelPrec };
        for (Matrix m : allMatrices) {
            for (int i = 0; i < m.getNumRows(); i++) {
                for (int j = 0; j < m.getNumCols(); j++) {
                    double value = m.get(i, j);
                    if (value > 0 && !Double.isInfinite(value) && !Double.isNaN(value) && value > tolerance) {
                        return false;
                    }
                }
            }
        }
        return true;
    }
}
