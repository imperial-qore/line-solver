/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.nc.analyzers;

import java.util.ArrayList;
import java.util.List;

import jline.api.lossn.Lossn_erlangfp;
import jline.api.lossn.Lossn_mci;
import jline.api.lossn.Lossn_manjunath;
import jline.api.lossn.Lossn_rec;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.nc.NCResult;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Solver_nc_lossn_analyzer {
    private Solver_nc_lossn_analyzer() {}

    /**
     * Analyzes open loss networks with FCR.
     *
     * Handles open queueing networks with a single multiclass Delay node inside a
     * Finite Capacity Region (FCR) with DROP policy.
     *
     * The FCR admission rule is A n &lt;= C on the per-class occupancy vector n of
     * the region, where the rows of A are assembled from every constraint the
     * region declares: the global job cap, the memory budget weighted by the
     * per-class sizes, the per-class job caps, and any explicit linear constraint
     * set with FiniteCapacityRegion.setConstraint. Rows left unbounded are
     * DROPPED rather than given a surrogate capacity.
     *
     * Method selection (options.method):
     *   "rec"             - MDD-rec (Lossn_rec): the same constant as the exact
     *                       sum over the admissible set, obtained by one memoised
     *                       walk of the decision diagram holding it. Places no
     *                       integrality demand on A or C.
     *   "exact" (default) - Manjunath-Sikdar transform (Lossn_manjunath): the
     *                       normalization constant is obtained exactly as a
     *                       multidimensional contour integral evaluated by
     *                       residues. Requires integer A and C.
     *   "erlangfp"        - Erlang fixed-point (reduced-load) approximation.
     *   "mci"             - Monte Carlo importance-sampling summation
     *                       (Ross-Wang 1992): estimates the normalization
     *                       constant g(C) and class blocking with confidence
     *                       intervals (options.samples/options.seed).
     *
     * The default is the residue transform on an integral region and MDD-rec on a
     * fractional one. It used to fall back to "erlangfp" there, an approximation,
     * because the residue argument counts whole units; MDD-rec needs only that
     * the admissible set be finite and bounded per coordinate, which it still is,
     * so the fractional case is now exact as well. Historically the default fell
     * back to "erlangfp" when the region declares fractional
     * class sizes or capacities, since the residue argument counts whole units.
     */
    public static NCResult solver_nc_lossn_analyzer(NetworkStruct sn, SolverOptions options) {
        long Tstart = System.nanoTime();
        int K = sn.nclasses;
        int M = sn.nstations;

        // 1. Locate the delay station inside the region
        Matrix regionMatrix = sn.region.get(0);
        int delayIdx = -1;
        for (int i = 0; i < M; i++) {
            boolean hasConstraint = false;
            for (int r = 0; r < K; r++) {
                if (regionMatrix.get(i, r) >= 0) {
                    hasConstraint = true;
                    break;
                }
            }
            if (regionMatrix.getNumCols() > K && regionMatrix.get(i, K) >= 0) {
                hasConstraint = true;
            }
            if (hasConstraint) {
                delayIdx = i;
                break;
            }
        }
        if (delayIdx < 0) {
            throw new RuntimeException("solver_nc_lossn_analyzer: no station found inside the "
                    + "finite capacity region");
        }

        // 2. Offered load per class. A route carries nu_r = arrival rate times mean
        // holding time INSIDE the region, i.e. the visit ratio at the delay divided
        // by its service rate. The bare arrival rate would be right only for unit
        // mean service times, and reports the wrong blocking otherwise.
        double[] nu = new double[K];
        double[] lambda = new double[K];
        double[] mu = new double[K];
        for (int r = 0; r < K; r++) {
            int sourceIdx = (int) sn.refstat.get(r);
            lambda[r] = sn.rates.get(sourceIdx, r);
            double vRatio = 1.0;
            if (sn.chains != null && sn.visits != null) {
                for (int c = 0; c < sn.chains.getNumRows(); c++) {
                    if (sn.chains.get(c, r) == 0) {
                        continue;
                    }
                    Matrix vc = sn.visits.get(c);
                    if (vc != null) {
                        double vref = vc.get(sourceIdx, r);
                        if (vref > 0) {
                            vRatio = vc.get(delayIdx, r) / vref;
                        }
                    }
                    break;
                }
            }
            mu[r] = sn.rates.get(delayIdx, r);
            if (mu[r] > 0) {
                nu[r] = lambda[r] * vRatio / mu[r];
            }
        }

        // 3. Assemble the admission constraints A n <= C_vec of the region
        List<double[]> rows = new ArrayList<double[]>();
        List<Double> rhs = new ArrayList<Double>();
        lossnRegionConstraints(sn, 0, delayIdx, K, rows, rhs);
        int J = rows.size();
        Matrix A = new Matrix(J, K);
        A.fill(0.0);
        Matrix C_vec = new Matrix(J, 1);
        for (int j = 0; j < J; j++) {
            for (int r = 0; r < K; r++) {
                A.set(j, r, rows.get(j)[r]);
            }
            C_vec.set(j, 0, rhs.get(j).doubleValue());
        }

        // 4. Method selection, following the reference's precedence: an explicit
        // "mci" wins, then an explicit "erlangfp", then the exact transform. Only
        // "default" falls back, and only on a fractional region: answering an
        // explicit "exact" with the Erlang approximation would report an
        // approximation under the name of an exact method.
        boolean hasMci = false;
        boolean hasErlangfp = false;
        boolean hasRec = false;
        boolean hasExact = false;
        if (options != null && options.method != null) {
            String[] tokens = options.method.toLowerCase().split("[./]");
            for (int t = 0; t < tokens.length; t++) {
                if (tokens[t].equals("mci")) {
                    hasMci = true;
                } else if (tokens[t].equals("erlangfp")) {
                    hasErlangfp = true;
                } else if (tokens[t].equals("rec")) {
                    hasRec = true;
                } else if (tokens[t].equals("exact") || tokens[t].equals("manjunath") || tokens[t].equals("ms")) {
                    hasExact = true;
                }
            }
        }
        boolean isIntegral = true;
        for (int j = 0; j < J; j++) {
            if (Math.abs(C_vec.get(j, 0) - Math.rint(C_vec.get(j, 0))) > 1e-9) {
                isIntegral = false;
            }
            for (int r = 0; r < K; r++) {
                if (Math.abs(A.get(j, r) - Math.rint(A.get(j, r))) > 1e-9) {
                    isIntegral = false;
                }
            }
        }

        String chosen;
        if (hasMci) {
            chosen = "mci";
        } else if (hasErlangfp) {
            chosen = "erlangfp";
        } else if (hasRec) {
            chosen = "rec";
        } else if (hasExact) {
            chosen = "exact";
        } else {
            // the residue argument counts whole units; MDD-rec does not
            chosen = isIntegral ? "exact" : "rec";
        }

        Matrix QLen;
        Matrix Loss;
        int niter;
        double lG = Double.NaN;   // normalization constant (finite for exact and mci)
        String actualMethod;
        if (chosen.equals("mci")) {
            int nsamples = (options != null && options.samples > 0) ? options.samples : 100000;
            long seed = (options != null) ? (long) options.seed : -1L;
            Ret.lossnMCI mciResult = Lossn_mci.lossn_mci(new Matrix(nu), A, C_vec,
                    nsamples, null, seed, 0.05);
            QLen = mciResult.qLen;
            Loss = mciResult.lossProb;
            lG = mciResult.lG;
            // The sampler does not iterate; the reference reports the realised
            // sample count in the slot the other methods use for the iteration
            // count, so a caller can tell how much work produced the estimate.
            niter = mciResult.nsamples;
            actualMethod = "lossn.mci";
        } else if (chosen.equals("rec")) {
            Lossn_rec.LossnRecResult recResult = Lossn_rec.lossn_rec(new Matrix(nu), A, C_vec);
            QLen = new Matrix(recResult.QLen);
            Loss = new Matrix(recResult.Loss);
            lG = recResult.lG;
            niter = recResult.niter;
            actualMethod = "lossn.rec";
        } else if (chosen.equals("exact")) {
            Ret.lossnManjunath msResult = Lossn_manjunath.lossn_manjunath(new Matrix(nu), A, C_vec);
            QLen = msResult.qLen;
            Loss = msResult.lossProb;
            lG = msResult.lG;
            // The transform is direct, so the iteration count is 1 by definition.
            niter = msResult.niter;
            actualMethod = "lossn.exact";
        } else {
            Ret.lossnErlangFP lossnResult = Lossn_erlangfp.lossn_erlangfp(new Matrix(nu), A, C_vec);
            QLen = lossnResult.qLen;
            Loss = lossnResult.lossProb;
            niter = lossnResult.niter;
            actualMethod = "lossn.erlangfp";
        }

        // 5. Convert to standard outputs. QLen is the CARRIED load E[n_r], so the
        // carried throughput follows from Little's law at the infinite server.
        Matrix Q = new Matrix(M, K);
        Q.fill(0.0);
        Matrix U = new Matrix(M, K);
        U.fill(0.0);
        Matrix T = new Matrix(M, K);
        T.fill(0.0);
        Matrix R = new Matrix(M, K);
        R.fill(0.0);
        Matrix X = new Matrix(1, K);
        X.fill(0.0);

        for (int r = 0; r < K; r++) {
            double Xc = lambda[r] * (1.0 - Loss.get(r));   // carried (accepted) rate
            X.set(0, r, Xc);
            T.set(delayIdx, r, Xc);
            // The source emits the accepted (post-drop) rate into the network so
            // that the routing-based arrival-rate computation (snGetArvRFromTput)
            // yields a non-zero ArvR at the delay, consistent with the
            // flow-conserving departure throughput and with the simulated rate
            // from SolverJMT.
            int sourceIdx = (int) sn.refstat.get(r);
            T.set(sourceIdx, r, Xc);
            Q.set(delayIdx, r, QLen.get(r));               // mean number in the region
            // An infinite server never queues, so the response time is the service
            // time and the "utilization" is the mean number of busy servers.
            if (mu[r] > 0) {
                R.set(delayIdx, r, 1.0 / mu[r]);
            }
            U.set(delayIdx, r, QLen.get(r));
        }

        Matrix CN = new Matrix(1, K);
        CN.fill(0.0);

        NCResult res = new NCResult();
        res.QN = Q;
        res.UN = U;
        res.RN = R;
        res.TN = T;
        res.XN = X;
        res.CN = CN;
        res.lG = lG;
        res.it = niter;
        res.iter = niter;
        res.method = actualMethod;
        res.runtime = (double) (System.nanoTime() - Tstart) / 1000000000.0;

        return res;
    }

    /**
     * Rows of the admission rule A n &lt;= C for region f, in the order in which the
     * simulation engines test them.
     *
     * Every row is a function of the per-class occupancy of the region only, which
     * is what the FiniteCapacityRegion API can express, so no row can distinguish
     * stations inside the region. Rows the region leaves unbounded are DROPPED
     * rather than given a surrogate capacity: a 1e6 stand-in would silently turn
     * an unconstrained dimension into a truncation at 1e6 and make the exact
     * transform allocate a dimension of a million coefficients for a constraint
     * that does not exist.
     */
    private static void lossnRegionConstraints(NetworkStruct sn, int f, int stationIdx, int K,
                                               List<double[]> rows, List<Double> rhs) {
        Matrix regionMatrix = sn.region.get(f);

        // Global job cap: sum_r n_r <= globalMaxJobs
        if (regionMatrix.getNumCols() > K) {
            double globalMax = regionMatrix.get(stationIdx, K);
            if (globalMax >= 0) {
                double[] row = new double[K];
                for (int r = 0; r < K; r++) {
                    row[r] = 1.0;
                }
                rows.add(row);
                rhs.add(Double.valueOf(globalMax));
            }
        }

        // Memory budget: sum_r classSize_r n_r <= globalMaxMemory. The class sizes
        // are the row, so this is the one row that can legitimately be fractional.
        if (sn.regionmaxmem != null && sn.regionmaxmem.size() > f) {
            Matrix memMatrix = sn.regionmaxmem.get(f);
            if (memMatrix != null && stationIdx < memMatrix.getNumRows()) {
                double maxmem = memMatrix.get(stationIdx, 0);
                if (maxmem >= 0) {
                    double[] row = new double[K];
                    for (int r = 0; r < K; r++) {
                        row[r] = (sn.regionsz != null && f < sn.regionsz.getNumRows()
                                && r < sn.regionsz.getNumCols()) ? sn.regionsz.get(f, r) : 1.0;
                    }
                    rows.add(row);
                    rhs.add(Double.valueOf(maxmem));
                }
            }
        }

        // Per-class job caps: n_r <= classMaxJobs_r, already folded with
        // classMaxMemory when the region was added.
        for (int r = 0; r < K; r++) {
            double classMax = regionMatrix.get(stationIdx, r);
            if (classMax >= 0) {
                double[] row = new double[K];
                row[r] = 1.0;
                rows.add(row);
                rhs.add(Double.valueOf(classMax));
            }
        }

        // Explicit linear constraints from FiniteCapacityRegion.setConstraint
        if (sn.regionlincon != null && sn.regionlincon.containsKey(Integer.valueOf(f))) {
            MatrixCell lincon = sn.regionlincon.get(Integer.valueOf(f));
            if (lincon != null && lincon.size() >= 2) {
                Matrix linA = lincon.get(0);
                Matrix linB = lincon.get(1);
                if (linA != null && linB != null) {
                    for (int k = 0; k < linA.getNumRows(); k++) {
                        double[] row = new double[K];
                        for (int r = 0; r < K && r < linA.getNumCols(); r++) {
                            row[r] = linA.get(k, r);
                        }
                        rows.add(row);
                        rhs.add(Double.valueOf(k < linB.getNumRows() ? linB.get(k, 0) : 0.0));
                    }
                }
            }
        }

        if (rows.isEmpty()) {
            throw new RuntimeException("solver_nc_lossn_analyzer: the finite capacity region "
                    + "declares no bounded constraint, so it admits every arrival and is not a "
                    + "loss network; give it a global job cap, a memory budget, a per-class cap or "
                    + "an explicit linear constraint");
        }
    }
}
