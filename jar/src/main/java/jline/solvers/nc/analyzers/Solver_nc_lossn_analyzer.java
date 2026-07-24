/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.nc.analyzers;

import jline.api.lossn.Lossn_erlangfp;
import jline.api.lossn.Lossn_mci;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.nc.NCResult;
import jline.util.matrix.Matrix;

public final class Solver_nc_lossn_analyzer {
    private Solver_nc_lossn_analyzer() {}

    /**
     * Analyzes open loss networks with FCR.
     *
     * Method selection (options.method):
     *   "erlangfp" (default) - Erlang fixed-point (reduced-load) approximation.
     *   "mci"                - Monte Carlo importance-sampling summation
     *                          (Ross-Wang 1992): estimates the normalization
     *                          constant g(C) and class blocking with confidence
     *                          intervals (options.samples/options.seed).
     */
    public static NCResult solver_nc_lossn_analyzer(NetworkStruct sn, SolverOptions options) {
        long Tstart = System.nanoTime();
        int K = sn.nclasses;
        int M = sn.nstations;

        double[] nu = new double[K];
        for (int r = 0; r < K; r++) {
            int sourceIdx = (int) sn.refstat.get(r);
            nu[r] = sn.rates.get(sourceIdx, r);
        }

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
            if (regionMatrix.get(i, K) >= 0) {
                hasConstraint = true;
            }
            if (hasConstraint) {
                delayIdx = i;
                break;
            }
        }

        double globalMax = regionMatrix.get(delayIdx, K);
        double[] classMax = new double[K];
        for (int r = 0; r < K; r++) {
            classMax[r] = regionMatrix.get(delayIdx, r);
        }

        double globalMaxVal = (globalMax < 0) ? 1e6 : globalMax;
        for (int r = 0; r < K; r++) {
            if (classMax[r] < 0) {
                classMax[r] = 1e6;
            }
        }

        int J = K + 1;
        Matrix A = new Matrix(J, K);
        A.fill(0.0);
        for (int r = 0; r < K; r++) {
            A.set(0, r, 1.0);
        }
        for (int r = 0; r < K; r++) {
            A.set(r + 1, r, 1.0);
        }

        Matrix C_vec = new Matrix(J, 1);
        C_vec.set(0, 0, globalMaxVal);
        for (int r = 0; r < K; r++) {
            C_vec.set(r + 1, 0, classMax[r]);
        }

        // Method selection
        boolean useMCI = false;
        if (options != null && options.method != null) {
            String[] tokens = options.method.toLowerCase().split("[./]");
            for (int t = 0; t < tokens.length; t++) {
                if (tokens[t].equals("mci")) {
                    useMCI = true;
                }
            }
        }

        Matrix QLen;
        int niter;
        double lG = Double.NaN;   // normalization constant (finite only for mci)
        String actualMethod;
        if (useMCI) {
            int nsamples = (options != null && options.samples > 0) ? options.samples : 100000;
            long seed = (options != null) ? (long) options.seed : -1L;
            Ret.lossnMCI mciResult = Lossn_mci.lossn_mci(new Matrix(nu), A, C_vec,
                    nsamples, null, seed, 0.05);
            QLen = mciResult.qLen;
            lG = mciResult.lG;
            niter = mciResult.nsamples;
            actualMethod = "mci";
        } else {
            Ret.lossnErlangFP lossnResult = Lossn_erlangfp.lossn_erlangfp(new Matrix(nu), A, C_vec);
            QLen = lossnResult.qLen;
            niter = lossnResult.niter;
            actualMethod = "erlangfp";
        }

        Matrix Q = new Matrix(M, K);
        Q.fill(0.0);
        Matrix U = new Matrix(M, K);
        U.fill(0.0);
        Matrix T = new Matrix(M, K);
        T.fill(0.0);
        Matrix R = new Matrix(M, K);
        R.fill(0.0);

        for (int r = 0; r < K; r++) {
            double effTput = QLen.get(r);
            T.set(delayIdx, r, effTput);
            // Source emits the accepted (post-drop) rate into the network so that
            // the routing-based arrival-rate computation (snGetArvRFromTput) yields
            // a non-zero ArvR at the delay, consistent with the flow-conserving
            // departure throughput and with the simulated rate from SolverJMT.
            int sourceIdx = (int) sn.refstat.get(r);
            T.set(sourceIdx, r, effTput);
            double mu_r = sn.rates.get(delayIdx, r);
            Q.set(delayIdx, r, effTput / mu_r);
            R.set(delayIdx, r, 1.0 / mu_r);
            U.set(delayIdx, r, effTput / mu_r);
        }

        Matrix X = new Matrix(1, K);
        for (int r = 0; r < K; r++) {
            X.set(0, r, QLen.get(r));
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
}
