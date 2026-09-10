/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.handlers;

import java.util.Arrays;
import java.util.List;

import jline.api.pfqn.lcfs.Pfqn_lcfsqn_mva;
import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.util.matrix.Matrix;

public final class Solver_mva_lcfsqn {
    private Solver_mva_lcfsqn() {}

    /**
     * Specialized MVA solver for LCFS + LCFS-PR 2-station networks.
     */
    public static MVAResult solver_mva_lcfsqn(NetworkStruct sn,
                                              SolverOptions options,
                                              int lcfsStat,
                                              int lcfsprStat) {
        long startTime = System.nanoTime();
        int M = sn.nstations;
        int nclasses = sn.nclasses;
        Matrix njobs = sn.njobs;

        Matrix alpha = new Matrix(1, nclasses);
        Matrix beta = new Matrix(1, nclasses);

        Matrix rates = sn.rates;
        for (int r = 0; r < nclasses; r++) {
            if (njobs.get(r) > 0) {
                double mu_lcfs = rates.get(lcfsStat, r);
                double mu_lcfspr = rates.get(lcfsprStat, r);

                if (mu_lcfs <= 0 || !Double.isFinite(mu_lcfs)) {
                    throw new RuntimeException("Invalid service rate at LCFS station for class " + r + ".");
                }
                if (mu_lcfspr <= 0 || !Double.isFinite(mu_lcfspr)) {
                    throw new RuntimeException("Invalid service rate at LCFS-PR station for class " + r + ".");
                }

                alpha.set(r, 1.0 / mu_lcfs);
                beta.set(r, 1.0 / mu_lcfspr);
            }
        }

        Matrix N = njobs;

        jline.api.pfqn.lcfs.LcfsqnMvaResult mvaResult = Pfqn_lcfsqn_mva.pfqn_lcfsqn_mva(alpha, beta, N);
        Matrix T_lcfs = mvaResult.getT();
        Matrix Q_lcfs = mvaResult.getQ();
        Matrix U_lcfs = mvaResult.getU();

        Matrix Q = new Matrix(M, nclasses);
        Matrix U = new Matrix(M, nclasses);
        Matrix T = new Matrix(M, nclasses);
        Matrix R = new Matrix(M, nclasses);
        Matrix X = new Matrix(1, nclasses);
        Matrix C = new Matrix(1, nclasses);

        for (int r = 0; r < nclasses; r++) {
            Q.set(lcfsStat, r, Q_lcfs.get(0, r));
            Q.set(lcfsprStat, r, Q_lcfs.get(1, r));
        }

        for (int r = 0; r < nclasses; r++) {
            U.set(lcfsStat, r, U_lcfs.get(0, r));
            U.set(lcfsprStat, r, U_lcfs.get(1, r));
        }

        for (int r = 0; r < nclasses; r++) {
            if (njobs.get(r) > 0) {
                X.set(r, T_lcfs.get(0, r));
                T.set(lcfsStat, r, T_lcfs.get(0, r));
                T.set(lcfsprStat, r, T_lcfs.get(0, r));
            }
        }

        List<Integer> stationsForR = Arrays.asList(lcfsStat, lcfsprStat);
        for (Integer k : stationsForR) {
            for (int r = 0; r < nclasses; r++) {
                if (T.get(k, r) > 0) {
                    R.set(k, r, Q.get(k, r) / T.get(k, r));
                }
            }
        }

        for (int r = 0; r < nclasses; r++) {
            if (njobs.get(r) > 0) {
                C.set(r, R.get(lcfsStat, r) + R.get(lcfsprStat, r));
            }
        }

        double runtime = (System.nanoTime() - startTime) / 1_000_000_000.0;

        MVAResult result = new MVAResult();
        result.QN = Q;
        result.UN = U;
        result.RN = R;
        result.TN = T;
        result.CN = C;
        result.XN = X;
        result.logNormConstAggr = Double.NaN;
        result.runtime = runtime;
        // suppress unused warning
        if (options == null) { /* no-op */ }

        return result;
    }
}
