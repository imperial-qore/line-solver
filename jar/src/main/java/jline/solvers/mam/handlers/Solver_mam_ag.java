/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mam.handlers;

import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.mam.MAMOptions;
import jline.solvers.mam.MAMResult;
import jline.util.matrix.Matrix;

import static jline.io.InputOutput.line_warning;

/**
 * RCAT-based solver for SolverMAM.
 */
public final class Solver_mam_ag {
    private Solver_mam_ag() {}

    private static final String MFILENAME = "solver_mam_ag";

    public static MAMResult solver_mam_ag(NetworkStruct sn, SolverOptions options) {
        int M = sn.nstations;
        int K = sn.nclasses;

        int maxStates = (options instanceof MAMOptions) ? ((MAMOptions) options).maxStates : 100;

        double tol = options.iter_tol > 0 ? options.iter_tol : 1e-6;
        int maxiter = options.iter_max > 0 ? options.iter_max : 1000;

        RCATModel rcat = Solver_mam_build_ag.solver_mam_build_ag(sn, maxStates);

        int numProcesses = rcat.N.length;
        int numActions = rcat.actionMap.size();

        if (numProcesses == 0 || (numActions == 0 && numProcesses > 1)) {
            line_warning(MFILENAME, "Network could not be mapped to RCAT format.");
            return createEmptyMAMResult(M, K);
        }

        String method = options.method;
        if ("default".equals(method)) {
            method = "inap";
        }

        INAPResult inapResult;
        if ("inap".equals(method)) {
            inapResult = Solver_mam_inap.solver_mam_inap(rcat, tol, maxiter, options.seed, "inap");
        } else if ("inapplus".equals(method)) {
            inapResult = Solver_mam_inap.solver_mam_inap(rcat, tol, maxiter, options.seed, "inapplus");
        } else if ("inapinf".equals(method)) {
            // Matrix-geometric INAP: solve isolated open components exactly on
            // the infinite state space (geometric tail), no truncation.
            boolean[] isOpenProc = new boolean[numProcesses];
            for (int p = 0; p < numProcesses; p++) {
                for (int ist = 0; ist < M && !isOpenProc[p]; ist++) {
                    for (int r = 0; r < K; r++) {
                        if ((int) rcat.processMap.get(ist, r) == p) {
                            isOpenProc[p] = Double.isInfinite(sn.njobs.get(r));
                            break;
                        }
                    }
                }
            }
            inapResult = Solver_mam_inap.solver_mam_inapinf(rcat, isOpenProc, tol, maxiter);
        } else if ("exact".equals(method)) {
            line_warning(MFILENAME, "AutoCAT (exact method) is not yet available in JAR. Falling back to INAP.");
            inapResult = Solver_mam_inap.solver_mam_inap(rcat, tol, maxiter, options.seed, "inap");
        } else {
            line_warning(MFILENAME, "Unknown method: " + method + ". Using INAP.");
            inapResult = Solver_mam_inap.solver_mam_inap(rcat, tol, maxiter, options.seed, "inap");
        }

        MetricsResult metrics = Solver_mam_metrics.solver_mam_metrics(sn, inapResult, rcat);

        MAMResult result = new MAMResult();
        result.QN = metrics.QN;
        result.UN = metrics.UN;
        result.RN = metrics.RN;
        result.TN = metrics.TN;
        result.CN = metrics.CN;
        result.XN = metrics.XN;
        result.iter = inapResult.iter;
        result.method = method;
        result.actionRates = inapResult.x;
        result.equilibrium = inapResult.pi;
        result.generators = inapResult.Q;

        return result;
    }

    private static MAMResult createEmptyMAMResult(int M, int K) {
        MAMResult result = new MAMResult();
        result.QN = new Matrix(M, K);
        result.UN = new Matrix(M, K);
        result.RN = new Matrix(M, K);
        result.TN = new Matrix(M, K);
        result.CN = new Matrix(1, K);
        result.XN = new Matrix(1, K);
        result.iter = 0;
        result.method = "inap";
        return result;
    }
}
