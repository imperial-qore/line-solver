/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.ag.handlers;

import jline.solvers.mam.handlers.MetricsResult;

import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.ag.AGOptions;
import jline.solvers.ag.AGResult;
import jline.solvers.ag.AgExec;
import jline.util.matrix.Matrix;

import static jline.io.InputOutput.line_warning;

/**
 * Entry point of the agent-based (RCAT) solver.
 *
 * <p>Builds the agents out of the network, runs the reversed-rate fixed point
 * through the requested execution backend, and converts the converged agents
 * into the mean measures.</p>
 */
public final class Solver_ag {
    private Solver_ag() {}

    private static final String MFILENAME = "solver_ag";

    public static AGResult solver_ag(NetworkStruct sn, SolverOptions options) {
        int M = sn.nstations;
        int K = sn.nclasses;

        int maxStates = (options instanceof AGOptions) ? ((AGOptions) options).maxStates : 100;

        double tol = options.iter_tol > 0 ? options.iter_tol : 1e-6;
        int maxiter = options.iter_max > 0 ? options.iter_max : 1000;

        RCATModel rcat = Solver_ag_build.solver_ag_build(sn, maxStates);

        int numProcesses = rcat.N.length;
        int numActions = rcat.actionMap.size();

        if (numProcesses == 0 || (numActions == 0 && numProcesses > 1)) {
            line_warning(MFILENAME, "Network could not be mapped to RCAT format.");
            return createEmptyAGResult(M, K);
        }

        String method = options.method;
        if ("default".equals(method)) {
            method = "inap";
        }

        AGOptions agOptions = (options instanceof AGOptions) ? (AGOptions) options : null;
        if (agOptions != null && AGOptions.EXEC_CLUSTER.equals(agOptions.exec)
                && "inapinf".equals(method)) {
            // The remote worker implements the FINITE agent solve. 'inapinf'
            // replaces it with the matrix-geometric treatment of an open agent --
            // Neuts' R matrix and the scalar-tail detection that precedes it --
            // which the worker does not carry, and answering with the finite
            // solve instead would silently change the method.
            throw new RuntimeException("The 'cluster' execution backend does not carry the "
                    + "'inapinf' agent solve (the matrix-geometric tail of an open agent runs "
                    + "on the coordinator only). Use exec 'serial' or 'parallel' with 'inapinf', "
                    + "or method 'inap'/'inapplus' with 'cluster'.");
        }

        INAPResult inapResult;
        AgExec exec = AgExec.create(agOptions);
        try {
        if ("inap".equals(method)) {
            inapResult = Solver_ag_inap.solver_ag_inap(rcat, tol, maxiter, options.seed, "inap", exec);
        } else if ("inapplus".equals(method)) {
            inapResult = Solver_ag_inap.solver_ag_inap(rcat, tol, maxiter, options.seed, "inapplus", exec);
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
            inapResult = Solver_ag_inap.solver_ag_inapinf(rcat, isOpenProc, tol, maxiter);
        } else if ("exact".equals(method)) {
            line_warning(MFILENAME, "AutoCAT (exact method) is not available. Falling back to INAP.");
            inapResult = Solver_ag_inap.solver_ag_inap(rcat, tol, maxiter, options.seed, "inap", exec);
        } else {
            line_warning(MFILENAME, "Unknown method: " + method + ". Using INAP.");
            inapResult = Solver_ag_inap.solver_ag_inap(rcat, tol, maxiter, options.seed, "inap", exec);
        }
        } finally {
            if (exec != null) exec.close();
        }

        MetricsResult metrics = Solver_ag_metrics.solver_ag_metrics(sn, inapResult, rcat);

        AGResult result = new AGResult();
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

    private static AGResult createEmptyAGResult(int M, int K) {
        AGResult result = new AGResult();
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
