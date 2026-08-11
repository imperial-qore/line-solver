/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mam.handlers;

import java.util.ArrayList;
import java.util.List;

import jline.api.mam.Ldqbd;
import jline.api.mam.LdqbdOptions;
import jline.api.mam.LdqbdResult;
import jline.api.mam.Map_mean;
import jline.api.mam.Map_pie;
import jline.io.InputOutput;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.JobClass;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.solvers.mam.MAMResult;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Solver for single-class closed queueing networks using Level-Dependent QBD.
 *
 * Uses Level-Dependent Quasi-Birth-Death (LD-QBD) process to compute
 * performance metrics for single-class closed queueing networks
 * consisting of a Delay (infinite server) and a Queue (FCFS).
 *
 * Exactness: exact for exponential service at any number of servers, and for
 * PH service at a single server. For PH service with c &gt; 1 servers it is an
 * approximation: the c parallel PH servers are collapsed into one PH process
 * scaled by min(n,c), which ignores the phase of each individual busy server
 * (the exact chain tracks the multiset of the min(n,c) in-service phases).
 *
 * The LD-QBD approach models the system where:
 *   - Level n = number of jobs at the queue (0 &lt;= n &lt;= N)
 *   - Jobs at delay = N - n
 *   - Transition rates depend on the current level
 *
 * Supports PH-type service distributions (Exp, Erlang, HyperExp, etc.)
 */
public final class Solver_mam_ldqbd {
    private Solver_mam_ldqbd() {}

    private static final String MFILENAME = "solver_mam_ldqbd";

    public static MAMResult solver_mam_ldqbd(NetworkStruct sn, SolverOptions options) {
        int M = sn.nstations;
        int K = sn.nclasses;

        // Check: single-class closed network
        if (K != 1) {
            InputOutput.line_error(MFILENAME, "LDQBD method requires a single-class model.");
            return createEmptyResult(M, K);
        }

        int N = (int) sn.njobs.get(0, 0);
        if (!Double.isFinite((double) N) || N <= 0) {
            InputOutput.line_error(MFILENAME, "LDQBD method requires a closed model with finite population.");
            return createEmptyResult(M, K);
        }

        // Check: must have exactly one delay and one queue
        int nDelay = 0;
        int nQueue = 0;
        int delayIdx = -1;
        int queueIdx = -1;

        for (int i = 0; i < M; i++) {
            SchedStrategy sched = sn.sched.get(sn.stations.get(i));
            if (sched == SchedStrategy.INF) {
                nDelay++;
                delayIdx = i;
            } else if (sched == SchedStrategy.FCFS) {
                nQueue++;
                queueIdx = i;
            }
        }

        if (nDelay != 1 || nQueue != 1 || M != 2) {
            InputOutput.line_error(MFILENAME, "LDQBD method requires exactly one Delay and one Queue station.");
            return createEmptyResult(M, K);
        }

        // Get service parameters
        Matrix rates = sn.rates;
        Matrix nservers = sn.nservers;

        // see _kb/06-solver-catalog.md for rationale
        double lambda_d = rates.get(delayIdx, 0);
        double lambda_eff = lambda_d * sn.rt.get(delayIdx, queueIdx);

        // Get station and job class objects for proc access
        Station queueStation = sn.stations.get(queueIdx);
        JobClass jobClass = sn.jobclasses.get(0);

        // Queue service process
        if (sn.proc == null || sn.proc.get(queueStation) == null
                || sn.proc.get(queueStation).get(jobClass) == null) {
            throw new RuntimeException("No service process for queue station");
        }
        MatrixCell PH_queue = sn.proc.get(queueStation).get(jobClass);
        int nServers = (int) nservers.get(queueIdx, 0);

        // Check if queue service is exponential (1x1 matrix) or PH
        if (PH_queue.get(0) == null) {
            throw new RuntimeException("No D0 matrix in service process");
        }
        Matrix D0 = PH_queue.get(0);
        boolean isExponential = D0.getNumRows() == 1 && D0.getNumCols() == 1;
        double mu;
        int nPhases;

        if (isExponential) {
            mu = -D0.get(0, 0);
            nPhases = 1;
        } else {
            nPhases = D0.getNumRows();
            mu = 1.0 / Map_mean.map_mean(PH_queue);
        }

        // Build LD-QBD matrices
        List<Matrix> Q0 = new ArrayList<Matrix>();
        List<Matrix> Q1 = new ArrayList<Matrix>();
        List<Matrix> Q2 = new ArrayList<Matrix>();

        if (isExponential) {
            // Exponential service case (scalar matrices)
            for (int n = 0; n < N; n++) {
                Matrix arrRate = new Matrix(1, 1);
                arrRate.set(0, 0, (N - n) * lambda_eff);
                Q0.add(arrRate);
            }

            for (int n = 0; n <= N; n++) {
                double arrivalRate = (N - n) * lambda_eff;
                double departureRate = (n > 0) ? Math.min(n, nServers) * mu : 0.0;
                Matrix local = new Matrix(1, 1);
                local.set(0, 0, -(arrivalRate + departureRate));
                Q1.add(local);
            }

            for (int n = 1; n <= N; n++) {
                Matrix depRate = new Matrix(1, 1);
                depRate.set(0, 0, Math.min(n, nServers) * mu);
                Q2.add(depRate);
            }
        } else {
            // PH-type service case (matrix-valued)
            if (PH_queue.get(1) == null) {
                throw new RuntimeException("No D1 matrix in PH service process");
            }
            Matrix D1 = PH_queue.get(1);
            Matrix alpha = Map_pie.map_pie(PH_queue);

            // Level 0 -> 1: arrival starts service in some phase
            Q0.add(alpha.scale(N * lambda_eff));

            // Level n -> n+1 (n >= 1): arrival, preserve phase
            for (int n = 1; n < N; n++) {
                Q0.add(Matrix.eye(nPhases).scale((N - n) * lambda_eff));
            }

            // Level 0: 1x1 (no service, only arrivals)
            Matrix Q1_0 = new Matrix(1, 1);
            Q1_0.set(0, 0, -(N * lambda_eff));
            Q1.add(Q1_0);

            // Level n >= 1: phase transitions within level
            for (int n = 1; n <= N; n++) {
                double arrivalRate = (N - n) * lambda_eff;
                int c_n = Math.min(n, nServers);
                Matrix local = D0.scale((double) c_n).sub(Matrix.eye(nPhases).scale(arrivalRate));
                Q1.add(local);
            }

            // Level 1 -> 0: service completion, go to empty state
            int c_1 = Math.min(1, nServers);
            Matrix Q2_1 = D1.scale((double) c_1).mult(Matrix.ones(nPhases, 1));
            Q2.add(Q2_1);

            // Level n -> n-1 (n >= 2): service completion, next job starts
            for (int n = 2; n <= N; n++) {
                int c_n = Math.min(n, nServers);
                Q2.add(D1.scale((double) c_n));
            }
        }

        // Solve LD-QBD
        LdqbdOptions ldqbdOptions = new LdqbdOptions(options.tol, options.iter_max, false);

        LdqbdResult ldqbdResult = Ldqbd.ldqbd(Q0, Q1, Q2, ldqbdOptions);
        Matrix pi_ldqbd = ldqbdResult.getPi();

        // Compute performance metrics from steady-state distribution
        double mean_queue = 0.0;
        for (int n = 0; n <= N; n++) {
            mean_queue += n * pi_ldqbd.get(0, n);
        }

        double mean_delay = N - mean_queue;

        // Queue flow X = mean_delay * lambda_eff (flow balance into the queue)
        double X = mean_delay * lambda_eff;

        // Mean service time at queue
        double mean_service = isExponential ? 1.0 / mu : Map_mean.map_mean(PH_queue);

        // see _kb/06-solver-catalog.md for rationale
        double util_queue;
        if (nServers == 1) {
            util_queue = 1.0 - pi_ldqbd.get(0, 0);
        } else {
            double u = 0.0;
            for (int n = 1; n <= N; n++) {
                u += ((double) Math.min(n, nServers) / nServers) * pi_ldqbd.get(0, n);
            }
            util_queue = u;
        }

        // Response time at queue (using Little's law: R = Q/X)
        double R_queue = (X > 0) ? mean_queue / X : 0.0;

        // Response time at delay (constant for infinite server)
        double R_delay = 1.0 / lambda_d;

        // Populate output matrices
        MAMResult result = new MAMResult();
        result.QN = new Matrix(M, K);
        result.UN = new Matrix(M, K);
        result.RN = new Matrix(M, K);
        result.TN = new Matrix(M, K);
        result.CN = new Matrix(1, K);
        result.XN = new Matrix(1, K);

        // see _kb/06-solver-catalog.md for rationale
        result.QN.set(delayIdx, 0, mean_delay);
        result.UN.set(delayIdx, 0, mean_delay);
        result.RN.set(delayIdx, 0, R_delay);
        result.TN.set(delayIdx, 0, mean_delay * lambda_d);

        // Queue station metrics
        result.QN.set(queueIdx, 0, mean_queue);
        result.UN.set(queueIdx, 0, util_queue);
        result.RN.set(queueIdx, 0, R_queue);
        result.TN.set(queueIdx, 0, X);

        // System-level metrics
        result.XN.set(0, 0, X);
        result.CN.set(0, 0, R_delay + R_queue);

        result.iter = 1;
        result.method = "ldqbd";

        return result;
    }

    private static MAMResult createEmptyResult(int M, int K) {
        MAMResult result = new MAMResult();
        result.QN = new Matrix(M, K);
        result.UN = new Matrix(M, K);
        result.RN = new Matrix(M, K);
        result.TN = new Matrix(M, K);
        result.CN = new Matrix(1, K);
        result.XN = new Matrix(1, K);
        result.iter = 0;
        result.method = "ldqbd";
        return result;
    }
}
