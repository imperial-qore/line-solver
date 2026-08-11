/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mam.handlers;

import java.util.ArrayList;
import java.util.List;

import jline.api.mam.Map_mean;
import jline.api.mam.Map_pie;
import jline.api.sn.SnHasLoadDependence;
import jline.io.InputOutput;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * LD-QBD block exposure, flattening and metric mapping for the SolverENV
 * state-vector analyzer's MAM backend.
 *
 * <p>Mirrors matlab/src/solvers/MAM/solver_mam_ldqbd.m (block construction and
 * the {@code ld} struct), solver_mam_ldqbd_flatten.m and solver_mam_ldqbd_avg.m.
 * Handles single-class Delay+Queue (closed) or Source+Queue (open) models, with
 * exact M/M/c boundary, PH service, and load-dependent scaling. The open regime
 * truncates the level space (options.cutoff or a negligible-tail bound).
 *
 * <p>The existing steady-state {@link Solver_mam_ldqbd} closed-path solver is
 * left untouched; this class provides only the block/flatten/avg pieces the
 * state-vector analyzer needs.
 */
public final class Solver_mam_ldqbd_statevec {
    private Solver_mam_ldqbd_statevec() {}

    private static final String MFILENAME = "solver_mam_ldqbd";

    /** Block-tridiagonal LD-QBD representation plus the parameters the avg mapping needs. */
    public static final class Ld {
        public List<Matrix> Q0;   // size Nlev : upward (arrival) blocks
        public List<Matrix> Q1;   // size Nlev+1 : local blocks
        public List<Matrix> Q2;   // size Nlev : downward (departure) blocks
        public int Nlev;
        public int nPhases;
        public boolean isPH;
        public boolean isOpen;
        public int queueIdx;
        public int refIdx;
        public int M;
        public int nServers;
        public double mean_service;
        public boolean hasLLD;
        public double lambda_eff;
        public double delayRate;  // NaN for open
        public double N;          // Inf for open
    }

    /** Flattened generator together with the queue level of each flat state. */
    public static final class Flat {
        public Matrix Q;
        public int[] levelOf;
    }

    /** Per-(station,class) mean metrics derived from a flat LD-QBD distribution. */
    public static final class Avg {
        public Matrix QN;
        public Matrix UN;
        public Matrix RN;
        public Matrix TN;
    }

    /**
     * Build the LD-QBD blocks and parameters (the {@code ld} struct) for a
     * single-class Delay/Queue (closed) or Source/Queue (open) model.
     */
    public static Ld solver_mam_ldqbd_ld(NetworkStruct sn, SolverOptions options) {
        int M = sn.nstations;
        int K = sn.nclasses;
        if (K != 1) {
            InputOutput.line_error(MFILENAME, "LDQBD method requires a single-class model.");
        }
        double Npop = sn.njobs.get(0);
        boolean isOpen = !Double.isFinite(Npop);

        int nDelay = 0, nQueue = 0, nSource = 0;
        int delayIdx = -1, queueIdx = -1, srcIdx = -1;
        for (int i = 0; i < M; i++) {
            SchedStrategy sched = sn.sched.get(sn.stations.get(i));
            if (sched == SchedStrategy.INF) { nDelay++; delayIdx = i; }
            else if (sched == SchedStrategy.FCFS) { nQueue++; queueIdx = i; }
            else if (sched == SchedStrategy.EXT) { nSource++; srcIdx = i; }
        }
        if (isOpen) {
            if (nSource != 1 || nQueue != 1 || M != 2) {
                InputOutput.line_error(MFILENAME, "Open LDQBD method requires exactly one Source and one Queue station.");
            }
        } else {
            if (nDelay != 1 || nQueue != 1 || M != 2) {
                InputOutput.line_error(MFILENAME, "Closed LDQBD method requires exactly one Delay and one Queue station.");
            }
        }

        Matrix rates = sn.rates;
        Matrix nservers = sn.nservers;
        Station queueStation = sn.stations.get(queueIdx);
        JobClass jobClass = sn.jobclasses.get(0);
        MatrixCell PH_queue = sn.proc.get(queueStation).get(jobClass);
        int nServers = (int) nservers.get(queueIdx, 0);

        Matrix D0 = PH_queue.get(0);
        Matrix D1 = (PH_queue.size() > 1) ? PH_queue.get(1) : null;
        boolean isPH = !(D0.getNumRows() == 1 && D0.getNumCols() == 1);
        double mu;
        int nPhases;
        double mean_service;
        Matrix alpha = null;
        if (!isPH) {
            mu = -D0.get(0, 0);
            nPhases = 1;
            mean_service = 1.0 / mu;
        } else {
            mu = Double.NaN;
            nPhases = D0.getNumRows();
            alpha = Map_pie.map_pie(PH_queue);
            mean_service = Map_mean.map_mean(PH_queue);
        }

        // Load-dependent per-level service factor (else min(n,c)).
        boolean hasLLD = SnHasLoadDependence.snHasLoadDependence(sn)
                && sn.lldscaling != null && !sn.lldscaling.isEmpty()
                && sn.lldscaling.getNumRows() > queueIdx
                && lldRowHasNonUnit(sn.lldscaling, queueIdx);
        double[] lld = null;
        int lldlimit = 0;
        double sfMax;
        if (hasLLD) {
            lldlimit = sn.lldscaling.getNumCols();
            lld = new double[lldlimit];
            for (int j = 0; j < lldlimit; j++) lld[j] = sn.lldscaling.get(queueIdx, j);
            sfMax = lld[lldlimit - 1];
        } else {
            sfMax = nServers;
        }

        Matrix rt = sn.rt;
        int Nlev;
        double lambda_eff;
        double[] arrRate;   // arrRate[n] = rate out of level n
        if (isOpen) {
            MatrixCell arrProc = sn.proc.get(sn.stations.get(srcIdx)).get(jobClass);
            if (arrProc.get(0).getNumRows() > 1) {
                InputOutput.line_error(MFILENAME, "Open LDQBD method currently supports Poisson (exponential) arrivals only; the Source uses a MAP/MMPP process.");
            }
            double lambda = rates.get(srcIdx, 0);
            lambda_eff = lambda * rt.get(srcIdx, queueIdx);
            double rho = lambda_eff * mean_service / sfMax;
            if (rho >= 1) {
                InputOutput.line_error(MFILENAME, String.format("Open LDQBD method requires a stable queue (rho = %.4f >= 1). Increase service capacity or reduce the arrival rate.", rho));
            }
            if (options.cutoff != null && !options.cutoff.isEmpty() && Double.isFinite(options.cutoff.get(0))) {
                Nlev = Math.max(nServers + 1, (int) Math.round(options.cutoff.get(0)));
            } else {
                double tailTol = 1e-10;
                Nlev = nServers + (int) Math.ceil(Math.log(tailTol) / Math.log(rho));
                Nlev = Math.min(Math.max(Nlev, nServers + 10), 100000);
            }
            arrRate = new double[Nlev + 1];
            for (int n = 0; n <= Nlev; n++) arrRate[n] = lambda_eff;
            arrRate[Nlev] = 0;
        } else {
            int N = (int) Npop;
            double lambda_d = rates.get(delayIdx, 0);
            lambda_eff = lambda_d * rt.get(delayIdx, queueIdx);
            Nlev = N;
            arrRate = new double[Nlev + 1];
            for (int n = 0; n <= Nlev; n++) arrRate[n] = (N - n) * lambda_eff;
        }

        // Per-level service factor sf(n), n = 1..Nlev (stored at index n-1).
        double[] sf = new double[Nlev];
        for (int n = 1; n <= Nlev; n++) {
            if (hasLLD) sf[n - 1] = lld[Math.min(n, lldlimit) - 1];
            else sf[n - 1] = Math.min(n, nServers);
        }

        List<Matrix> Q0 = new ArrayList<Matrix>();
        List<Matrix> Q1 = new ArrayList<Matrix>();
        List<Matrix> Q2 = new ArrayList<Matrix>();
        if (!isPH) {
            for (int n = 0; n < Nlev; n++) {
                Matrix b = new Matrix(1, 1); b.set(0, 0, arrRate[n]); Q0.add(b);
            }
            for (int n = 0; n <= Nlev; n++) {
                double dep = (n > 0) ? sf[n - 1] * mu : 0.0;
                Matrix b = new Matrix(1, 1); b.set(0, 0, -(arrRate[n] + dep)); Q1.add(b);
            }
            for (int n = 1; n <= Nlev; n++) {
                Matrix b = new Matrix(1, 1); b.set(0, 0, sf[n - 1] * mu); Q2.add(b);
            }
        } else {
            Q0.add(alpha.scale(arrRate[0]));                              // level 0 -> 1
            for (int n = 1; n < Nlev; n++) {
                Q0.add(Matrix.eye(nPhases).scale(arrRate[n]));           // level n -> n+1
            }
            Matrix Q1_0 = new Matrix(1, 1); Q1_0.set(0, 0, -arrRate[0]); // level 0
            Q1.add(Q1_0);
            for (int n = 1; n <= Nlev; n++) {
                Q1.add(D0.scale(sf[n - 1]).sub(Matrix.eye(nPhases).scale(arrRate[n])));
            }
            Q2.add(D1.scale(sf[0]).mult(Matrix.ones(nPhases, 1)));        // level 1 -> 0
            for (int n = 2; n <= Nlev; n++) {
                Q2.add(D1.scale(sf[n - 1]));                              // level n -> n-1
            }
        }

        Ld ld = new Ld();
        ld.Q0 = Q0; ld.Q1 = Q1; ld.Q2 = Q2;
        ld.Nlev = Nlev; ld.nPhases = nPhases; ld.isPH = isPH; ld.isOpen = isOpen;
        ld.queueIdx = queueIdx; ld.M = M; ld.nServers = nServers;
        ld.mean_service = mean_service; ld.hasLLD = hasLLD; ld.lambda_eff = lambda_eff;
        if (isOpen) {
            ld.refIdx = srcIdx; ld.delayRate = Double.NaN; ld.N = Double.POSITIVE_INFINITY;
        } else {
            ld.refIdx = delayIdx; ld.delayRate = rates.get(delayIdx, 0); ld.N = Npop;
        }
        return ld;
    }

    private static boolean lldRowHasNonUnit(Matrix lldscaling, int row) {
        for (int j = 0; j < lldscaling.getNumCols(); j++) {
            if (lldscaling.get(row, j) != 1) return true;
        }
        return false;
    }

    /** Assemble a dense generator from the block-tridiagonal LD-QBD representation. */
    public static Flat solver_mam_ldqbd_flatten(Ld ld) {
        int Nlev = ld.Nlev;
        int[] levelSize = new int[Nlev + 1];
        for (int n = 0; n <= Nlev; n++) levelSize[n] = ld.Q1.get(n).getNumRows();
        int[] levelStart = new int[Nlev + 2];
        for (int n = 0; n <= Nlev; n++) levelStart[n + 1] = levelStart[n] + levelSize[n];
        int dim = levelStart[Nlev + 1];

        Matrix Q = new Matrix(dim, dim);
        Q.zero();
        int[] levelOf = new int[dim];
        for (int n = 0; n <= Nlev; n++) {
            int r0 = levelStart[n];
            for (int a = 0; a < levelSize[n]; a++) levelOf[r0 + a] = n;
            setBlock(Q, r0, r0, ld.Q1.get(n));                 // within-level
            if (n < Nlev) setBlock(Q, r0, levelStart[n + 1], ld.Q0.get(n));   // upward
            if (n >= 1) setBlock(Q, r0, levelStart[n - 1], ld.Q2.get(n - 1)); // downward
        }
        Flat flat = new Flat();
        flat.Q = Q; flat.levelOf = levelOf;
        return flat;
    }

    private static void setBlock(Matrix Q, int r0, int c0, Matrix block) {
        for (int a = 0; a < block.getNumRows(); a++) {
            for (int b = 0; b < block.getNumCols(); b++) {
                Q.set(r0 + a, c0 + b, block.get(a, b));
            }
        }
    }

    /** Map a flat LD-QBD distribution to per-(station,class) mean metrics. */
    public static Avg solver_mam_ldqbd_avg(Ld ld, Matrix piflat, int[] levelOf) {
        int Nlev = ld.Nlev;
        int M = ld.M;
        int qi = ld.queueIdx;
        int ri = ld.refIdx;
        int c = ld.nServers;

        double[] pf = piflat.toArray1D();
        double psum = 0.0;
        for (int i = 0; i < pf.length; i++) { if (pf[i] < 0) pf[i] = 0; psum += pf[i]; }
        if (psum > 0) for (int i = 0; i < pf.length; i++) pf[i] /= psum;

        double[] pLevel = new double[Nlev + 1];
        for (int i = 0; i < pf.length; i++) pLevel[levelOf[i]] += pf[i];

        double mean_queue = 0.0;
        for (int n = 0; n <= Nlev; n++) mean_queue += n * pLevel[n];

        double util;
        if (ld.hasLLD || c == 1) {
            util = 1 - pLevel[0];
        } else {
            util = 0;
            for (int n = 1; n <= Nlev; n++) util += (Math.min(n, c) / (double) c) * pLevel[n];
        }

        Matrix QN = new Matrix(M, 1); QN.zero();
        Matrix UN = new Matrix(M, 1); UN.zero();
        Matrix RN = new Matrix(M, 1); RN.zero();
        Matrix TN = new Matrix(M, 1); TN.zero();

        if (ld.isOpen) {
            double X = ld.lambda_eff * (1 - pLevel[Nlev]);
            double R_queue = (X > 0) ? mean_queue / X : 0.0;
            QN.set(ri, 0, 0); UN.set(ri, 0, 0); RN.set(ri, 0, 0); TN.set(ri, 0, X);
            QN.set(qi, 0, mean_queue); UN.set(qi, 0, util); RN.set(qi, 0, R_queue); TN.set(qi, 0, X);
        } else {
            double mean_delay = ld.N - mean_queue;
            double X = mean_delay * ld.lambda_eff;
            double R_queue = (X > 0) ? mean_queue / X : 0.0;
            double R_delay = 1.0 / ld.delayRate;
            QN.set(ri, 0, mean_delay); UN.set(ri, 0, mean_delay); RN.set(ri, 0, R_delay); TN.set(ri, 0, X);
            QN.set(qi, 0, mean_queue); UN.set(qi, 0, util / c); RN.set(qi, 0, R_queue); TN.set(qi, 0, X);
        }
        Avg avg = new Avg();
        avg.QN = QN; avg.UN = UN; avg.RN = RN; avg.TN = TN;
        return avg;
    }
}
