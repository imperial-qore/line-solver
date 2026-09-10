/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mam.handlers;

import java.util.ArrayList;
import java.util.List;

import jline.api.mam.LdqbdMphc;
import jline.api.mam.Map_mean;
import jline.api.mam.Map_pie;
import jline.api.sn.SnHasLoadDependence;
import jline.io.InputOutput;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * LD-QBD block construction, flattening and metric mapping.
 *
 * <p>Mirrors matlab/src/solvers/MAM/solver_mam_ldqbd.m (block construction and
 * the {@code ld} struct), solver_mam_ldqbd_flatten.m and solver_mam_ldqbd_avg.m.
 * Handles single-class Delay+Queue (closed) or Source+Queue (open) models, with
 * exact M/M/c boundary, PH service at any number of servers (the busy-server
 * phase multiset of {@link LdqbdMphc}), and load-dependent scaling. The open
 * regime truncates the level space (options.cutoff or a negligible-tail bound).
 *
 * <p>This is the CANONICAL builder for both consumers: the steady-state
 * {@link Solver_mam_ldqbd} solver and the SolverENV state-vector analyzer's MAM
 * backend. Keeping one construction is what stops the two from drifting -- the
 * steady-state path carried its own copy until 2026-08-18 and silently ignored
 * sn.lldscaling and the open regime as a result.
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
        /** Per-level service factor sf(n), sf[n-1] holding level n. */
        public double[] sf;
        /**
         * The capacity that normalizes the utilization: max(c, max(alpha)), the
         * LARGEST factor the load-dependence table declares rather than the
         * saturated one, since a non-monotone alpha peaks in the middle. Same
         * rule as CTMC's ceff, which is what makes the two report the same
         * number.
         */
        public double utilPeak;
        public double lambda_eff;
        public double delayRate;  // NaN for open
        public double N;          // Inf for open
        /**
         * The station alternates OFF -&gt; setup -&gt; busy -&gt; delay-off around the
         * service, so the chain carries phases the block builder has no place
         * for. When this is set the blocks above describe a server that is
         * ALWAYS warm and must not be used: the closed regime hands the whole
         * chain to {@code Qbd_setupdelayoff_closed} instead.
         */
        public boolean hasSetup;
        public double alpharate;
        public double alphascv;
        public double betarate;
        public double betascv;
        /** Exponential service rate; NaN under phase-type service. */
        public double mu;
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
        double utilPeak;
        if (hasLLD) {
            lldlimit = sn.lldscaling.getNumCols();
            lld = new double[lldlimit];
            for (int j = 0; j < lldlimit; j++) lld[j] = sn.lldscaling.get(queueIdx, j);
            sfMax = lld[lldlimit - 1];
            utilPeak = nServers;
            for (int j = 0; j < lldlimit; j++) utilPeak = Math.max(utilPeak, lld[j]);
        } else {
            sfMax = nServers;
            utilPeak = nServers;
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
            // PH service: the level carries the MULTISET of the phases the min(n,c)
            // busy servers sit in, which is exact at any number of servers. At c == 1
            // the multiset is just the phase, so this reproduces the single-server
            // blocks (sf(n)*D0 - arr*I, sf(n)*D1) entry for entry.
            LdqbdMphc.Blocks blk =
                    LdqbdMphc.ldqbd_mphc(D0, D1, alpha, nServers, arrRate, sf);
            Q0 = blk.Q0; Q1 = blk.Q1; Q2 = blk.Q2;
        }

        Ld ld = new Ld();
        ld.Q0 = Q0; ld.Q1 = Q1; ld.Q2 = Q2;
        ld.Nlev = Nlev; ld.nPhases = nPhases; ld.isPH = isPH; ld.isOpen = isOpen;
        ld.queueIdx = queueIdx; ld.M = M; ld.nServers = nServers;
        ld.mean_service = mean_service; ld.hasLLD = hasLLD; ld.lambda_eff = lambda_eff;
        ld.sf = sf; ld.utilPeak = utilPeak; ld.mu = mu;
        // SETUP AND DELAY-OFF. Refused BY NAME outside the closed, single-server,
        // exponential, load-independent case rather than answered as if the
        // server were always warm, which is what this solver did until 2026-09
        // and is BUG-78.
        ld.hasSetup = sn.hassetup != null && sn.hassetup.getNumRows() > queueIdx
                && sn.hassetup.get(queueIdx, 0) == 1.0;
        ld.alpharate = Double.NaN; ld.alphascv = Double.NaN;
        ld.betarate = Double.NaN; ld.betascv = Double.NaN;
        if (ld.hasSetup) {
            if (isOpen) {
                throw new RuntimeException("Open LDQBD does not model a setup/delay-off server; "
                        + "use method='dec.source', whose qbd_setupdelayoff covers the open case.");
            }
            if (isPH || nServers > 1) {
                throw new RuntimeException("Closed LDQBD models a setup/delay-off server with "
                        + "exponential service at a single server only; this station has "
                        + "phase-type service or several servers.");
            }
            if (hasLLD) {
                throw new RuntimeException("Closed LDQBD models a setup/delay-off server at its "
                        + "nominal rate only; this station also declares a load-dependent scaling.");
            }
            Station setupStation = sn.stations.get(queueIdx);
            boolean got = false;
            if (setupStation instanceof Queue) {
                Queue setupQueue = (Queue) setupStation;
                for (int k = 0; k < K && !got; k++) {
                    Object setupDist = setupQueue.getSetupTime(sn.jobclasses.get(k));
                    Object delayOffDist = setupQueue.getDelayOffTime(sn.jobclasses.get(k));
                    if (setupDist != null && delayOffDist != null) {
                        try {
                            java.lang.reflect.Method getMean = setupDist.getClass().getMethod("getMean");
                            java.lang.reflect.Method getSCV = setupDist.getClass().getMethod("getSCV");
                            ld.alpharate = 1.0 / ((Number) getMean.invoke(setupDist)).doubleValue();
                            ld.alphascv = ((Number) getSCV.invoke(setupDist)).doubleValue();
                            ld.betarate = 1.0 / ((Number) getMean.invoke(delayOffDist)).doubleValue();
                            ld.betascv = ((Number) getSCV.invoke(delayOffDist)).doubleValue();
                            got = true;
                        } catch (Exception ex) {
                            got = false;
                        }
                    }
                }
            }
            ld.hasSetup = got;
        }
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

        double[] pf = piflat.toArray1D();
        double psum = 0.0;
        for (int i = 0; i < pf.length; i++) { if (pf[i] < 0) pf[i] = 0; psum += pf[i]; }
        if (psum > 0) for (int i = 0; i < pf.length; i++) pf[i] /= psum;

        double[] pLevel = new double[Nlev + 1];
        for (int i = 0; i < pf.length; i++) pLevel[levelOf[i]] += pf[i];

        double mean_queue = 0.0;
        for (int n = 0; n <= Nlev; n++) mean_queue += n * pLevel[n];

        // Utilization is the fraction of the station's PEAK capacity in use,
        // sum_n p(n)*sf(n)/utilPeak, the work-based convention CTMC, MVA, NC and
        // serial SSA all report. Without load dependence sf(n) = min(n,c) and
        // utilPeak = c, so this is the average fraction of c servers in use; at
        // c = 1 that is sf(n) = 1 for every n >= 1 and the sum collapses to
        // 1 - p(0).
        double util = 0;
        for (int n = 1; n <= Nlev; n++) util += (ld.sf[n - 1] / ld.utilPeak) * pLevel[n];

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
            // util is already per-server: the /utilPeak is inside the sum above
            QN.set(qi, 0, mean_queue); UN.set(qi, 0, util); RN.set(qi, 0, R_queue); TN.set(qi, 0, X);
        }
        Avg avg = new Avg();
        avg.QN = QN; avg.UN = UN; avg.RN = RN; avg.TN = TN;
        return avg;
    }
}
