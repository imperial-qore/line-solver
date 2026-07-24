/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mam.analyzers;

import java.util.Map;

import jline.api.qsys.Qsys_mmck;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.solvers.mam.MAMResult;
import jline.util.matrix.Matrix;

/**
 * Exact fast-path for the finite-capacity M/M/c/K queue: a Source-Queue(-Sink)
 * model whose station has a finite total capacity, Poisson arrivals and a
 * shared exponential service rate. The decomposition methods ignore the finite
 * capacity entirely (returning the unstable/unconstrained answer), while the
 * truncated Erlang closed form is exact. Mirrors the MATLAB
 * mam_detect_mmck/qsys_mmck dispatch (solver_mam.m) and the python-native MAM
 * finite-capacity branch.
 */
public final class Solver_mam_mmck_exact {
    private Solver_mam_mmck_exact() {}

    /**
     * @return an exact MAMResult, or null when the model is not M/M/c/K
     */
    public static MAMResult solver_mam_mmck_exact(NetworkStruct sn) {
        int K = sn.nclasses;
        int M = sn.nstations;
        if (M != 2) {
            return null;
        }
        for (int r = 0; r < K; r++) {
            if (!Double.isInfinite(sn.njobs.get(r))) {
                return null; // closed/mixed models not covered
            }
        }
        int sourceIst = -1;
        int queueIst = -1;
        for (int i = 0; i < M; i++) {
            SchedStrategy sched = sn.sched.get(sn.stations.get(i));
            if (sched == SchedStrategy.EXT) {
                sourceIst = i;
            } else if (sched == SchedStrategy.FCFS || sched == SchedStrategy.HOL) {
                queueIst = i;
            }
        }
        if (sourceIst < 0 || queueIst < 0) {
            return null;
        }
        // see _kb/06-solver-catalog.md for rationale
        if (sn.droprule != null) {
            Map<JobClass, DropStrategy> dropMap = sn.droprule.get(sn.stations.get(queueIst));
            if (dropMap != null) {
                for (DropStrategy ds : dropMap.values()) {
                    if (ds == DropStrategy.Retrial || ds == DropStrategy.RetrialWithLimit) {
                        return null;
                    }
                }
            }
        }

        double capK = sn.cap.get(queueIst);
        // see _kb/06-solver-catalog.md for rationale
        if (Double.isInfinite(capK) || capK <= 0 || capK >= 1e6) {
            return null;
        }
        double sD = sn.nservers.get(queueIst);
        if (Double.isInfinite(sD) || sD < 1) {
            return null;
        }
        int c = (int) sD;
        if (capK < c) {
            return null;
        }

        // Every active class must be Poisson-in / Exp-out with a shared rate
        double muShared = Double.NaN;
        double lambdaTot = 0.0;
        double[] lam = new double[K];
        for (int r = 0; r < K; r++) {
            JobClass jc = sn.jobclasses.get(r);
            Map<JobClass, ProcessType> srcMap = sn.procid.get(sn.stations.get(sourceIst));
            Map<JobClass, ProcessType> qMap = sn.procid.get(sn.stations.get(queueIst));
            ProcessType arv = srcMap != null ? srcMap.get(jc) : null;
            ProcessType svc = qMap != null ? qMap.get(jc) : null;
            double lr = sn.rates.get(sourceIst, r);
            if (Double.isNaN(lr) || lr <= 0) {
                continue; // disabled class
            }
            if (arv != ProcessType.EXP || svc != ProcessType.EXP) {
                return null;
            }
            double mr = sn.rates.get(queueIst, r);
            if (Double.isNaN(mr) || mr <= 0) {
                return null;
            }
            if (Double.isNaN(muShared)) {
                muShared = mr;
            } else if (Math.abs(mr - muShared) > 1e-9 * Math.max(1.0, muShared)) {
                return null; // per-class service rates differ
            }
            lam[r] = lr;
            lambdaTot += lr;
        }
        if (Double.isNaN(muShared) || lambdaTot <= 0) {
            return null;
        }

        Qsys_mmck.Result ex = Qsys_mmck.qsys_mmck(lambdaTot, muShared, c, (int) capK);

        Matrix QN = new Matrix(M, K);
        Matrix UN = new Matrix(M, K);
        Matrix RN = new Matrix(M, K);
        Matrix TN = new Matrix(M, K);
        Matrix CN = new Matrix(M, K);
        Matrix XN = new Matrix(1, K);
        for (int r = 0; r < K; r++) {
            if (lam[r] <= 0) continue;
            double share = lam[r] / lambdaTot;
            double tputR = ex.throughput * share;
            TN.set(sourceIst, r, tputR);
            TN.set(queueIst, r, tputR);
            XN.set(0, r, tputR);
            QN.set(queueIst, r, ex.meanQueueLength * share);
            UN.set(queueIst, r, ex.utilization * share);
            RN.set(queueIst, r, ex.meanSojournTime);
            CN.set(queueIst, r, ex.meanSojournTime);
        }

        MAMResult result = new MAMResult();
        result.QN = QN;
        result.UN = UN;
        result.RN = RN;
        result.TN = TN;
        result.CN = CN;
        result.XN = XN;
        result.method = "exact:mmck";
        result.iter = 0;
        return result;
    }
}
