/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 *
 * Exact MAP/MAP/1 solution for a single-class, single-server, open
 * Source -> FCFS Queue -> Sink model. The decomposition methods approximate
 * this queue (fitting the arrival to a simpler process or treating the service
 * as a renewal phase-type); when the arrival or service is a genuinely
 * correlated (non-renewal) MAP, this routine returns the exact mean queue
 * length via the matrix-geometric MAP/MAP/1 solver (Q_CT_MAP_MAP_1).
 */
package jline.solvers.mam.analyzers;

import java.util.Map;

import jline.api.mam.Map_lambda;
import jline.api.mam.Map_pie;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.lib.qmam.MAPMAP1Result;
import jline.lib.qmam.Q_CT_MAP_MAP_1;
import jline.solvers.mam.MAMResult;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Solver_mam_mapmap1_exact {

    private Solver_mam_mapmap1_exact() {}

    /**
     * Returns an exact MAMResult when {@code sn} is a single-class single-server
     * open MAP/MAP/1 queue with a correlated (non-renewal) arrival or service;
     * otherwise returns {@code null} so the caller falls back to decomposition.
     */
    public static MAMResult solver_mam_mapmap1_exact(NetworkStruct sn) {
        int M = sn.nstations;
        int K = sn.nclasses;
        if (K != 1) {
            return null;
        }
        if (!Double.isInfinite(sn.njobs.get(0))) {
            return null;   // open model only
        }

        int sourceIdx = -1;
        int queueIdx = -1;
        for (int i = 0; i < M; i++) {
            SchedStrategy sched = sn.sched.get(sn.stations.get(i));
            if (sched == SchedStrategy.EXT) {
                sourceIdx = i;
            } else if (sched == SchedStrategy.FCFS) {
                queueIdx = i;
            } else {
                return null;   // any other station type breaks the single-queue assumption
            }
        }
        if (sourceIdx < 0 || queueIdx < 0) {
            return null;
        }
        if ((int) sn.nservers.get(queueIdx, 0) != 1) {
            return null;
        }

        Matrix[] arr = readMap(sn, sourceIdx);
        Matrix[] svc = readMap(sn, queueIdx);
        if (arr == null || svc == null) {
            return null;
        }
        Matrix Da0 = arr[0], Da1 = arr[1];
        Matrix Ds0 = svc[0], Ds1 = svc[1];

        // see _kb/06-solver-catalog.md for rationale
        if (!isMarkovianMap(Da0, Da1) || !isMarkovianMap(Ds0, Ds1)) {
            return null;
        }

        // Only fire when the decomposition methods are inexact: a genuinely
        // correlated (non-renewal) MAP arrival or service.
        if (isRenewalMap(Da0, Da1) && isRenewalMap(Ds0, Ds1)) {
            return null;
        }

        double lambda = Map_lambda.map_lambda(Da0, Da1);
        double mu = Map_lambda.map_lambda(Ds0, Ds1);
        if (!(lambda < mu)) {
            return null;   // unstable or degenerate: leave to the fallback path
        }

        MAPMAP1Result mm1 = Q_CT_MAP_MAP_1.qCtMapMap1(Da0, Da1, Ds0, Ds1);
        Matrix ql = mm1.getQueueLength();
        double EN = 0.0;
        int n = ql.getNumElements();
        for (int kk = 0; kk < n; kk++) {
            EN += kk * ql.get(kk);
        }
        double rho = lambda / mu;

        MAMResult result = new MAMResult();
        Matrix QN = new Matrix(M, K);
        Matrix UN = new Matrix(M, K);
        Matrix RN = new Matrix(M, K);
        Matrix TN = new Matrix(M, K);
        TN.set(sourceIdx, 0, lambda);
        TN.set(queueIdx, 0, lambda);
        QN.set(queueIdx, 0, EN);
        UN.set(queueIdx, 0, rho);
        RN.set(queueIdx, 0, EN / lambda);
        Matrix CN = new Matrix(1, K);
        CN.set(0, 0, EN / lambda);
        Matrix XN = new Matrix(1, K);
        XN.set(0, 0, lambda);

        result.QN = QN;
        result.UN = UN;
        result.RN = RN;
        result.TN = TN;
        result.CN = CN;
        result.XN = XN;
        result.iter = 1;
        result.method = "exact.mapmap1";
        return result;
    }

    private static Matrix[] readMap(NetworkStruct sn, int idx) {
        Station station = sn.stations.get(idx);
        JobClass jobClass = sn.jobclasses.get(0);
        Map<JobClass, MatrixCell> procMap = sn.proc.get(station);
        if (procMap == null) {
            return null;
        }
        MatrixCell cell = procMap.get(jobClass);
        if (cell == null) {
            return null;
        }
        Matrix D0 = cell.get(0);
        Matrix D1 = cell.get(1);
        if (D0 == null || D1 == null) {
            return null;
        }
        return new Matrix[]{D0, D1};
    }

    /**
     * True when (D0, D1) is a genuine MAP rather than a RAP or an ME process:
     * non-negative off-diagonal rates in D0, non-negative rates in D1, and
     * (D0 + D1) an infinitesimal generator (zero row sums). A RAP or an ME
     * violates the sign conditions while still defining a valid point process.
     */
    private static boolean isMarkovianMap(Matrix D0, Matrix D1) {
        int ns = D0.getNumRows();
        double scale = 1.0;
        for (int i = 0; i < ns; i++) {
            for (int j = 0; j < ns; j++) {
                scale = Math.max(scale, Math.abs(D0.get(i, j)));
                scale = Math.max(scale, Math.abs(D1.get(i, j)));
            }
        }
        double tol = 1e-9 * scale;
        for (int i = 0; i < ns; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < ns; j++) {
                if (i != j && D0.get(i, j) < -tol) {
                    return false;
                }
                if (D1.get(i, j) < -tol) {
                    return false;
                }
                rowSum += D0.get(i, j) + D1.get(i, j);
            }
            if (Math.abs(rowSum) > tol) {
                return false;
            }
        }
        return true;
    }

    private static boolean isRenewalMap(Matrix D0, Matrix D1) {
        int ns = D0.getNumRows();
        if (ns == 1) {
            return true;
        }
        Matrix tExit = D1.mult(Matrix.ones(ns, 1));
        Matrix alpha = Map_pie.map_pie(D0, D1);   // 1 x ns
        double diff = 0.0, ref = 0.0;
        for (int i = 0; i < ns; i++) {
            for (int j = 0; j < ns; j++) {
                double d = D1.get(i, j) - tExit.get(i, 0) * alpha.get(0, j);
                diff += d * d;
                ref += D1.get(i, j) * D1.get(i, j);
            }
        }
        return Math.sqrt(diff) < 1e-9 * Math.max(1.0, Math.sqrt(ref));
    }
}
