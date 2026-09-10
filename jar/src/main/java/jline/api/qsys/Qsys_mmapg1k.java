/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.qsys;

import java.util.List;

import jline.api.mam.Map_prob;
import jline.util.matrix.Matrix;

import static jline.io.InputOutput.line_error;

/**
 * Exact per-class throughput and loss ratio of an MMAP[K]/G/1/K queue with tail
 * drop: marked Markovian arrivals, an arbitrary service law F common to all
 * classes, and a buffer of K packets counting the one in transmission.
 *
 * <p>Port of {@code matlab/src/api/qsys/qsys_mmapg1k.m}. Two classes of equal
 * arrival rate but different interarrival variability or autocorrelation receive
 * DIFFERENT loss ratios, which is the effect that motivates [1].
 * Aggregate-only finite-buffer analyses cannot express it: they return a single
 * blocking probability p and set T_k = lambda_k (1-p), making the loss ratio
 * identical across classes by construction.
 *
 * <p>METHOD. The aggregate MAP {D0, sum_k D1c(k)} drives {@link Qsys_mapg1k},
 * whose embedded chain returns the joint law of buffer level and MAP phase. A
 * class-k arrival leaves phase i at rate (D1c(k) e)_i, so the rate of class-k
 * arrivals meeting a full buffer is pKvec D1c(k) e, and
 * <pre>
 *   lambda_k = pi D1c(k) e,   L_k = (pKvec D1c(k) e)/lambda_k.
 * </pre>
 * This is exact: no independence between classes is assumed and no PASTA
 * argument is used, the phase resolution of pKvec doing that work.
 *
 * <p>Assumes a single server and a service law that is iid and independent of
 * class. Per-class service would make the departure rate depend on which class
 * holds the server, which this model does not represent; build the
 * arrival-weighted mixture instead (see
 * {@code jline.solvers.mam.handlers.Mam_svc_mixture}).
 *
 * <p>References:
 * <br>[1] Chydzinski, A. Per-Flow Throughput of a FIFO Buffer. Applied System
 * Innovation 2026, 9, 112.
 */
public final class Qsys_mmapg1k {

    private Qsys_mmapg1k() {
    }

    /** qsys_mmapg1k with the reference defaults tol = 1e-12, nmax = 200000. */
    public static QsysMmapG1kResult qsys_mmapg1k(Matrix D0, List<Matrix> D1c, QsysServiceLaw svc,
                                                 int K) {
        return qsys_mmapg1k(D0, D1c, svc, K, Qsys_mapg1k.DEFAULT_TOL, Qsys_mapg1k.DEFAULT_NMAX);
    }

    /**
     * MMAP[K]/G/1/K with tail drop.
     *
     * @param D0      M x M hidden transition matrix of the arrival MMAP
     * @param D1c     per-class arrival matrices; D0 + sum_k D1c(k) irreducible
     * @param svc     service law, common to every class
     * @param K       buffer size in packets, K &gt;= 1
     * @param tol     uniformization truncation tolerance
     * @param nmaxCap cap on the uniformization order
     */
    public static QsysMmapG1kResult qsys_mmapg1k(Matrix D0, List<Matrix> D1c, QsysServiceLaw svc,
                                                 int K, double tol, int nmaxCap) {
        if (D1c == null || D1c.isEmpty()) {
            line_error("qsys_mmapg1k", "D1C must hold at least one per-class D1 matrix.");
        }
        int R = D1c.size();
        int M = D0.getNumRows();
        Matrix D1 = new Matrix(M, M);
        for (int k = 0; k < R; k++) {
            Matrix Dk = D1c.get(k);
            if (Dk.getNumRows() != M || Dk.getNumCols() != M) {
                line_error("qsys_mmapg1k", "D1C{" + (k + 1) + "} must be " + M + "x" + M + ".");
            }
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    D1.set(i, j, D1.get(i, j) + Dk.get(i, j));
                }
            }
        }

        QsysMapG1kResult r = Qsys_mapg1k.qsys_mapg1k(D0, D1, svc, K, tol, nmaxCap);
        Matrix pit = Map_prob.map_prob(D0, D1);

        Matrix lam = new Matrix(1, R);
        Matrix Lk = new Matrix(1, R);
        Matrix Tk = new Matrix(1, R);
        for (int k = 0; k < R; k++) {
            Matrix Dk = D1c.get(k);
            double lamk = 0.0;
            double blocked = 0.0;
            for (int i = 0; i < M; i++) {
                double rowsum = 0.0;
                for (int j = 0; j < M; j++) {
                    rowsum += Dk.get(i, j);
                }
                lamk += pit.get(i) * rowsum;
                blocked += r.pKvec.get(0, i) * rowsum;
            }
            lam.set(0, k, lamk);
            double lossk = (lamk > 0) ? blocked / lamk : 0.0;
            Lk.set(0, k, lossk);
            Tk.set(0, k, lamk * (1.0 - lossk));
        }

        QsysMmapG1kResult out = new QsysMmapG1kResult();
        out.throughput = Tk;
        out.lossRatio = Lk;
        out.lambda = lam;
        double sumL = 0.0;
        double sumT = 0.0;
        for (int k = 0; k < R; k++) {
            sumL += lam.get(0, k);
            sumT += Tk.get(0, k);
        }
        out.lambdaAggregate = sumL;
        out.throughputAggregate = sumT;
        out.lossAggregate = r.lossProbability;
        out.p0 = r.p0;
        out.pK = r.pK;
        out.pKvec = r.pKvec;
        out.plevel = r.plevel;
        out.meanQueueLength = r.meanQueueLength;
        out.meanServiceTime = r.meanServiceTime;
        out.utilization = r.utilization;
        out.rho = r.rho;
        return out;
    }
}
