/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.qsys;

import java.util.List;

import jline.api.mam.Map_lambda;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import static jline.io.InputOutput.line_error;

/**
 * Per-flow throughput and loss ratio of a FIFO buffer with tail drop fed by N
 * flows of arbitrary and mutually different statistical character.
 *
 * <p>Port of {@code matlab/src/api/qsys/qsys_mapg1k_perflow.m}. Flow n is
 * described by its own MAP, so two flows may share an arrival rate and still
 * differ in the shape and autocorrelation of their interarrival times.
 *
 * <p>METHOD. The exact model of N flows would need a chain tracking the
 * modulating state of every flow jointly with the buffer occupancy, hence
 * prod_n M_n (K+1) states, already out of reach at N = 10, M_n = 3, K = 10.
 * Instead one model per flow is solved: flow n is kept exactly as MAP_n while
 * the other N-1 flows are replaced by a single Poisson stream of rate
 * lambda - lambda_n, which the Palm-Khinchin limiting theorem justifies for a
 * superposition of many point processes. The substitution is applied N times,
 * once per flow, so no flow is ever the one being Poissonized when its own
 * throughput is computed. Superposing MAP_n with the Poisson background yields
 * <pre>
 *   D0 = D0n - lambdaBar_n I,   D1 = D1n + lambdaBar_n I           ([1], eq. 5)
 * </pre>
 * which is passed to {@link Qsys_mapg1k}. The sweep is O(N (K M)^3) against the
 * O(M^(3N) K^3) of the exact joint model: linear rather than exponential in N.
 *
 * <p>Accuracy. [1] reports errors against simulation of the exact model below
 * about 8 per cent for N &gt;= 9 with K &gt;= 20, falling to 2.1 per cent at K = 50
 * and 0.5 per cent at K = 100. Errors are largest when flows are few, highly
 * variable, and the buffer is small.
 *
 * <p>References:
 * <br>[1] Chydzinski, A. Per-Flow Throughput of a FIFO Buffer. Applied System
 * Innovation 2026, 9, 112, Theorem 1.
 */
public final class Qsys_mapg1k_perflow {

    private Qsys_mapg1k_perflow() {
    }

    /** qsys_mapg1k_perflow with the reference defaults tol = 1e-12, nmax = 200000. */
    public static QsysMapG1kPerflowResult qsys_mapg1k_perflow(List<MatrixCell> maps,
                                                              QsysServiceLaw svc, int K) {
        return qsys_mapg1k_perflow(maps, svc, K, Qsys_mapg1k.DEFAULT_TOL, Qsys_mapg1k.DEFAULT_NMAX);
    }

    /**
     * Per-flow analysis of a MAP-superposition FIFO buffer.
     *
     * @param maps    per-flow MAPs, each a MatrixCell holding {D0n, D1n}
     * @param svc     service law
     * @param K       buffer size in packets, K &gt;= 1
     * @param tol     uniformization truncation tolerance
     * @param nmaxCap cap on the uniformization order
     */
    public static QsysMapG1kPerflowResult qsys_mapg1k_perflow(List<MatrixCell> maps,
                                                              QsysServiceLaw svc, int K,
                                                              double tol, int nmaxCap) {
        if (maps == null || maps.isEmpty()) {
            line_error("qsys_mapg1k_perflow", "at least one flow is required.");
        }
        int N = maps.size();
        double[] lam = new double[N];
        double lamTot = 0.0;
        for (int n = 0; n < N; n++) {
            MatrixCell mn = maps.get(n);
            if (mn == null || mn.size() < 2) {
                line_error("qsys_mapg1k_perflow",
                        "MAPS{" + (n + 1) + "} must hold a {D0,D1} pair.");
            }
            lam[n] = Map_lambda.map_lambda(mn.get(0), mn.get(1));
            lamTot += lam[n];
        }

        Matrix Tn = new Matrix(1, N);
        Matrix Ln = new Matrix(1, N);
        Matrix p0 = new Matrix(1, N);
        Matrix pK = new Matrix(1, N);
        double Smean = Double.NaN;
        for (int n = 0; n < N; n++) {
            Matrix D0n = maps.get(n).get(0);
            Matrix D1n = maps.get(n).get(1);
            double lamBar = lamTot - lam[n];
            int Mn = D0n.getNumRows();
            Matrix D0 = new Matrix(Mn, Mn);
            Matrix D1 = new Matrix(Mn, Mn);
            for (int i = 0; i < Mn; i++) {
                for (int j = 0; j < Mn; j++) {
                    D0.set(i, j, D0n.get(i, j) - (i == j ? lamBar : 0.0));
                    D1.set(i, j, D1n.get(i, j) + (i == j ? lamBar : 0.0));
                }
            }
            QsysMapG1kResult r = Qsys_mapg1k.qsys_mapg1k(D0, D1, svc, K, tol, nmaxCap);
            Smean = r.meanServiceTime;
            p0.set(0, n, r.p0);
            pK.set(0, n, r.pK);
            // Eq. (20): the aggregate departure rate (1-p0)/S less the background
            // throughput lamBar*(1-pK), the background loss being pK by PASTA.
            double t = (1.0 - r.p0) / Smean + r.pK * lamBar - lamBar;
            Tn.set(0, n, t);
            Ln.set(0, n, 1.0 - t / lam[n]);
        }

        QsysMapG1kPerflowResult out = new QsysMapG1kPerflowResult();
        out.throughput = Tn;
        out.lossRatio = Ln;
        Matrix lamM = new Matrix(1, N);
        double sumT = 0.0;
        double weighted = 0.0;
        for (int n = 0; n < N; n++) {
            lamM.set(0, n, lam[n]);
            sumT += Tn.get(0, n);
            weighted += Ln.get(0, n) * lam[n];
        }
        out.lambda = lamM;
        out.lambdaAggregate = lamTot;
        out.throughputAggregate = sumT;
        out.lossAggregate = weighted / lamTot;
        out.p0 = p0;
        out.pK = pK;
        out.meanServiceTime = Smean;
        out.rho = lamTot * Smean;
        return out;
    }
}
