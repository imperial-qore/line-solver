/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.handlers;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.function.DoubleUnaryOperator;

import jline.api.mam.Map_count_idc;
import jline.api.mam.Map_idc;
import jline.api.mam.Map_lambda;
import jline.api.mam.Map_pie;
import jline.api.npfqn.Npfqn_feedback_elim;
import jline.api.npfqn.IdcFunction;
import jline.api.npfqn.Npfqn_traffic_idc;
import jline.api.qsys.Qsys_gig1_rq;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Station;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.solvers.mva.SolverMVA;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Robust Queueing Network Analyzer (RQNA) based on indices of dispersion.
 *
 * <p>Approximates the steady-state performance of a single-class open queueing
 * network of single-server FCFS queues with Markovian routing and general
 * (non-renewal) external arrival and (non-exponential) service processes.
 *
 * <p>Reference: W. Whitt and W. You (2018), "A Robust Queueing Network Analyzer
 * Based on Indices of Dispersion", INFORMS J. on Computing. Implements
 * Algorithm 1 (traffic-rate equations, limiting variability equations,
 * time-dependent IDC equations with default alpha/beta corrections, and the
 * robust-queueing workload approximation), plus the near-immediate feedback
 * elimination of Algorithm 2 / Section 4.2 (default on).
 */
public final class Solver_rqna {
    private Solver_rqna() {}

    public static MVAResult solver_rqna(NetworkStruct sn, SolverOptions options) {
        // One predicate for the gate and the run: SolverMVA.supportsModelMethod asks
        // the same question before the report offers "rqna", so the sentence a
        // caller reads here is the sentence that kept the row off the report.
        String rqnaReason = SolverMVA.singleClassOpenReason(sn, "rqna");
        if (!rqnaReason.isEmpty()) {
            throw new RuntimeException(rqnaReason);
        }
        for (int r = 0; r < sn.njobs.length(); r++) {
            if (Double.isFinite(sn.njobs.get(r))) {
                throw new RuntimeException("RQNA supports open networks only (no closed classes).");
            }
        }

        int M = sn.nstations;
        double tol = options.tol;

        Matrix Q = new Matrix(M, 1);
        Matrix U = new Matrix(M, 1);
        Matrix R = new Matrix(M, 1);
        Matrix T = new Matrix(M, 1);
        Matrix AN = new Matrix(M, 1);
        Matrix WN = new Matrix(M, 1);
        Matrix X = new Matrix(1, 1);

        // ----- identify source and queueing stations -----
        boolean[] isSource = new boolean[M];
        boolean[] schedInf = new boolean[M];
        List<Integer> qstatList = new ArrayList<Integer>();
        int src = -1;
        for (int i = 0; i < M; i++) {
            int nd = (int) sn.stationToNode.get(i);
            isSource[i] = (sn.nodetype.get(nd) == NodeType.Source);
            schedInf[i] = (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF);
            if (isSource[i]) {
                if (src < 0) src = i;
            } else {
                qstatList.add(i);
            }
        }
        if (src < 0) {
            throw new RuntimeException("RQNA requires an open network with a Source station.");
        }
        int nq = qstatList.size();
        int[] qstat = new int[nq];
        for (int a = 0; a < nq; a++) qstat[a] = qstatList.get(a);

        // see _kb/06-solver-catalog.md for rationale
        Matrix rtS = sn.rt;

        MatrixCell arvMAP = sn.proc.get(sn.stations.get(src)).get(sn.jobclasses.get(0));
        final Matrix arvD0 = arvMAP.get(0);
        final Matrix arvD1 = arvMAP.get(1);
        double lambda_src = Map_lambda.map_lambda(arvD0, arvD1);
        double c2_src = Map_idc.map_idc(arvD0, arvD1);
        final double arvIdcInf = c2_src;

        // ----- per-queue data -----
        final double[] mu = new double[nq];
        final double[] cs2 = new double[nq];
        final double[] lambda0 = new double[nq];
        final MatrixCell[] svcMAP = new MatrixCell[nq];
        final double[] qsplit = new double[nq];
        Matrix P = new Matrix(nq, nq);
        for (int a = 0; a < nq; a++) {
            int ia = qstat[a];
            mu[a] = sn.rates.get(ia, 0);
            cs2[a] = sn.scv.get(ia, 0);
            svcMAP[a] = sn.proc.get(sn.stations.get(ia)).get(sn.jobclasses.get(0));
            qsplit[a] = rtS.get(src, ia);
            lambda0[a] = lambda_src * qsplit[a];
            for (int b = 0; b < nq; b++) {
                P.set(a, b, rtS.get(ia, qstat[b]));
            }
        }

        // external arrival IDC seen by each queue = split of the source process
        double[] c2a0 = new double[nq];
        for (int a = 0; a < nq; a++) {
            c2a0[a] = qsplit[a] * c2_src + (1.0 - qsplit[a]);
        }
        final MatrixCell arvMAPf = arvMAP;
        final double[] qsplitf = qsplit;
        IdcFunction a0IdcFun = new IdcFunction() {
            @Override
            public double[] eval(double t) {
                double idc = Map_count_idc.map_count_idc(arvMAPf, t);
                double[] out = new double[qsplitf.length];
                for (int a = 0; a < qsplitf.length; a++) {
                    out[a] = qsplitf[a] * idc + (1.0 - qsplitf[a]);
                }
                return out;
            }
        };
        Npfqn_traffic_idc.PerQueueIdc sIdcFun = new Npfqn_traffic_idc.PerQueueIdc() {
            @Override
            public double[] evalVec(double[] perQueueTimes) {
                return svcIdc(svcMAP, perQueueTimes);
            }
        };

        // ----- correction toggles -----
        boolean useAlpha = (options.config == null || options.config.rqna_alpha == null)
                || options.config.rqna_alpha.booleanValue();
        boolean useBeta = (options.config == null || options.config.rqna_beta == null)
                || options.config.rqna_beta.booleanValue();
        boolean doElim = (options.config == null || options.config.rqna_feedback_elim == null)
                || options.config.rqna_feedback_elim.booleanValue();

        final Npfqn_traffic_idc ctx = Npfqn_traffic_idc.npfqn_traffic_idc(
                lambda0, P, c2a0, a0IdcFun, mu, cs2, sIdcFun, useAlpha, useBeta);
        double[] lambda = ctx.lambda;
        double[] rho = ctx.rho;

        // ----- per-queue robust-queueing workload and performance -----
        for (int a = 0; a < nq; a++) {
            int ia = qstat[a];
            T.set(ia, 0, lambda[a]);
            AN.set(ia, 0, lambda[a]);
            if (lambda[a] <= 0) {
                continue;
            }
            if (schedInf[ia]) {
                // infinite-server (delay) station: no waiting
                U.set(ia, 0, lambda[a] / mu[a]);
                Q.set(ia, 0, lambda[a] / mu[a]);
                R.set(ia, 0, 1.0 / mu[a]);
                WN.set(ia, 0, 0.0);
                continue;
            }
            double phat = 0.0;
            if (doElim) {
                phat = phatFeedback(P, rho, a);
            }
            double Ra;
            if (phat > tol) {
                Ra = elimResponse(a, P, rho, lambda, mu, cs2, svcMAP, lambda0,
                        arvMAP, qsplit, arvIdcInf, useAlpha, useBeta);
            } else {
                final int aa = a;
                DoubleUnaryOperator IaFun_a = new DoubleUnaryOperator() {
                    @Override
                    public double applyAsDouble(double x) {
                        return ctx.iaFun(x)[aa];
                    }
                };
                double[] zwqx = Qsys_gig1_rq.qsys_gig1_rq(rho[a], mu[a], cs2[a], IaFun_a);
                Ra = zwqx[1] + 1.0 / mu[a];   // per-visit response = waiting + service
            }
            R.set(ia, 0, Ra);
            U.set(ia, 0, rho[a]);
            Q.set(ia, 0, lambda[a] * Ra);        // number in system (Little, incl. service)
            WN.set(ia, 0, Math.max(0.0, Ra - 1.0 / mu[a]));
        }

        T.set(src, 0, lambda_src);
        AN.set(src, 0, 0.0);
        Matrix CN = R.sumCols();
        X.set(0, 0, lambda_src);
        Q.removeNaN();
        U.removeNaN();
        R.removeNaN();
        CN.removeNaN();

        MVAResult result = new MVAResult();
        result.QN = Q;
        result.UN = U;
        result.RN = R;
        result.TN = T;
        result.AN = AN;
        result.WN = WN;
        result.CN = CN;
        result.XN = X;
        result.iter = 1;
        result.logNormConstAggr = Double.NaN;
        return result;
    }

    // Service IDC vector I_{s,a}(t); perQueueTimes[a] is the evaluation time for
    // queue a (the RQNA evaluates each queue's service IDC at rho_a * t).
    private static double[] svcIdc(MatrixCell[] svcMAP, double[] perQueueTimes) {
        int nq = svcMAP.length;
        double[] out = new double[nq];
        for (int a = 0; a < nq; a++) {
            out[a] = Map_count_idc.map_count_idc(svcMAP[a], perQueueTimes[a]);
        }
        return out;
    }

    // Near-immediate feedback probability at station a (Whitt-You flows paper
    // eq. 3.8/3.9, H={a}): probability of returning to a before visiting any
    // station with strictly higher traffic intensity.
    /**
     * Near-immediate feedback probability at station a (Whitt-You eq. 3.8/3.9,
     * H={a}).
     *
     * <p>Delegated to {@link Npfqn_feedback_elim} so that the solver and the API
     * function cannot drift apart: they answer the same question, and a private
     * copy of the rule here is how the two came to differ on ties in the first
     * place.
     *
     * @param P   the routing matrix over queueing stations
     * @param rho the traffic intensity of each station
     * @param a   the station
     * @return the near-immediate feedback probability
     */
    private static double phatFeedback(Matrix P, double[] rho, int a) {
        double[] phat = (double[]) Npfqn_feedback_elim
                .npfqn_feedback_elim(P, rho, null, null, false).get("feedbackProb");
        return phat[a];
    }

    // Near-immediate feedback elimination at station a (Whitt-You Algorithm 2 /
    // Corollary 4.2): censor the equal/lower-rho cloud to instantaneous
    // switches, retain strictly higher-rho stations as real queues, remove the
    // self-return via immediate-feedback elimination with geometric-sum service,
    // re-solve the IDC equations on the reduced network and apply RQ at a.
    private static double elimResponse(int a, Matrix P, double[] rho, double[] lambda,
            final double[] mu, double[] cs2, MatrixCell[] svcMAP, double[] lambda0,
            final MatrixCell arvMAP, final double[] qsplit, final double arvIdcInf,
            boolean useAlpha, boolean useBeta) {
        int nq = P.getNumRows();
        List<Integer> HcL = new ArrayList<Integer>();
        List<Integer> HiL = new ArrayList<Integer>();
        for (int i = 0; i < nq; i++) {
            if (i == a) continue;
            if (rho[i] <= rho[a] + 1e-9) HcL.add(i);
            else HiL.add(i);
        }
        final int[] R = new int[1 + HiL.size()];
        R[0] = a;
        for (int i = 0; i < HiL.size(); i++) R[i + 1] = HiL.get(i);
        int m = R.length;
        final int[] Hc = new int[HcL.size()];
        for (int i = 0; i < HcL.size(); i++) Hc[i] = HcL.get(i);
        int hc = Hc.length;

        Matrix Pred = new Matrix(m, m);
        final Matrix G = new Matrix(hc, m);   // first-passage Hc-station -> R
        if (hc == 0) {
            for (int x = 0; x < m; x++)
                for (int y = 0; y < m; y++)
                    Pred.set(x, y, P.get(R[x], R[y]));
        } else {
            Matrix Phh = new Matrix(hc, hc);
            for (int x = 0; x < hc; x++)
                for (int y = 0; y < hc; y++)
                    Phh.set(x, y, P.get(Hc[x], Hc[y]));
            Matrix Fhc = Matrix.eye(hc).sub(Phh).inv();
            Matrix PhR = new Matrix(hc, m);
            for (int x = 0; x < hc; x++)
                for (int y = 0; y < m; y++)
                    PhR.set(x, y, P.get(Hc[x], R[y]));
            Matrix Gm = Fhc.mult(PhR);
            for (int x = 0; x < hc; x++)
                for (int y = 0; y < m; y++)
                    G.set(x, y, Gm.get(x, y));
            Matrix PRh = new Matrix(m, hc);
            for (int x = 0; x < m; x++)
                for (int y = 0; y < hc; y++)
                    PRh.set(x, y, P.get(R[x], Hc[y]));
            Matrix censored = PRh.mult(Fhc).mult(PhR);
            for (int x = 0; x < m; x++)
                for (int y = 0; y < m; y++)
                    Pred.set(x, y, P.get(R[x], R[y]) + censored.get(x, y));
        }
        double phat = Pred.get(0, 0);
        phat = Math.min(Math.max(phat, 0.0), 1.0 - 1e-9);
        if (phat > 0) {
            for (int y = 0; y < m; y++) Pred.set(0, y, Pred.get(0, y) / (1.0 - phat));
        }
        Pred.set(0, 0, 0.0);

        // first-passage external arrival rate and IDC into each retained station
        final double[] lam0R = new double[m];
        for (int rr = 0; rr < m; rr++) {
            lam0R[rr] = lambda0[R[rr]];
            for (int ii = 0; ii < hc; ii++) {
                lam0R[rr] += lambda0[Hc[ii]] * G.get(ii, rr);
            }
        }
        double[] c2a0R = new double[m];
        for (int rr = 0; rr < m; rr++) {
            if (lam0R[rr] <= 0) continue;
            double s = lambda0[R[rr]] * (qsplit[R[rr]] * arvIdcInf + (1.0 - qsplit[R[rr]]));
            for (int ii = 0; ii < hc; ii++) {
                double g = G.get(ii, rr);
                double ci = qsplit[Hc[ii]] * arvIdcInf + (1.0 - qsplit[Hc[ii]]);
                s += lambda0[Hc[ii]] * g * (g * ci + (1.0 - g));
            }
            c2a0R[rr] = s / lam0R[rr];
        }
        final MatrixCell arvMAPff = arvMAP;
        final int[] Rf = R;
        final int[] Hcf = Hc;
        final double[] lambda0f = lambda0;
        final double[] qsplitf = qsplit;
        IdcFunction a0IdcR = new IdcFunction() {
            @Override
            public double[] eval(double t) {
                double idc = Map_count_idc.map_count_idc(arvMAPff, t);
                int mm = Rf.length;
                double[] Ir = new double[mm];
                for (int rr = 0; rr < mm; rr++) {
                    if (lam0R[rr] <= 0) { Ir[rr] = 1.0; continue; }
                    double c2R = qsplitf[Rf[rr]] * idc + (1.0 - qsplitf[Rf[rr]]);
                    double num = lambda0f[Rf[rr]] * c2R;
                    for (int ii = 0; ii < Hcf.length; ii++) {
                        double g = G.get(ii, rr);
                        double c2H = qsplitf[Hcf[ii]] * idc + (1.0 - qsplitf[Hcf[ii]]);
                        num += lambda0f[Hcf[ii]] * g * (g * c2H + (1.0 - g));
                    }
                    Ir[rr] = num / lam0R[rr];
                }
                return Ir;
            }
        };

        // service data on R; station a gets the geometric-sum (folded) service
        double[] muR = new double[m];
        double[] cs2R = new double[m];
        final MatrixCell[] svcR = new MatrixCell[m];
        for (int rr = 0; rr < m; rr++) {
            muR[rr] = mu[R[rr]];
            cs2R[rr] = cs2[R[rr]];
            svcR[rr] = svcMAP[R[rr]];
        }
        svcR[0] = geomMap(svcMAP[a], phat);
        muR[0] = (1.0 - phat) * mu[a];
        cs2R[0] = phat + (1.0 - phat) * cs2[a];
        Npfqn_traffic_idc.PerQueueIdc sIdcR = new Npfqn_traffic_idc.PerQueueIdc() {
            @Override
            public double[] evalVec(double[] perQueueTimes) {
                return svcIdc(svcR, perQueueTimes);
            }
        };

        final Npfqn_traffic_idc ctxR = Npfqn_traffic_idc.npfqn_traffic_idc(
                lam0R, Pred, c2a0R, a0IdcR, muR, cs2R, sIdcR, useAlpha, useBeta);
        DoubleUnaryOperator IaFunA = new DoubleUnaryOperator() {
            @Override
            public double applyAsDouble(double x) {
                return ctxR.iaFun(x)[0];
            }
        };
        double[] zwqx = Qsys_gig1_rq.qsys_gig1_rq(rho[a], muR[0], cs2R[0], IaFunA);
        return (1.0 - phat) * zwqx[1] + 1.0 / mu[a];   // per-visit adjustment (visits 1/(1-phat))
    }

    // Geometric random sum of i.i.d. PH service times with success prob (1-p):
    // PH(alpha,T) -> PH(alpha, T + p*t0*alpha), t0 = -T*e. Yields the
    // immediate-feedback-eliminated service (Whitt-You Section 4.1).
    private static MatrixCell geomMap(MatrixCell map, double p) {
        Matrix D0 = map.get(0);
        int n = D0.getNumRows();
        Matrix e = Matrix.ones(n, 1);
        Matrix t0 = D0.mult(e).scale(-1.0);
        Matrix al = Map_pie.map_pie(map.get(0), map.get(1));   // 1 x n
        Matrix D0m = D0.add(p, t0.mult(al));
        Matrix D1m = t0.mult(al).scale(1.0 - p);
        return new MatrixCell(D0m, D1m);
    }
}
