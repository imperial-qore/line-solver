/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.handlers;

import java.util.ArrayList;
import java.util.List;

import jline.api.mam.Map_lambda;
import jline.api.mam.Map_scv;
import jline.api.npfqn.Npfqn_traffic_rqt;
import jline.api.qsys.Qsys_gigk_rqt;
import jline.api.qsys.Qsys_gigk_rqt_gamma;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.solvers.mva.SolverMVA;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import org.apache.commons.math3.util.FastMath;

/**
 * Robust Queueing Network Analyzer (RQNA) of Robust Queueing Theory.
 *
 * <p>Estimates the steady-state performance of a single-class open network of
 * FCFS queues with Markovian routing by replacing the stochastic primitives with
 * polyhedral uncertainty sets and taking a worst-case view of each node in
 * isolation. The algorithm is Section 7.2 of the reference: the external streams
 * get Gamma_a = sigma_a, the effective arrival process at each node follows from
 * the network characterization of Theorem 10, the service variability parameter
 * from the adaptation of Section 7.1, and the system time at each node is the
 * worst-case bound of Theorem 3. The published step 3, path enumeration, is not
 * needed here: LINE aggregates per-node system times into per-class response
 * times through the visit ratios.
 *
 * <p>The adaptation is regressed against simulation in heavy traffic, so accuracy
 * degrades at low utilization: on M/M/1 the error is about 5% at rho=0.9 but over
 * 50% at rho=0.5.
 *
 * <p>Reference: C. Bandi, D. Bertsimas, N. Youssef (2015), "Robust Queueing
 * Theory", Operations Research 63(3), 676-700.
 */
public final class Solver_rqt {
    private Solver_rqt() {}

    public static MVAResult solver_rqt(NetworkStruct sn, SolverOptions options) {
        // One predicate for the gate and the run: SolverMVA.supportsModelMethod asks
        // the same question before the report offers "rqt", so the sentence a caller
        // reads here is the sentence that kept the row off the report.
        String rqtReason = SolverMVA.singleClassOpenReason(sn, "rqt");
        if (!rqtReason.isEmpty()) {
            throw new RuntimeException(rqtReason);
        }
        for (int r = 0; r < sn.njobs.length(); r++) {
            if (Double.isFinite(sn.njobs.get(r))) {
                throw new RuntimeException("RQT supports open networks only (no closed classes).");
            }
        }

        int M = sn.nstations;
        Matrix Q = new Matrix(M, 1);
        Matrix U = new Matrix(M, 1);
        Matrix R = new Matrix(M, 1);
        Matrix T = new Matrix(M, 1);
        Matrix AN = new Matrix(M, 1);
        Matrix WN = new Matrix(M, 1);
        Matrix X = new Matrix(1, 1);

        // ----- configuration -----
        String regime = "independent";
        boolean useExact = false;
        double alphaAcfg = 2.0;
        double alphaScfg = 2.0;
        if (options.config != null) {
            if (options.config.rqt_regime != null) {
                regime = options.config.rqt_regime;
            }
            if (options.config.rqt_exact != null) {
                useExact = options.config.rqt_exact.booleanValue();
            }
            if (options.config.rqt_alpha_a != null) {
                alphaAcfg = options.config.rqt_alpha_a.doubleValue();
            }
            if (options.config.rqt_alpha_s != null) {
                alphaScfg = options.config.rqt_alpha_s.doubleValue();
            }
        }

        // ----- identify source and queueing stations -----
        boolean[] schedInf = new boolean[M];
        List<Integer> qstatList = new ArrayList<Integer>();
        int src = -1;
        for (int i = 0; i < M; i++) {
            int nd = (int) sn.stationToNode.get(i);
            boolean isSource = (sn.nodetype.get(nd) == NodeType.Source);
            schedInf[i] = (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF);
            if (isSource) {
                if (src < 0) {
                    src = i;
                }
            } else {
                qstatList.add(i);
            }
        }
        if (src < 0) {
            throw new RuntimeException("RQT requires an open network with a Source station.");
        }
        int nq = qstatList.size();
        int[] qstat = new int[nq];
        for (int a = 0; a < nq; a++) {
            qstat[a] = qstatList.get(a);
        }

        Matrix rtS = sn.rt;

        // ----- external arrival process -----
        MatrixCell arvMAP = sn.proc.get(sn.stations.get(src)).get(sn.jobclasses.get(0));
        double lambdaSrc = Map_lambda.map_lambda(arvMAP.get(0), arvMAP.get(1));
        double sigmaASrc = FastMath.sqrt(Map_scv.map_scv(arvMAP.get(0), arvMAP.get(1))) / lambdaSrc;

        // ----- per-node primitives -----
        double[] mu = new double[nq];
        double[] sigmaS = new double[nq];
        int[] nserv = new int[nq];
        double[] alphaS = new double[nq];
        double[] lambda0 = new double[nq];
        double[] Gamma0 = new double[nq];
        double[] alpha0 = new double[nq];
        Matrix F = new Matrix(nq, nq);
        for (int a = 0; a < nq; a++) {
            int ia = qstat[a];
            mu[a] = sn.rates.get(ia, 0);
            sigmaS[a] = FastMath.sqrt(sn.scv.get(ia, 0)) / mu[a];
            double ns = sn.nservers.get(ia, 0);
            nserv[a] = (Double.isFinite(ns) && ns > 0) ? (int) ns : 1;
            alphaS[a] = alphaScfg;
            alpha0[a] = alphaAcfg;
            // the source stream reaches node a thinned by q, Theorem 6
            double q = rtS.get(src, ia);
            lambda0[a] = lambdaSrc * q;
            if (q > 0) {
                Gamma0[a] = sigmaASrc * FastMath.pow(1.0 / q, 1.0 / alphaAcfg);
            }
            for (int b = 0; b < nq; b++) {
                F.set(a, b, rtS.get(ia, qstat[b]));
            }
        }

        // ----- effective arrival processes, Theorem 10 -----
        double[][] eff = Npfqn_traffic_rqt.npfqn_traffic_rqt(lambda0, Gamma0, alpha0, F);
        double[] lambda = eff[0];
        double[] gammaA = eff[1];
        double[] alphaA = eff[2];

        // ----- per-node worst-case analysis -----
        for (int a = 0; a < nq; a++) {
            int ia = qstat[a];
            T.set(ia, 0, lambda[a]);
            AN.set(ia, 0, lambda[a]);
            if (lambda[a] <= 0) {
                continue;
            }
            if (schedInf[ia]) {
                U.set(ia, 0, lambda[a] / mu[a]);
                Q.set(ia, 0, lambda[a] / mu[a]);
                R.set(ia, 0, 1.0 / mu[a]);
                WN.set(ia, 0, 0.0);
                continue;
            }
            double rho = lambda[a] / (nserv[a] * mu[a]);
            double gammaS = Qsys_gigk_rqt_gamma.qsys_gigk_rqt_gamma(rho, mu[a], gammaA[a], sigmaS[a],
                    nserv[a], alphaA[a], regime);
            double[] wrs = Qsys_gigk_rqt.qsys_gigk_rqt(lambda[a], mu[a], gammaA[a], gammaS,
                    nserv[a], alphaA[a], alphaS[a]);
            double Ra = useExact ? wrs[2] : wrs[0];
            R.set(ia, 0, Ra);
            U.set(ia, 0, rho);
            Q.set(ia, 0, lambda[a] * Ra);   // Little's law, number in system
            WN.set(ia, 0, FastMath.max(0.0, Ra - 1.0 / mu[a]));
        }

        T.set(src, 0, lambdaSrc);
        AN.set(src, 0, 0.0);
        Matrix CN = R.sumCols();
        X.set(0, 0, lambdaSrc);
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
}
