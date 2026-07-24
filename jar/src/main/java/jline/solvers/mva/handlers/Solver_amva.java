/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.handlers;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import jline.GlobalConstants;
import jline.api.pfqn.mva.Pfqn_ab_amva;
import jline.api.pfqn.mva.Pfqn_schmidt_amva;
import jline.api.pfqn.mva.Pfqn_aql;
import jline.api.pfqn.mva.Pfqn_bs;
import jline.api.pfqn.mva.Pfqn_conwayms;
import jline.api.pfqn.mva.Pfqn_linearizermx;
import jline.api.pfqn.mva.Pfqn_sqni;
import jline.api.sn.SnDeaggregateChainResults;
import jline.api.sn.SnGetDemandsChain;
import jline.api.sn.SnGetProductFormChainParams;
import jline.api.sn.SnHasClassSwitching;
import jline.api.sn.SnHasHomogeneousScheduling;
import jline.api.sn.SnHasLoadDependence;
import jline.api.sn.SnHasMultiServer;
import jline.api.sn.SnHasOpenClasses;
import jline.api.sn.SnHasProductForm;
import jline.api.sn.SnHasProductFormNotHetFCFS;
import jline.io.InputOutput;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.util.matrix.Matrix;

/**
 * Handler for the solver_amva function.
 */
public final class Solver_amva {
    private Solver_amva() {}

    public static MVAResult solver_amva(NetworkStruct sn, SolverOptions options) {
        Ret.snGetDemands chainReturn = SnGetDemandsChain.snGetDemandsChain(sn);
        Matrix Lchain = chainReturn.Dchain;
        Matrix STchain = chainReturn.STchain;
        Matrix Vchain = chainReturn.Vchain;
        Matrix alpha = chainReturn.alpha;
        Matrix Nchain = chainReturn.Nchain;

        if (options.config.np_priority == null) {
            options.config.np_priority = "default";
        }
        if (options.config.multiserver == null) {
            options.config.multiserver = "default";
        }
        if (options.config.highvar == null) {
            options.config.highvar = "default";
        }
        String mth = options.method;
        if ("amva.qli".equals(mth)) {
            options.method = "qli";
        } else if ("amva.qd".equals(mth) || "amva.qdamva".equals(mth) || "qdamva".equals(mth)) {
            options.method = "qd";
        } else if ("amva.lin".equals(mth)) {
            options.method = "lin";
        } else if ("amva.gflin".equals(mth)) {
            options.method = "gflin";
        } else if ("amva.egflin".equals(mth)) {
            options.method = "egflin";
        } else if ("amva.qdlin".equals(mth)) {
            options.method = "qdlin";
        } else if ("amva.fli".equals(mth)) {
            options.method = "fli";
        } else if ("amva.bs".equals(mth)) {
            options.method = "bs";
        } else if ("default".equals(mth) || "amva".equals(mth)) {
            double NchainSum = 0.0;
            boolean flag = false;
            int i = 0;
            while (i < Nchain.getNumRows()) {
                int j = 0;
                while (j < Nchain.getNumCols()) {
                    NchainSum += Nchain.get(i, j);
                    if (Nchain.get(i, j) < 1) flag = true;
                    j++;
                }
                i++;
            }
            if (NchainSum <= 2 || flag) {
                options.method = "qd";
            } else {
                options.method = "egflin";
                int ii = 0;
                while (ii < sn.nstations) {
                    if (sn.nservers.get(ii) > 1 && sn.nservers.get(ii) < Integer.MAX_VALUE) {
                        options.method = "lin";
                    }
                    ii++;
                }
            }
        }

        // trivial models
        if (SnHasHomogeneousScheduling.snHasHomogeneousScheduling(sn, SchedStrategy.INF)) {
            options.config.multiserver = "default";
            return Solver_amvald.solver_amvald(sn, options);
        }

        ArrayList<Integer> queueIdx = new ArrayList<Integer>();
        ArrayList<Integer> delayIdx = new ArrayList<Integer>();
        ArrayList<Integer> sourceIdx = new ArrayList<Integer>();
        for (int i = 0; i < sn.nodetype.size(); i++) {
            NodeType nt = sn.nodetype.get(i);
            if (nt == NodeType.Source) sourceIdx.add(i);
            else if (nt == NodeType.Queue) queueIdx.add(i);
            else if (nt == NodeType.Delay) delayIdx.add(i);
            else if (nt == NodeType.Join) queueIdx.add(i);
        }

        // Run amva method
        int M = sn.nstations;
        int C = sn.nchains;
        Matrix V = new Matrix(M, C);
        for (Integer s : sourceIdx) {
            int i = (int) sn.nodeToStation.get(s);
            for (int c = 0; c < sn.nchains; c++) {
                double rateSum = 0.0;
                for (Integer sp : sourceIdx) {
                    for (int r = 0; r < sn.chains.getNumCols(); r++) {
                        if (sn.chains.get(c, r) == 0.0) continue;
                        if (!Double.isNaN(sn.rates.get((int) sn.nodeToStation.get(sp), r))) {
                            rateSum += sn.rates.get((int) sn.nodeToStation.get(sp), r);
                        }
                    }
                }
                if (rateSum > 0) V.set(i, c, 1);
            }
        }
        Matrix Q = new Matrix(M, C);
        Matrix U = new Matrix(M, C);
        Matrix X = null;
        int totiter = 0;

        // see _kb/06-solver-catalog.md for rationale
        boolean cond = SnHasProductFormNotHetFCFS.snHasProductFormNotHetFCFS(sn)
                && !SnHasLoadDependence.snHasLoadDependence(sn)
                && (sn.cdscaling == null || sn.cdscaling.isEmpty())
                && (sn.jdscaling == null || sn.jdscaling.isEmpty())
                && (!SnHasOpenClasses.snHasOpenClasses(sn)
                    || (SnHasProductForm.snHasProductForm(sn) && SnHasOpenClasses.snHasOpenClasses(sn) && "lin".equals(options.method)));
        if (cond) {
            Ret.snGetProductFormParams ret = SnGetProductFormChainParams.snGetProductFormChainParams(sn);
            Matrix L0 = new Matrix(ret.D);
            Matrix L = ret.D;
            Matrix N = ret.N;
            Matrix Z = ret.Z;
            Matrix Z0 = new Matrix(ret.Z);
            Matrix nservers = ret.S;
            Matrix lambda = ret.lambda;

            // see _kb/06-solver-catalog.md for rationale
            Matrix Q0 = null;
            if (options.init_sol != null && !options.init_sol.isEmpty()) {
                Matrix Qi = options.init_sol;
                if (Qi.getNumRows() == sn.nstations && Qi.getNumCols() == C
                        && queueIdx.size() == L.getNumRows()) {
                    Q0 = new Matrix(queueIdx.size(), C);
                    int ri = 0;
                    for (Integer ii : queueIdx) {
                        int st = (int) sn.nodeToStation.get(ii);
                        for (int j = 0; j < C; j++) {
                            Q0.set(ri, j, Qi.get(st, j));
                        }
                        ri++;
                    }
                } else if (Qi.getNumRows() == L.getNumRows() && Qi.getNumCols() == C) {
                    Q0 = Qi;
                }
                if (Q0 != null && Q0.elementMin() < 0) {
                    Q0 = null;
                }
            }

            int vIdx = 0;
            for (int i = 0; i < sn.nodetype.size(); i++) {
                NodeType nt = sn.nodetype.get(i);
                if (nt != NodeType.Queue && nt != NodeType.Delay) continue;
                for (int j = 0; j < V.getNumCols(); j++) {
                    V.set((int) sn.nodeToStation.get(i), j, ret.V.get(vIdx, j));
                }
                vIdx++;
            }
            String ms = options.config.multiserver;
            if ("default".equals(ms) || "seidmann".equals(ms) || "rolia".equals(ms)) {
                Matrix nserversRep = nservers.columnMajorOrder().repmat(1, C);
                int i = 0;
                while (i < L.getNumRows()) {
                    int j = 0;
                    while (j < L.getNumCols()) {
                        L.set(i, j, L.get(i, j) / nserversRep.get(i, j));
                        j++;
                    }
                    i++;
                }
                int j = 0;
                while (j < L.getNumRows()) {
                    int k = 0;
                    while (k < Z.getNumCols()) {
                        Z.set(0, k, Z.get(0, k) + L0.get(j, k) * (nservers.get(j) - 1) / nservers.get(j));
                        k++;
                    }
                    j++;
                }
            } else if ("softmin".equals(ms)) {
                return Solver_amvald.solver_amvald(sn, options);
            }

            String method = options.method;
            if ("sqni".equals(method)) {
                if (sn.nstations == 2) {
                    Pfqn_sqni.PfqnSqniResult result = Pfqn_sqni.pfqn_sqni(N, L, Z);
                    Q = result.Q.copy();
                    U = result.U.copy();
                    X = result.X.copy();
                    totiter = 1;
                } else {
                    InputOutput.line_error(InputOutput.mfilename(new Object() {}),
                            "SQNI cannot handle more than a single queue with a delay.");
                }
                SchedStrategy[] schd = new SchedStrategy[queueIdx.size()];
                int i = 0;
                while (i < queueIdx.size()) {
                    schd[i] = sn.sched.get(sn.stations.get(queueIdx.get(i)));
                    i++;
                }
                Ret.pfqnAMVA bsret = Pfqn_bs.pfqn_bs(L, N, Z, options.tol, options.iter_max, Q0, schd);
                X = bsret.X;
                int iResult = 0;
                for (Integer ii : queueIdx) {
                    int j = 0;
                    while (j < C) {
                        Q.set((int) sn.nodeToStation.get(ii), j, bsret.Q.get(iResult, j));
                        U.set((int) sn.nodeToStation.get(ii), j, bsret.U.get(iResult, j));
                        j++;
                    }
                    iResult++;
                }
                totiter = bsret.totiter;
            } else if ("bs".equals(method)) {
                SchedStrategy[] schd = new SchedStrategy[queueIdx.size()];
                int i = 0;
                while (i < queueIdx.size()) {
                    schd[i] = sn.sched.get(sn.stations.get(queueIdx.get(i)));
                    i++;
                }
                Ret.pfqnAMVA bsret = Pfqn_bs.pfqn_bs(L, N, Z, options.tol, options.iter_max, Q0, schd);
                X = bsret.X;
                int iResult = 0;
                for (Integer ii : queueIdx) {
                    int j = 0;
                    while (j < C) {
                        Q.set((int) sn.nodeToStation.get(ii), j, bsret.Q.get(iResult, j));
                        U.set((int) sn.nodeToStation.get(ii), j, bsret.U.get(iResult, j));
                        j++;
                    }
                    iResult++;
                }
                totiter = bsret.totiter;
            } else if ("aql".equals(method)) {
                if (SnHasMultiServer.snHasMultiServer(sn)) {
                    InputOutput.line_error(InputOutput.mfilename(new Object() {}),
                            "AQL cannot handle multi-server stations. Try with the \"default\" or \"lin\" methods.");
                }
                Ret.pfqnAMVA aqlret = Pfqn_aql.pfqn_aql(L, N, Z, options.tol, options.iter_max, Q0);
                X = aqlret.X;
                int idxResult = 0;
                for (Integer ii : queueIdx) {
                    int j = 0;
                    while (j < C) {
                        Q.set((int) sn.nodeToStation.get(ii), j, aqlret.Q.get(idxResult, j));
                        U.set((int) sn.nodeToStation.get(ii), j, aqlret.U.get(idxResult, j));
                        j++;
                    }
                    idxResult++;
                }
                totiter = aqlret.totiter;
            } else if ("ab".equals(method)) {
                List<SchedStrategy> schedStrategies = new ArrayList<SchedStrategy>();
                for (Station station : sn.stations) {
                    if (station instanceof Queue) {
                        schedStrategies.add(((Queue) station).getSchedStrategy());
                    } else {
                        schedStrategies.add(SchedStrategy.INF);
                    }
                }
                Ret.pfqnAMVAMS abres = Pfqn_ab_amva.ab_amva(STchain, N, V, sn.nservers, schedStrategies, false, "ab");
                X = abres.X;
                int idxResult = 0;
                for (Integer ii : queueIdx) {
                    int j = 0;
                    while (j < C) {
                        Q.set((int) sn.nodeToStation.get(ii), j, abres.Q.get(idxResult, j));
                        U.set((int) sn.nodeToStation.get(ii), j, abres.U.get(idxResult, j));
                        j++;
                    }
                    idxResult++;
                }
            } else if ("schmidt".equals(method)) {
                List<SchedStrategy> schedStrategies = new ArrayList<SchedStrategy>();
                for (Station station : sn.stations) {
                    if (station instanceof Queue) {
                        schedStrategies.add(((Queue) station).getSchedStrategy());
                    } else {
                        schedStrategies.add(SchedStrategy.INF);
                    }
                }
                Ret.pfqnAMVASchmidt schmidtres = Pfqn_schmidt_amva.pfqn_schmidt(sn.rates, sn.njobs, sn.nservers, V, schedStrategies);
                X = schmidtres.X;
                int idxResult = 0;
                for (Integer ii : queueIdx) {
                    int j = 0;
                    while (j < C) {
                        Q.set((int) sn.nodeToStation.get(ii), j, schmidtres.Q.get(idxResult, j));
                        U.set((int) sn.nodeToStation.get(ii), j, schmidtres.U.get(idxResult, j));
                        j++;
                    }
                    idxResult++;
                }
            } else if ("schmidt-ext".equals(method)) {
                List<SchedStrategy> schedStrategies = new ArrayList<SchedStrategy>();
                for (Station station : sn.stations) {
                    if (station instanceof Queue) {
                        schedStrategies.add(((Queue) station).getSchedStrategy());
                    } else {
                        schedStrategies.add(SchedStrategy.INF);
                    }
                }
                Ret.pfqnAMVASchmidt schmidtExtRes = Pfqn_schmidt_amva.pfqn_schmidt_ext(sn.rates, sn.njobs, sn.nservers, V, schedStrategies);
                X = schmidtExtRes.X;
                int idxResult = 0;
                for (Integer ii : queueIdx) {
                    int j = 0;
                    while (j < C) {
                        Q.set((int) sn.nodeToStation.get(ii), j, schmidtExtRes.Q.get(idxResult, j));
                        U.set((int) sn.nodeToStation.get(ii), j, schmidtExtRes.U.get(idxResult, j));
                        j++;
                    }
                    idxResult++;
                }
            } else if ("lin".equals(method) || "gflin".equals(method) || "egflin".equals(method)) {
                if (nservers.elementMax() == 1.0) {
                    SchedStrategy[] schdi = new SchedStrategy[queueIdx.size()];
                    int idx = 0;
                    for (Integer ii : queueIdx) {
                        schdi[idx] = sn.sched.get(sn.stations.get((int) sn.nodeToStation.get(ii)));
                        idx++;
                    }
                    Ret.pfqnAMVA res = Pfqn_linearizermx.pfqn_linearizermx(lambda, L, N, Z, nservers, schdi,
                            options.tol, options.iter_max, options.method, Q0);
                    int iRes = 0;
                    for (Integer ii : queueIdx) {
                        int j = 0;
                        while (j < C) {
                            Q.set((int) sn.nodeToStation.get(ii), j, res.Q.get(iRes, j));
                            U.set((int) sn.nodeToStation.get(ii), j, res.U.get(iRes, j));
                            j++;
                        }
                        iRes++;
                    }
                    X = res.X;
                    totiter = res.totiter;
                } else {
                    String mss = options.config.multiserver;
                    if ("conway".equals(mss)) {
                        SchedStrategy[] schdi = new SchedStrategy[queueIdx.size()];
                        int idx = 0;
                        for (Integer ii : queueIdx) {
                            schdi[idx] = sn.sched.get(sn.stations.get((int) sn.nodeToStation.get(ii)));
                            idx++;
                        }
                        Ret.pfqnAMVAMS res = Pfqn_conwayms.pfqn_conwayms(L, N, Z, nservers.toIntArray1D(),
                                schdi, options.tol, options.iter_max, Q0);
                        int iRes = 0;
                        for (Integer ii : queueIdx) {
                            int j = 0;
                            while (j < C) {
                                Q.set((int) sn.nodeToStation.get(ii), j, res.Q.get(iRes, j));
                                U.set((int) sn.nodeToStation.get(ii), j, res.U.get(iRes, j));
                                j++;
                            }
                            iRes++;
                        }
                        X = res.X;
                        totiter = res.totiter;
                    } else if ("erlang".equals(mss)) {
                        options.config.multiserver = "default";
                        return Solver_amvald.solver_amvald(sn, options);
                    } else if ("krzesinski".equals(mss)) {
                        SchedStrategy[] schdi = new SchedStrategy[queueIdx.size()];
                        int idx = 0;
                        for (Integer ii : queueIdx) {
                            schdi[idx] = sn.sched.get(sn.stations.get((int) sn.nodeToStation.get(ii)));
                            idx++;
                        }
                        Ret.pfqnAMVA res = Pfqn_linearizermx.pfqn_linearizermx(lambda, L, N, Z, nservers, schdi,
                                options.tol, options.iter_max, "default", Q0);
                        int iRes = 0;
                        for (Integer ii : queueIdx) {
                            int j = 0;
                            while (j < C) {
                                Q.set((int) sn.nodeToStation.get(ii), j, res.Q.get(iRes, j));
                                U.set((int) sn.nodeToStation.get(ii), j, res.U.get(iRes, j));
                                j++;
                            }
                            iRes++;
                        }
                        X = res.X;
                        totiter = res.totiter;
                    } else if ("default".equals(mss) || "softmin".equals(mss) || "seidmann".equals(mss) || "rolia".equals(mss) || "suri".equals(mss)) {
                        return Solver_amvald.solver_amvald(sn, options);
                    }
                }
            } else {
                String mss2 = options.config.multiserver;
                if ("conway".equals(mss2) || "erlang".equals(mss2) || "krzesinski".equals(mss2)) {
                    options.config.multiserver = "default";
                }
                return Solver_amvald.solver_amvald(sn, options);
            }

            // Compute performance at delay, then unapply seidmann if needed
            for (int i = 0; i < Z0.getNumRows(); i++) {
                if (!delayIdx.isEmpty()) {
                    Matrix mult = X.repmat(delayIdx.size(), 1).elementMult(Z, null);
                    int zidx = 0;
                    for (Integer d : delayIdx) {
                        for (int j = 0; j < Q.getNumCols(); j++) {
                            Q.set((int) sn.nodeToStation.get(d), j, mult.get(zidx, j));
                            U.set((int) sn.nodeToStation.get(d), j, mult.get(zidx, j));
                        }
                        zidx++;
                    }
                }
                String ms2 = options.config.multiserver;
                if ("default".equals(ms2) || "seidmann".equals(ms2)) {
                    if (!"ab".equals(options.method) && !options.method.startsWith("schmidt")) {
                        int j = 0;
                        while (j < L.getNumRows()) {
                            if (i == 0 && nservers.get(j) > 1) {
                                int k = 0;
                                while (k <= j) {
                                    if (k >= queueIdx.size()) break;
                                    int jq = queueIdx.get(k);
                                    int l = 0;
                                    while (l < Q.getNumCols()) {
                                        Q.set(jq, l,
                                                Q.get(jq, l) + (L0.get(j, l) * (nservers.get(j) - 1) / nservers.get(j)) * X.get(l));
                                        l++;
                                    }
                                    k++;
                                }
                            }
                            j++;
                        }
                    }
                }
            }
            Matrix T = V.elementMult(X.repmat(M, 1), null);
            Matrix R = new Matrix(Q.getNumRows(), Q.getNumCols());
            for (int i = 0; i < R.getNumRows(); i++) {
                for (int j = 0; j < R.getNumCols(); j++) {
                    if (V.get(i, j) <= GlobalConstants.Zero) {
                        R.set(i, j, 0.0);
                    } else {
                        R.set(i, j, Q.get(i, j) / T.get(i, j));
                    }
                }
            }
            Matrix Cm = new Matrix(N.getNumRows(), N.getNumCols());
            for (int i = 0; i < Cm.getNumRows(); i++) {
                for (int j = 0; j < Cm.getNumCols(); j++) {
                    Cm.set(i, j, N.get(i, j) / X.get(i, j) - Z.get(i, j));
                }
            }
            double lG = Double.NaN;
            MVAResult result = new MVAResult();
            result.QN = Q;
            result.UN = U;
            result.RN = R;
            result.TN = T;
            result.CN = Cm;
            result.XN = X;
            result.logNormConstAggr = lG;
            result.iter = totiter;
            result.method = options.method;
            if (SnHasClassSwitching.snHasClassSwitching(sn)) {
                Ret.snDeaggregateChainResults ret1 = SnDeaggregateChainResults.snDeaggregateChainResults(sn, Lchain,
                        null, STchain, Vchain, alpha, null, null, R, T, null, X);
                result.QN = ret1.Q;
                result.UN = ret1.U;
                result.RN = ret1.R;
                result.TN = ret1.T;
                result.CN = ret1.C;
                result.XN = ret1.X;
            }
            return result;
        } else {
            String mss = options.config.multiserver;
            if ("conway".equals(mss) || "erlang".equals(mss) || "krzesinski".equals(mss)) {
                options.config.multiserver = "default";
            }
            return Solver_amvald.solver_amvald(sn, options);
        }
    }
}
