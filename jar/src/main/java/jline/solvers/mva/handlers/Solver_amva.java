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
import jline.api.pfqn.mva.Pfqn_qsa;
import jline.api.pfqn.mva.Pfqn_bs;
import jline.api.pfqn.mva.Pfqn_lcp;
import jline.api.pfqn.mva.Pfqn_chow;
import jline.api.pfqn.mva.Pfqn_pam;
import jline.api.pfqn.mva.Pfqn_clust;
import jline.api.pfqn.mva.Pfqn_dmlin;
import jline.api.pfqn.mva.Pfqn_conwayms;
import jline.api.pfqn.mva.Pfqn_linearizermx;
import jline.api.pfqn.mva.Pfqn_sqni;
import jline.api.pfqn.mva.Pfqn_scat;
import jline.api.pfqn.mva.Pfqn_tay;
import jline.api.sn.SnDeaggregateChainResults;
import jline.api.sn.SnGetDemandsChain;
import jline.api.sn.SnGetProductFormChainParams;
import jline.api.sn.SnHasClassSwitching;
import jline.api.sn.SnHasHomogeneousScheduling;
import jline.api.sn.SnInterlockChain;
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
import jline.solvers.mva.SolverMVA;
import jline.util.matrix.Matrix;

/**
 * Handler for the solver_amva function.
 */
public final class Solver_amva {
    private Solver_amva() {}

    /**
     * Conservative form of the branch test below: true when this model may be solved by the
     * product-form AMVA kernels (the linearizer family and relatives) rather than by
     * {@link Solver_amvald}. The mixed case is reported true for every resolved method, while
     * the branch itself takes it only for "lin", so a caller reading this is never told that
     * Solver_amvald will run when it might not.
     *
     * <p>Read by {@link jline.solvers.mva.analyzers.Solver_mva_analyzer#mvaCarriesInterlock},
     * which needs to know whether a supplied interlock matrix would be honoured in place or
     * would force the model onto another algorithm.</p>
     */
    public static boolean amvaUsesProductFormKernels(NetworkStruct sn) {
        return SnHasProductFormNotHetFCFS.snHasProductFormNotHetFCFS(sn)
                && !SnHasLoadDependence.snHasLoadDependence(sn)
                && (sn.cdscaling == null || sn.cdscaling.isEmpty())
                && (sn.jdscaling == null || sn.jdscaling.isEmpty())
                && (!SnHasOpenClasses.snHasOpenClasses(sn) || SnHasProductForm.snHasProductForm(sn));
    }

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
        } else if ("amva.qsa".equals(mth)) {
            options.method = "qsa";
        } else if ("amva.aql".equals(mth)) {
            options.method = "aql";
        } else if ("amva.tay".equals(mth)) {
            options.method = "tay";
        } else if ("amva.scat".equals(mth)) {
            options.method = "scat";
        } else if ("amva.ab".equals(mth)) {
            options.method = "ab";
        } else if ("amva.schmidt".equals(mth)) {
            options.method = "schmidt";
        } else if ("amva.schmidt-ext".equals(mth)) {
            options.method = "schmidt-ext";
        } else if ("amva.lcp".equals(mth)) {
            options.method = "lcp";
        } else if ("amva.chow".equals(mth)) {
            options.method = "chow";
        } else if ("amva.pamb".equals(mth)) {
            options.method = "pamb";
        } else if ("amva.pami".equals(mth)) {
            options.method = "pami";
        } else if ("amva.pamt".equals(mth)) {
            options.method = "pamt";
        } else if ("amva.clust".equals(mth)) {
            options.method = "clust";
        } else if ("amva.dmlin".equals(mth)) {
            options.method = "dmlin";
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

        // The closed-population AMVA family (Bard-Schweitzer, SQNI, Tay, SCAT, AQL,
        // QSA, Bard LCP, Chow SA, Hsieh-Lam PAM, clustering, Improved Linearizer,
        // Akyildiz-Bolch, Schmidt) lives ONLY in the product-form branch below. The
        // same predicate the report gates on decides here, so a name the report
        // offers is a name that runs and a name it withholds errors rather than
        // falling through to Solver_amvald and returning the qd-family answer -- or
        // a table of zeros -- under a method the caller did not ask for.
        String amvaReason = SolverMVA.closedPopulationReason(sn, options.method);
        if (!amvaReason.isEmpty()) {
            InputOutput.line_error(InputOutput.mfilename(new Object() {}), amvaReason);
        }

        // trivial models
        if (SnHasHomogeneousScheduling.snHasHomogeneousScheduling(sn, SchedStrategy.INF)) {
            options.config.multiserver = "default";
            return Solver_amvald.solver_amvald(sn, options);
        }

        // Interlocked flow (Franks 1999, Eq. 4.7). options.config.interlock arrives
        // CLASS-indexed and is translated to the chain basis the handlers work in without
        // overwriting it: this method re-enters itself on the Conway fallback below, where the
        // class-level matrix has to survive. Only Solver_amvald carries the correction, so an
        // interlocked model goes there rather than to the closed-form or linearizermx branches,
        // which have no interlock term and would drop it silently.
        options.config.interlock_chain = null;
        if (options.config.interlock != null && !options.config.interlock.isEmpty()) {
            options.config.interlock_chain = SnInterlockChain.snInterlockChain(sn, options.config.interlock);
            if (options.config.interlock_chain != null) {
                return Solver_amvald.solver_amvald(sn, options);
            }
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

        // see _kb/06-solver-catalog.md for rationale. ab / schmidt / schmidt-ext ARE the
        // class-dependent FCFS algorithms and live only in the product-form branch below,
        // so the het-FCFS exclusion must not divert them.
        boolean hetFcfsOwn = "ab".equals(options.method) || "schmidt".equals(options.method)
                || "schmidt-ext".equals(options.method);
        boolean cond = (SnHasProductFormNotHetFCFS.snHasProductFormNotHetFCFS(sn)
                    || (hetFcfsOwn && SnHasProductFormNotHetFCFS.snHasProductFormNotHetFCFS(sn, false)))
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
            // 'rolia' is NOT in this set in MATLAB solver_amva.m:111 ({'default','seidmann'}).
            // Including it applied Seidmann's transform under a rule that does not ask for
            // it, i.e. solved a different model whenever a caller selected rolia.
            if ("default".equals(ms) || "seidmann".equals(ms)) {
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
                // pfqn_sqni is a closed form for one queueing station with a delay.
                int nInf = 0;
                for (Station st : sn.stations) {
                    if (sn.sched.get(st) == SchedStrategy.INF) nInf++;
                }
                if (sn.nstations != 2 || nInf != 1) {
                    InputOutput.line_error(InputOutput.mfilename(new Object() {}),
                            "SQNI is defined for a single queueing station with a delay. Try with the \"default\" or \"lin\" methods.");
                }
                Pfqn_sqni.PfqnSqniResult result = Pfqn_sqni.pfqn_sqni(N, L, Z);
                int iResult = 0;
                for (Integer ii : queueIdx) {
                    int j = 0;
                    while (j < C) {
                        Q.set((int) sn.nodeToStation.get(ii), j, result.Q.get(iResult, j));
                        U.set((int) sn.nodeToStation.get(ii), j, result.U.get(iResult, j));
                        j++;
                    }
                    iResult++;
                }
                X = result.X.copy();
                totiter = 1;
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
            } else if ("lcp".equals(method)) {
                // Bard LCP: the Schweitzer proportional term set to zero
                SchedStrategy[] schd = new SchedStrategy[queueIdx.size()];
                int i = 0;
                while (i < queueIdx.size()) {
                    schd[i] = sn.sched.get(sn.stations.get(queueIdx.get(i)));
                    i++;
                }
                Ret.pfqnAMVA lcpret = Pfqn_lcp.pfqn_lcp(L, N, Z, options.tol, options.iter_max, Q0, schd);
                X = lcpret.X;
                int iResult = 0;
                for (Integer ii : queueIdx) {
                    int j = 0;
                    while (j < C) {
                        Q.set((int) sn.nodeToStation.get(ii), j, lcpret.Q.get(iResult, j));
                        U.set((int) sn.nodeToStation.get(ii), j, lcpret.U.get(iResult, j));
                        j++;
                    }
                    iResult++;
                }
                totiter = lcpret.totiter;
            } else if ("chow".equals(method)) {
                // Chow Second Approximation: theta-terms taken off the LCP solution
                SchedStrategy[] schd = new SchedStrategy[queueIdx.size()];
                int i = 0;
                while (i < queueIdx.size()) {
                    schd[i] = sn.sched.get(sn.stations.get(queueIdx.get(i)));
                    i++;
                }
                Ret.pfqnAMVA chret = Pfqn_chow.pfqn_chow(L, N, Z, options.tol, options.iter_max, Q0, schd);
                X = chret.X;
                int iResult = 0;
                for (Integer ii : queueIdx) {
                    int j = 0;
                    while (j < C) {
                        Q.set((int) sn.nodeToStation.get(ii), j, chret.Q.get(iResult, j));
                        U.set((int) sn.nodeToStation.get(ii), j, chret.U.get(iResult, j));
                        j++;
                    }
                    iResult++;
                }
                totiter = chret.totiter;
            } else if ("pamb".equals(method) || "pami".equals(method) || "pamt".equals(method)) {
                // Hsieh-Lam proportional approximations, noniterative
                Ret.pfqnAMVA pamret = Pfqn_pam.pfqn_pam(L, N, Z, method);
                X = pamret.X;
                int iResult = 0;
                for (Integer ii : queueIdx) {
                    int j = 0;
                    while (j < C) {
                        Q.set((int) sn.nodeToStation.get(ii), j, pamret.Q.get(iResult, j));
                        U.set((int) sn.nodeToStation.get(ii), j, pamret.U.get(iResult, j));
                        j++;
                    }
                    iResult++;
                }
                totiter = 1;
            } else if ("clust".equals(method)) {
                // de Souza e Silva-Lavenberg-Muntz clustering approximation
                Ret.pfqnAMVA clret = Pfqn_clust.pfqn_clust(L, N, Z, options.tol, options.iter_max);
                X = clret.X;
                int iResult = 0;
                for (Integer ii : queueIdx) {
                    int j = 0;
                    while (j < C) {
                        Q.set((int) sn.nodeToStation.get(ii), j, clret.Q.get(iResult, j));
                        U.set((int) sn.nodeToStation.get(ii), j, clret.U.get(iResult, j));
                        j++;
                    }
                    iResult++;
                }
                totiter = clret.totiter;
            } else if ("dmlin".equals(method)) {
                // de Souza e Silva-Muntz Improved Linearizer
                Ret.pfqnAMVA dmret = Pfqn_dmlin.pfqn_dmlin(L, N, Z, options.tol, options.iter_max, Q0);
                X = dmret.X;
                int iResult = 0;
                for (Integer ii : queueIdx) {
                    int j = 0;
                    while (j < C) {
                        Q.set((int) sn.nodeToStation.get(ii), j, dmret.Q.get(iResult, j));
                        U.set((int) sn.nodeToStation.get(ii), j, dmret.U.get(iResult, j));
                        j++;
                    }
                    iResult++;
                }
                totiter = dmret.totiter;
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
            } else if ("qsa".equals(method)) {
                if (SnHasMultiServer.snHasMultiServer(sn)) {
                    InputOutput.line_error(InputOutput.mfilename(new Object() {}),
                            "QSA cannot handle multi-server stations. Try with the \"default\" or \"lin\" methods.");
                }
                SchedStrategy[] schdq = new SchedStrategy[queueIdx.size()];
                int iq = 0;
                while (iq < queueIdx.size()) {
                    schdq[iq] = sn.sched.get(sn.stations.get(queueIdx.get(iq)));
                    iq++;
                }
                Ret.pfqnAMVA qsaret = Pfqn_qsa.pfqn_qsa(L, N, Z, schdq, options.tol, options.iter_max, 3, Q0);
                X = qsaret.X;
                int idxQsa = 0;
                for (Integer ii : queueIdx) {
                    int j = 0;
                    while (j < C) {
                        Q.set((int) sn.nodeToStation.get(ii), j, qsaret.Q.get(idxQsa, j));
                        U.set((int) sn.nodeToStation.get(ii), j, qsaret.U.get(idxQsa, j));
                        j++;
                    }
                    idxQsa++;
                }
                totiter = qsaret.totiter;
            } else if ("tay".equals(method)) {
                if (SnHasMultiServer.snHasMultiServer(sn)) {
                    InputOutput.line_error(InputOutput.mfilename(new Object() {}),
                            "Tay's approximation is defined for single-server stations. Try with the \"default\" or \"lin\" methods.");
                }
                Pfqn_tay.Result tayret = Pfqn_tay.pfqn_tay(L, N, Z, options.tol, options.iter_max, Q0);
                X = tayret.X;
                int idxResult = 0;
                for (Integer ii : queueIdx) {
                    int j = 0;
                    while (j < C) {
                        Q.set((int) sn.nodeToStation.get(ii), j, tayret.Q.get(idxResult, j));
                        U.set((int) sn.nodeToStation.get(ii), j, tayret.U.get(idxResult, j));
                        j++;
                    }
                    idxResult++;
                }
                totiter = tayret.totiter;
            } else if ("scat".equals(method)) {
                // Neuse-Chandy SCAT: the Linearizer fixed point with a single Delta
                // refresh. Multiserver stations arrive here already Seidmann-scaled,
                // as they do for bs, so no separate guard is needed.
                SchedStrategy[] schdsc = new SchedStrategy[queueIdx.size()];
                int isc = 0;
                while (isc < queueIdx.size()) {
                    schdsc[isc] = sn.sched.get(sn.stations.get(queueIdx.get(isc)));
                    isc++;
                }
                Ret.pfqnAMVA scatret = Pfqn_scat.pfqn_scat(L, N, Z, schdsc, options.tol, options.iter_max, Q0);
                X = scatret.X;
                int idxScat = 0;
                for (Integer ii : queueIdx) {
                    int j = 0;
                    while (j < C) {
                        Q.set((int) sn.nodeToStation.get(ii), j, scatret.Q.get(idxScat, j));
                        U.set((int) sn.nodeToStation.get(ii), j, scatret.U.get(idxScat, j));
                        j++;
                    }
                    idxScat++;
                }
                totiter = scatret.totiter;
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
                // ONE ROW PER STATION, AND THE DELAYS COME FROM THE KERNEL.
                // STchain, V, sn.nservers and schedStrategies are all indexed by
                // STATION, delays included, so the result is too -- exactly as
                // solver_amva.m maps D_full = [Z0; L0] back with its nDelays
                // offset. Walking abres.Q from row 0 for the QUEUES alone put the
                // first delay's queue length on the first queueing station, and
                // left the delay itself to the X*Z0 rule below; on the repairmen
                // model that made sum(Q) 3.1857 against N = 3. A delay's
                // utilization is its queue length, which is what the reference
                // writes there.
                for (Integer ii : queueIdx) {
                    int ist = (int) sn.nodeToStation.get(ii);
                    for (int j = 0; j < C; j++) {
                        Q.set(ist, j, abres.Q.get(ist, j));
                        U.set(ist, j, abres.U.get(ist, j));
                    }
                }
                for (Integer dd : delayIdx) {
                    int ist = (int) sn.nodeToStation.get(dd);
                    for (int j = 0; j < C; j++) {
                        Q.set(ist, j, abres.Q.get(ist, j));
                        U.set(ist, j, abres.Q.get(ist, j));
                    }
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
                Matrix schmidtRates = schmidtChainRates(sn, STchain);
                Ret.pfqnAMVASchmidt schmidtres = (schmidtRates == null)
                        ? Pfqn_schmidt_amva.pfqn_schmidt(sn.rates, sn.njobs, sn.nservers, V, schedStrategies)
                        : Pfqn_schmidt_amva.pfqn_schmidt(schmidtRates, N, sn.nservers, V, schedStrategies);
                X = schmidtres.X;
                // ONE ROW PER STATION, AND THE DELAYS COME FROM THE KERNEL; see
                // the 'ab' arm above, which had the same off-by-nDelays mapping
                // and the same consequence for sum(Q).
                for (Integer ii : queueIdx) {
                    int ist = (int) sn.nodeToStation.get(ii);
                    for (int j = 0; j < C; j++) {
                        Q.set(ist, j, schmidtres.Q.get(ist, j));
                        U.set(ist, j, schmidtres.U.get(ist, j));
                    }
                }
                for (Integer dd : delayIdx) {
                    int ist = (int) sn.nodeToStation.get(dd);
                    for (int j = 0; j < C; j++) {
                        Q.set(ist, j, schmidtres.Q.get(ist, j));
                        U.set(ist, j, schmidtres.Q.get(ist, j));
                    }
                }
            } else if ("schmidt-ext".equals(method)) {
                // One predicate for the gate and the run, asked about the numbers
                // THIS arm passes -- the CHAIN populations under class switching and
                // the class ones otherwise, which is what Pfqn_schmidt_ext is handed
                // below. It forms its alpha correction from the network with one
                // class-r customer tagged, and an empty one has none to tag.
                List<Boolean> sxFcfs = new ArrayList<Boolean>();
                for (int ist = 0; ist < sn.nstations; ist++) {
                    sxFcfs.add(sn.sched.get(sn.stations.get(ist)) == SchedStrategy.FCFS);
                }
                String sxReason = SolverMVA.schmidtExtReason(
                        SnHasClassSwitching.snHasClassSwitching(sn) ? N : sn.njobs,
                        sxFcfs, "schmidt-ext");
                if (!sxReason.isEmpty()) {
                    InputOutput.line_error(InputOutput.mfilename(new Object() {}), sxReason);
                }
                List<SchedStrategy> schedStrategies = new ArrayList<SchedStrategy>();
                for (Station station : sn.stations) {
                    if (station instanceof Queue) {
                        schedStrategies.add(((Queue) station).getSchedStrategy());
                    } else {
                        schedStrategies.add(SchedStrategy.INF);
                    }
                }
                Matrix schmidtExtRates = schmidtChainRates(sn, STchain);
                Ret.pfqnAMVASchmidt schmidtExtRes = (schmidtExtRates == null)
                        ? Pfqn_schmidt_amva.pfqn_schmidt_ext(sn.rates, sn.njobs, sn.nservers, V, schedStrategies)
                        : Pfqn_schmidt_amva.pfqn_schmidt_ext(schmidtExtRates, N, sn.nservers, V, schedStrategies);
                X = schmidtExtRes.X;
                // ONE ROW PER STATION, AND THE DELAYS COME FROM THE KERNEL; see
                // the 'ab' arm above.
                for (Integer ii : queueIdx) {
                    int ist = (int) sn.nodeToStation.get(ii);
                    for (int j = 0; j < C; j++) {
                        Q.set(ist, j, schmidtExtRes.Q.get(ist, j));
                        U.set(ist, j, schmidtExtRes.U.get(ist, j));
                    }
                }
                for (Integer dd : delayIdx) {
                    int ist = (int) sn.nodeToStation.get(dd);
                    for (int j = 0; j < C; j++) {
                        Q.set(ist, j, schmidtExtRes.Q.get(ist, j));
                        U.set(ist, j, schmidtExtRes.Q.get(ist, j));
                    }
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
                        // options.method, NOT the literal "default": the method name selects
                        // the linearizer variant inside pfqn_linearizermx, so hardcoding it
                        // ran a different algorithm than the caller asked for. Mirrors
                        // MATLAB solver_amva.m:251.
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

            // Compute performance at delay, then unapply seidmann if needed.
            // The delay is charged the ORIGINAL think time Z0: Seidmann folds
            // L(m-1)/m of each multiserver station into Z, but that population is
            // in service at the station and is given back to it below, so charging
            // Z here as well counted it twice and sum(Q) exceeded N.
            // 'ab' is exempt: its kernel was handed the delay rows and answered
            // them, and the arm above has already written them. Recomputing them
            // as X*Z0 would replace a consistent pair (X, Q) with one that does
            // not sum to N.
            if (!delayIdx.isEmpty() && !"ab".equals(options.method)
                    && !options.method.startsWith("schmidt")) {
                Matrix mult = X.repmat(delayIdx.size(), 1).elementMult(Z0, null);
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
                // ab and schmidt were handed the original demands, so the
                // transform was never applied to them
                if (!"ab".equals(options.method) && !options.method.startsWith("schmidt")) {
                    int j = 0;
                    while (j < L.getNumRows()) {
                        // row j belongs to station nodeToStation(queueIdx(j)) alone;
                        // the old loop walked k = 0..j over RAW node indices, so it
                        // both sprayed the term over the earlier queues and indexed
                        // Q in node space instead of station space
                        if (nservers.get(j) > 1 && j < queueIdx.size()) {
                            int jq = (int) sn.nodeToStation.get(queueIdx.get(j));
                            int l = 0;
                            while (l < Q.getNumCols()) {
                                Q.set(jq, l,
                                        Q.get(jq, l) + (L0.get(j, l) * (nservers.get(j) - 1) / nservers.get(j)) * X.get(l));
                                l++;
                            }
                        }
                        j++;
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
            // Cycle time excludes the think time actually spent at the delays,
            // which is Z0 summed over them; with the Seidmann Z it disagreed with
            // the station residence times sum(R.*V) that Q now reports
            Matrix Cm = new Matrix(N.getNumRows(), N.getNumCols());
            for (int i = 0; i < Cm.getNumRows(); i++) {
                for (int j = 0; j < Cm.getNumCols(); j++) {
                    double z0tot = 0.0;
                    for (int zr = 0; zr < Z0.getNumRows(); zr++) {
                        z0tot += Z0.get(zr, j);
                    }
                    Cm.set(i, j, N.get(i, j) / X.get(i, j) - z0tot);
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
            // Nothing of the closed-population family can reach here: it is refused
            // above by SolverMVA.closedPopulationReason, the one predicate the
            // report also gates on. Keeping a second copy of that rule here is what
            // let the two drift, so that "bs", "sqni", "ab" and the two Schmidt arms
            // fell through to Solver_amvald and returned the qd-family answer under
            // their name while "aql", "qsa", "tay" and the Chapter-2 survey
            // algorithms errored.
            String mss = options.config.multiserver;
            if ("conway".equals(mss) || "erlang".equals(mss) || "krzesinski".equals(mss)) {
                options.config.multiserver = "default";
            }
            return Solver_amvald.solver_amvald(sn, options);
        }
    }

    /**
     * The chain-level "rates" matrix Pfqn_schmidt_amva wants, or null when the
     * model needs none.
     *
     * THE CLOSED-POPULATION AMVA FAMILY RECURS ON A CONSERVED POPULATION. Under
     * class switching a job CHANGES CLASS as it moves, so no per-class population
     * is conserved and the vector these kernels need is the CHAIN one -- which is
     * why the whole product-form branch above is built out of
     * SnGetProductFormChainParams and deaggregated at the end. The Schmidt arms
     * were the two that reached past it for sn.rates and sn.njobs, so on a
     * class-switching model they solved a different network and their per-class
     * answer was then read column by column as if it were per-chain and
     * deaggregated a second time.
     *
     * <p>The kernel's entry point takes SERVICE RATES and inverts them internally,
     * so the chain service times have to go in as reciprocals. That round trip is
     * not bit-exact, which is why it is taken ONLY under class switching: with one
     * class per chain the two bases carry the same numbers and sn.rates is the one
     * without a division in it.</p>
     *
     * @param sn      the network struct
     * @param STchain chain-level mean service times (M x C)
     * @return the reciprocal of STchain, or null when the model has no class switching
     */
    private static Matrix schmidtChainRates(NetworkStruct sn, Matrix STchain) {
        if (!SnHasClassSwitching.snHasClassSwitching(sn)) {
            return null;
        }
        Matrix out = new Matrix(STchain.getNumRows(), STchain.getNumCols());
        for (int i = 0; i < STchain.getNumRows(); i++) {
            for (int c = 0; c < STchain.getNumCols(); c++) {
                double st = STchain.get(i, c);
                out.set(i, c, (st > 0 && Double.isFinite(st)) ? 1.0 / st : 0.0);
            }
        }
        return out;
    }
}
