/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.handlers;

import jline.GlobalConstants;
import jline.api.pfqn.ld.*;
import jline.api.sn.SnDeaggregateChainResults;
import jline.api.sn.SnGetDemandsChain;
import jline.io.InputOutput;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.lang.nodes.Station;
import jline.util.SerializableFunction;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Handler for the solver_amvald function
 */
public final class Solver_amvald {

    private Solver_amvald() {
    }

    public static MVAResult solver_amvald(NetworkStruct sn, SolverOptions options) {
        if (options == null) {
            options = new SolverOptions(SolverType.MVA);
        }
        long startTime = System.nanoTime();

        Ret.snGetDemands res = SnGetDemandsChain.snGetDemandsChain(sn);
        Matrix Lchain = res.Dchain;
        Matrix STchain = res.STchain;
        Matrix Vchain = res.Vchain;
        Matrix alpha = res.alpha;
        Matrix Nchain = res.Nchain;
        Matrix SCVchain = res.SCVchain;
        Matrix refstatchain = res.refstatchain;

        int M = sn.nstations;
        int K = sn.nchains;
        double Nt = 0.0;
        for (int col = 0; col < Nchain.getNumCols(); col++) {
            if (Double.isFinite(Nchain.get(0, col))) Nt += Nchain.get(0, col);
        }
        double tol = options.iter_tol;
        Matrix nservers = sn.nservers;
        Matrix schedparam = sn.schedparam;
        Matrix lldscaling = sn.lldscaling;
        Map<Station, SerializableFunction<Matrix, Matrix>> cdscaling = sn.cdscaling;
        Map<Station, SerializableFunction<Matrix, Matrix>> jdscaling = sn.jdscaling;
        List<Station> stations = sn.stations;
        List jobClasses = sn.jobclasses;

        List<SchedStrategy> sched = new ArrayList<SchedStrategy>();
        for (int i = 0; i < M; i++) {
            sched.add(sn.sched.get(sn.stations.get(i)));
        }

        Matrix Uchain = new Matrix(M, K);
        Matrix Tchain = new Matrix(M, K);
        Matrix Rchain = new Matrix(M, K);
        Matrix Cchain_s = new Matrix(1, K);

        Matrix Qchain = (options.init_sol != null) ? options.init_sol.copy() : null;
        if (Qchain != null && !Qchain.isEmpty() && (Qchain.getNumRows() != M || Qchain.getNumCols() != K)) {
            // stale warm-start hint from a different station/chain basis
            Qchain = null;
        }
        if ((Qchain == null) || (Qchain.isEmpty())) {
            // balanced initialization to match MATLAB implementation
            Qchain = Matrix.ones(M, K);

            // Normalize by column sums and multiply by Nchain values
            for (int col = 0; col < K; col++) {
                double colSum = Qchain.sumCols(col);
                for (int row = 0; row < M; row++) {
                    Qchain.set(row, col, (Qchain.get(row, col) / colSum) * Nchain.get(col));
                }
            }

            // Set infinite values to 0 (open classes)
            Qchain.apply(GlobalConstants.Inf, 0.0, "equal");

            // For open classes, set reference station queue length to 0
            for (int r = 0; r < Nchain.getNumCols(); r++) {
                if (Double.isInfinite(Nchain.get(r))) {
                    Qchain.set((int) refstatchain.get(r), r, 0.0);
                }
            }
        }

        List<Integer> nnzclasses = new ArrayList<Integer>();
        Matrix Xchain = new Matrix(1, STchain.getNumCols());
        for (int r = 0; r < Nchain.getNumCols(); r++) {
            if (Double.isInfinite(Nchain.get(0, r))) {
                Xchain.set(0, r, 1.0 / STchain.get((int) refstatchain.get(r, 0), r));
            } else {
                Xchain.set(0, r, 1.0 / STchain.sumCols(r));
            }

            if (Nchain.get(0, r) > 0) nnzclasses.add(r);
        }

        for (int k = 0; k < M; k++) {
            for (Integer r : nnzclasses) {
                if (Double.isInfinite(nservers.get(k, 0))) {
                    Uchain.set(k, r, Vchain.get(k, r) * STchain.get(k, r) * Xchain.get(0, r));
                } else {
                    Uchain.set(k, r, (Vchain.get(k, r) * STchain.get(k, r) * Xchain.get(0, r)) / nservers.get(k, 0));
                }
            }
        }

        if (options.config.np_priority == null) options.config.np_priority = "default";
        if (options.config.multiserver == null) options.config.multiserver = "default";
        if (options.config.highvar == null) options.config.highvar = "default";

        if ("default".equals(options.method)) {
            if (Nt <= 2) options.method = "bs";
            else options.method = "lin";
        }

        // Use list JLineMatrix to represent 3-D Matrix
        List<Matrix> gamma = new ArrayList<Matrix>();
        Matrix tau = new Matrix(K, K);
        if ("default".equals(options.method) || "amva_lin".equals(options.method) || "lin".equals(options.method)
                || "amva_qdlin".equals(options.method) || "qdlin".equals(options.method)) {
            int i = 0;
            while (i < K) {
                gamma.add(new Matrix(K, M));
                i++;
            }
        } else {
            gamma.add(new Matrix(K, M));
        }

        /* Main Loop */
        double omicron = 0.5; // under-relaxation parameter
        double outer_iter = 0.0;
        int totiter = 0;
        Matrix QchainOuter_1 = Qchain.copy();
        Matrix XchainOuter_1 = Xchain.copy();
        Matrix UchainOuter_1 = Uchain.copy();
        // Initialized to STchain so an early wall-clock timeout (options.timeout)
        // that breaks before any forward evaluation still yields a valid result.
        Matrix STeff = STchain.copy();
        while ((outer_iter < 2 || Qchain.sub(QchainOuter_1).elementMaxAbs() > tol)
                && (outer_iter < FastMath.sqrt((double) options.iter_max))
                && !Solver.timeExceeded(startTime, options.timeout)) {
            outer_iter++;

            QchainOuter_1 = Qchain.copy();
            XchainOuter_1 = Xchain.copy();
            UchainOuter_1 = Uchain.copy();

            if (Double.isFinite(Nt) && Nt > 0) {
                if ("default".equals(options.method) || "lin".equals(options.method) || "qdlin".equals(options.method)) {
                    /* Iteration at population N-1_s */
                    for (int s = 0; s < K; s++) {
                        if (Double.isFinite(Nchain.get(0, s))) {
                            double iter_s = 0.0;
                            List<Integer> sList = new ArrayList<Integer>();
                            sList.add(s);
                            Matrix Nchain_s = Matrix.oner(Nchain, new ArrayList<Integer>(sList));
                            Matrix Qchain_s = Qchain.copy();
                            Qchain_s.scaleEq((Nt - 1) / Nt, Qchain_s);
                            Matrix Xchain_s = Xchain.copy();
                            Xchain_s.scaleEq((Nt - 1) / Nt, Xchain_s);
                            Matrix Uchain_s = Uchain.copy();
                            Uchain_s.scaleEq((Nt - 1) / Nt, Uchain_s);
                            Matrix Qchain_s_1 = Qchain_s.copy();
                            Matrix Xchain_s_1 = Xchain_s.copy();
                            Matrix Uchain_s_1 = Uchain_s.copy();

                            int count_inf = 0;
                            double Nt_s = 0.0;
                            Matrix deltaclass = new Matrix(Nchain_s.getNumRows(), Nchain_s.getNumCols());
                            List<Integer> ocl = new ArrayList<Integer>();
                            List<Integer> ccl = new ArrayList<Integer>();
                            List<Integer> nnzclasses_s = new ArrayList<Integer>();
                            for (int col = 0; col < Nchain_s.getNumCols(); col++) {
                                double nchainValue = Nchain_s.get(0, col);
                                if (Double.isInfinite(nchainValue)) {
                                    count_inf++;
                                    deltaclass.set(0, col, 1.0);
                                    ocl.add(col);
                                } else {
                                    Nt_s += nchainValue;
                                    deltaclass.set(0, col, (nchainValue - 1) / nchainValue);
                                    if (nchainValue > 0) ccl.add(col);
                                }

                                if (nchainValue > 0) nnzclasses_s.add(col);
                            }
                            double delta = (count_inf == Nchain_s.getNumCols()) ? 1.0 : (Nt_s - 1) / Nt_s;

                            Map<Integer, List<Integer>> nnzclasses_eprio = new HashMap<Integer, List<Integer>>(nnzclasses_s.size());
                            Map<Integer, List<Integer>> nnzclasses_hprio = new HashMap<Integer, List<Integer>>(nnzclasses_s.size());
                            Map<Integer, List<Integer>> nnzclasses_ehprio = new HashMap<Integer, List<Integer>>(nnzclasses_s.size());
                            for (Integer r : nnzclasses_s) {
                                double prio = sn.classprio.get(0, r);
                                List<Integer> eprio_list = new ArrayList<Integer>();
                                List<Integer> hprio_list = new ArrayList<Integer>();
                                for (int i = 0; i < sn.classprio.getNumCols(); i++) {
                                    if (Double.compare(prio, sn.classprio.get(0, i)) == 0) eprio_list.add(i);
                                    else if (Double.compare(prio, sn.classprio.get(0, i)) > 0) hprio_list.add(i);
                                }
                                List<Integer> eprio_common = new ArrayList<Integer>(nnzclasses_s);
                                List<Integer> hprio_common = new ArrayList<Integer>(nnzclasses_s);
                                eprio_common.retainAll(eprio_list);
                                hprio_common.retainAll(hprio_list);
                                nnzclasses_eprio.put(r, eprio_common);
                                nnzclasses_hprio.put(r, hprio_common);

                                Set<Integer> ehprio_common = new LinkedHashSet<Integer>(eprio_common);
                                ehprio_common.addAll(hprio_common);
                                nnzclasses_ehprio.put(r, new ArrayList<Integer>(ehprio_common));
                            }
                            while ((iter_s < 2 || Qchain_s.sub(Qchain_s_1).elementMaxAbs() > tol)
                                    && (iter_s <= FastMath.sqrt((double) options.iter_max))
                                    && !Solver.timeExceeded(startTime, options.timeout)) {
                                iter_s++;

                                Qchain_s_1 = Qchain_s.copy();
                                Xchain_s_1 = Xchain_s.copy();
                                Uchain_s_1 = Uchain_s.copy();

                                Pair<Matrix, Matrix> ret2 = solver_amvald_forward(gamma,
                                        tau,
                                        Qchain_s_1,
                                        Xchain_s_1,
                                        Uchain_s_1,
                                        STchain,
                                        Vchain,
                                        Nchain_s,
                                        SCVchain,
                                        Nt_s,
                                        delta,
                                        deltaclass,
                                        ocl,
                                        ccl,
                                        nnzclasses_s,
                                        nnzclasses_eprio,
                                        nnzclasses_hprio,
                                        nnzclasses_ehprio,
                                        M,
                                        K,
                                        nservers,
                                        schedparam,
                                        lldscaling,
                                        cdscaling,
                                        jdscaling,
                                        sched,
                                        stations,
                                        options);
                                totiter++;
                                if (totiter > options.iter_max) {
                                    break;
                                }
                                Matrix Wchain_s = ret2.getLeft();
                                Matrix STeff_s = ret2.getRight();

                                for (Integer r : nnzclasses) {
                                    if (Wchain_s.sumCols(r) == 0.0) {
                                        Xchain_s.remove(0, r);
                                    } else {
                                        if (Double.isInfinite(Nchain_s.get(0, r))) {
                                            double sumValue = 0.0;
                                            for (int row = 0; row < Vchain.getNumRows(); row++) {
                                                sumValue += Vchain.get(row, r) * Wchain_s.get(row, r);
                                            }
                                            Cchain_s.set(0, r, sumValue);
                                        } else if (Nchain.get(0, r) == 0.0) {
                                            Xchain_s.remove(0, r);
                                            Cchain_s.remove(0, r);
                                        } else {
                                            double sumValue = 0.0;
                                            for (int row = 0; row < Vchain.getNumRows(); row++) {
                                                sumValue += Vchain.get(row, r) * Wchain_s.get(row, r);
                                            }
                                            Cchain_s.set(0, r, sumValue);
                                            Xchain_s.set(0, r,
                                                    omicron * Nchain_s.get(0, r) / Cchain_s.get(0, r) + (1 - omicron) * Xchain_s_1.get(0, r));
                                        }
                                    }
                                    for (int k = 0; k < M; k++) {
                                        Qchain_s.set(k, r,
                                                omicron * Xchain_s.get(0, r) * Vchain.get(k, r) * Wchain_s.get(k, r) + (1 - omicron) * Qchain_s_1.get(k, r));
                                        Uchain_s.set(k, r,
                                                omicron * Vchain.get(k, r) * STeff_s.get(k, r) * Xchain_s.get(0, r) + (1 - omicron) * Uchain_s_1.get(k, r));
                                    }
                                }
                            }

                            if ("default".equals(options.method) || "lin".equals(options.method)) {
                                int k = 0;
                                while (k < M) {
                                    for (Integer r : nnzclasses) {
                                        if (Double.isFinite(Nchain.get(0, r)) && Nchain_s.get(0, r) > 0) {
                                            gamma.get(r).set(s, k,
                                                    Qchain_s_1.get(k, r) / Nchain_s.get(0, r) - QchainOuter_1.get(k, r) / Nchain.get(0, r));
                                        }
                                    }
                                    k++;
                                }
                            } else {
                                int k = 0;
                                while (k < M) {
                                    gamma.get(0).set(s, k,
                                            Qchain_s_1.sumRows(k) / (Nt - 1) - QchainOuter_1.sumRows(k) / Nt);
                                    k++;
                                }
                            }

                            for (Integer r : nnzclasses) {
                                tau.set(s, r, Xchain_s_1.get(0, r) - XchainOuter_1.get(0, r));
                            }
                        }
                    }
                }
            }

            double inner_iter = 0.0;
            Matrix Qchain_inner = Qchain.copy();
            Matrix Xchain_inner = Xchain.copy();
            Matrix Uchain_inner = Uchain.copy();

            int count_inf = 0;
            double Nt_inner = 0.0;
            Matrix deltaclass = new Matrix(Nchain.getNumRows(), Nchain.getNumCols());
            List<Integer> ocl = new ArrayList<Integer>();
            List<Integer> ccl = new ArrayList<Integer>();
            List<Integer> nnzclasses_inner = new ArrayList<Integer>();
            for (int col = 0; col < Nchain.getNumCols(); col++) {
                double nchainValue = Nchain.get(0, col);
                if (Double.isInfinite(nchainValue)) {
                    count_inf++;
                    deltaclass.set(0, col, 1.0);
                    ocl.add(col);
                } else {
                    Nt_inner += nchainValue;
                    deltaclass.set(0, col, (nchainValue - 1) / nchainValue);
                    if (nchainValue > 0) ccl.add(col);
                }

                if (nchainValue > 0) nnzclasses_inner.add(col);
            }
            double delta = (count_inf == Nchain.getNumCols()) ? 1.0 : (Nt_inner - 1) / Nt_inner;

            Map<Integer, List<Integer>> nnzclasses_eprio = new HashMap<Integer, List<Integer>>(nnzclasses_inner.size());
            Map<Integer, List<Integer>> nnzclasses_hprio = new HashMap<Integer, List<Integer>>(nnzclasses_inner.size());
            Map<Integer, List<Integer>> nnzclasses_ehprio = new HashMap<Integer, List<Integer>>(nnzclasses_inner.size());
            for (Integer r : nnzclasses_inner) {
                double prio = sn.classprio.get(0, r);
                List<Integer> eprio_list = new ArrayList<Integer>();
                List<Integer> hprio_list = new ArrayList<Integer>();
                for (int i = 0; i < sn.classprio.getNumCols(); i++) {
                    if (Double.compare(prio, sn.classprio.get(0, i)) == 0) eprio_list.add(i);
                    else if (Double.compare(prio, sn.classprio.get(0, i)) > 0) hprio_list.add(i);
                }
                List<Integer> eprio_common = new ArrayList<Integer>(nnzclasses_inner);
                List<Integer> hprio_common = new ArrayList<Integer>(nnzclasses_inner);
                eprio_common.retainAll(eprio_list);
                hprio_common.retainAll(hprio_list);
                nnzclasses_eprio.put(r, eprio_common);
                nnzclasses_hprio.put(r, hprio_common);

                Set<Integer> ehprio_common = new LinkedHashSet<Integer>(eprio_common);
                ehprio_common.addAll(hprio_common);
                nnzclasses_ehprio.put(r, new ArrayList<Integer>(ehprio_common));
            }
            while ((inner_iter < 2 || Qchain_inner.sub(Qchain).elementMaxAbs() > tol)
                    && (inner_iter <= FastMath.sqrt((double) options.iter_max))
                    && !Solver.timeExceeded(startTime, options.timeout)) {
                inner_iter++;

                Qchain_inner = Qchain.copy();
                Xchain_inner = Xchain.copy();
                Uchain_inner = Uchain.copy();

                Pair<Matrix, Matrix> ret2 = solver_amvald_forward(gamma,
                        tau,
                        Qchain_inner,
                        Xchain_inner,
                        Uchain_inner,
                        STchain,
                        Vchain,
                        Nchain,
                        SCVchain,
                        Nt_inner,
                        delta,
                        deltaclass,
                        ocl,
                        ccl,
                        nnzclasses_inner,
                        nnzclasses_eprio,
                        nnzclasses_hprio,
                        nnzclasses_ehprio,
                        M,
                        K,
                        nservers,
                        schedparam,
                        lldscaling,
                        cdscaling,
                        jdscaling,
                        sched,
                        stations,
                        options);
                totiter++;
                if (totiter > options.iter_max) {
                    break;
                }
                Matrix Wchain = ret2.getLeft();
                STeff = ret2.getRight();

                for (Integer r : nnzclasses) {
                    if (Wchain.sumCols(r) == 0.0) {
                        Xchain.remove(0, r);
                    } else {
                        if (Double.isInfinite(Nchain.get(0, r))) {
                            double sumValue = 0.0;
                            for (int i = 0; i < Vchain.getNumRows(); i++) sumValue += Vchain.get(i, r) * Wchain.get(i, r);
                            Cchain_s.set(0, r, sumValue);
                        } else if (Nchain.get(0, r) == 0.0) {
                            Xchain.remove(0, r);
                            Cchain_s.remove(0, r);
                        } else {
                            double sumValue = 0.0;
                            for (int i = 0; i < Vchain.getNumRows(); i++) sumValue += Vchain.get(i, r) * Wchain.get(i, r);
                            Cchain_s.set(0, r, sumValue);
                            Xchain.set(0, r,
                                    omicron * Nchain.get(0, r) / Cchain_s.get(0, r) + (1 - omicron) * Xchain_inner.get(0, r));
                        }
                    }
                    for (int k = 0; k < M; k++) {
                        Qchain.set(k, r,
                                omicron * Xchain.get(0, r) * Vchain.get(k, r) * Wchain.get(k, r) + (1 - omicron) * Qchain_inner.get(k, r));
                        Tchain.set(k, r, Xchain.get(0, r) * Vchain.get(k, r));
                        Uchain.set(k, r,
                                omicron * Vchain.get(k, r) * STeff.get(k, r) * Xchain.get(0, r) + (1 - omicron) * Uchain_inner.get(k, r));
                    }
                }
            }
        }

        for (int k = 0; k < M; k++) {
            for (int r = 0; r < K; r++) {
                if (Vchain.get(k, r) * STeff.get(k, r) > 0) {
                    SchedStrategy schedK = sn.sched.get(sn.stations.get(k));
                    if (schedK == SchedStrategy.FCFS || schedK == SchedStrategy.SIRO || schedK == SchedStrategy.PS
                            || schedK == SchedStrategy.LCFSPR || schedK == SchedStrategy.DPS || schedK == SchedStrategy.HOL
                            || schedK == SchedStrategy.FCFSPRIO) {
                        if (Uchain.sumRows(k) > 1) {
                            double sum_vchain_steff_xchain_k = 0.0;
                            int i = 0;
                            while (i < Vchain.getNumCols()) {
                                sum_vchain_steff_xchain_k += Vchain.get(k, i) * STeff.get(k, i) * Xchain.get(0, i);
                                i++;
                            }
                            Uchain.set(k, r,
                                    (Math.min(Uchain.sumRows(k), 1.0) * Vchain.get(k, r) * STeff.get(k, r) * Xchain.get(0, r)) / sum_vchain_steff_xchain_k);
                        }
                    }
                }
            }
        }

        for (int k = 0; k < M; k++) {
            for (int r = 0; r < K; r++) {
                if (Qchain.get(k, r) < GlobalConstants.Zero) {
                    Rchain.set(k, r, 0.0);
                } else {
                    Rchain.set(k, r, Qchain.get(k, r) / Tchain.get(k, r));
                }
            }
        }
        Xchain.apply(Double.POSITIVE_INFINITY, 0.0, "equal");
        Xchain.apply(Double.NaN, 0.0, "equal");
        Uchain.apply(Double.POSITIVE_INFINITY, 0.0, "equal");
        Uchain.apply(Double.NaN, 0.0, "equal");
        Rchain.apply(Double.POSITIVE_INFINITY, 0.0, "equal");
        Rchain.apply(Double.NaN, 0.0, "equal");

        for (int col = 0; col < K; col++) {
            if (Nchain.get(0, col) == 0.0) {
                Xchain.remove(0, col);
                for (int row = 0; row < M; row++) {
                    Uchain.remove(row, col);
                    Rchain.remove(row, col);
                    Tchain.remove(row, col);
                }
            }
        }

        Ret.snDeaggregateChainResults ret;
        if ((sn.lldscaling == null || sn.lldscaling.isEmpty()) && (sn.cdscaling == null || sn.cdscaling.size() == 0) && (sn.jdscaling == null || sn.jdscaling.size() == 0)) {
            ret = SnDeaggregateChainResults.snDeaggregateChainResults(sn, Lchain, null, STchain, Vchain, alpha, null, null, Rchain, Tchain, null, Xchain);
        } else {
            ret = SnDeaggregateChainResults.snDeaggregateChainResults(sn, Lchain, null, STchain, Vchain, alpha, null, Uchain, Rchain, Tchain, null, Xchain);
        }

        // see _kb/06-solver-catalog.md for rationale
        if (sn.cdscaling != null && sn.cdscaling.size() > 0) {
            boolean allFinite = true;
            for (int k = 0; k < sn.nclasses; k++) {
                if (Double.isInfinite(sn.njobs.get(k))) { allFinite = false; break; }
            }
            if (allFinite) {
                int Kcls = sn.nclasses; // ret.U and rates are class-indexed; K counts chains
                for (int ist = 0; ist < M; ist++) {
                    Station stat = sn.stations.get(ist);
                    SerializableFunction<Matrix, Matrix> beta = sn.cdscaling.get(stat);
                    if (beta == null) continue;
                    Matrix peakVec = sn.cdscalingpeak != null ? sn.cdscalingpeak.get(stat) : null;
                    for (int k = 0; k < Kcls; k++) {
                        double rate = sn.rates.get(ist, k);
                        double bmax = (peakVec != null) ? peakVec.get(0, k) : 1.0;
                        if (Double.isFinite(rate) && rate > 0 && bmax > 0) {
                            ret.U.set(ist, k, ret.T.get(ist, k) / rate / bmax);
                        } else {
                            ret.U.set(ist, k, 0.0);
                        }
                    }
                }
            }
        }

        // joint-dependence utilization normalization: U = T*S/peak using the
        // declared sn.jdscalingpeak, mirroring the class-dependence block.
        if (sn.jdscaling != null && sn.jdscaling.size() > 0) {
            boolean allFinite = true;
            for (int k = 0; k < sn.nclasses; k++) {
                if (Double.isInfinite(sn.njobs.get(k))) { allFinite = false; break; }
            }
            if (allFinite) {
                int Kcls = sn.nclasses;
                for (int ist = 0; ist < M; ist++) {
                    Station stat = sn.stations.get(ist);
                    SerializableFunction<Matrix, Matrix> eta = sn.jdscaling.get(stat);
                    if (eta == null) continue;
                    Matrix peakVec = sn.jdscalingpeak != null ? sn.jdscalingpeak.get(stat) : null;
                    for (int k = 0; k < Kcls; k++) {
                        double rate = sn.rates.get(ist, k);
                        double bmax = (peakVec != null) ? peakVec.get(0, k) : 1.0;
                        if (Double.isFinite(rate) && rate > 0 && bmax > 0) {
                            ret.U.set(ist, k, ret.T.get(ist, k) / rate / bmax);
                        } else {
                            ret.U.set(ist, k, 0.0);
                        }
                    }
                }
            }
        }

        List<Integer> ccl = new ArrayList<Integer>();
        for (int i = 0; i < Nchain.getNumCols(); i++) {
            if (Double.isFinite(Nchain.get(0, i))) ccl.add(i);
        }
        Matrix Nclosed = new Matrix(1, ccl.size());
        Matrix Xclosed = new Matrix(1, ccl.size());
        for (int i = 0; i < ccl.size(); i++) {
            Nclosed.set(0, i, Nchain.get(0, ccl.get(i)));
            Xclosed.set(0, i, Xchain.get(0, ccl.get(i)));
        }
        double lG = 0.0;
        for (int i = 0; i < ccl.size(); i++) {
            if (Xclosed.get(0, i) > options.tol) lG += -Nclosed.get(0, i) * FastMath.log(Xclosed.get(0, i));
        }

        long endTime = System.nanoTime();
        long runTime = endTime - startTime;

        MVAResult result = new MVAResult();
        result.method = options.method;
        result.QN = ret.Q;
        result.RN = ret.R;
        result.XN = ret.X;
        result.UN = ret.U;
        result.TN = ret.T;
        result.CN = ret.C;
        result.runtime = runTime / 1000000000.0;
        result.logNormConstAggr = lG;
        // see _kb/06-solver-catalog.md for rationale
        result.iter = totiter;
        // see _kb/06-solver-catalog.md for rationale
        result.converged = Boolean.valueOf(Qchain.sub(QchainOuter_1).elementMaxAbs() <= tol);
        return result;
    }

    /**
     * Column-wise mean of the per-chain correction {@code g}, as an M x 1 column.
     *
     * <p>MATLAB seeds the correction as the SCALAR {@code g = 0} and accumulates
     * one row per closed chain, so with no closed chains {@code mean(g,1)} is
     * {@code 0} and the correction vanishes. Here {@code g} is allocated with one
     * row per closed chain, which is a 0 x M matrix in that case, and averaging
     * zero rows would yield NaN (0/0) rather than zero. That NaN propagates
     * through msterm into STeff and Wchain, and is only caught by the terminal
     * non-finite sanitisation, surfacing as an all-zero result.
     */
    private static Matrix meanColOrZeros(Matrix g, int M) {
        if (g.getNumRows() == 0) {
            return new Matrix(M, 1);
        }
        return g.meanCol().transpose();
    }

    public static Pair<Matrix, Matrix> solver_amvald_forward(List<Matrix> gamma,
                                                             Matrix tau,
                                                             Matrix Qchain_in,
                                                             Matrix Xchain_in,
                                                             Matrix Uchain_in,
                                                             Matrix STchain_in,
                                                             Matrix Vchain_in,
                                                             Matrix Nchain_in,
                                                             Matrix SCVchain_in,
                                                             double Nt,
                                                             double delta,
                                                             Matrix deltaclass,
                                                             List<Integer> ocl,
                                                             List<Integer> ccl,
                                                             List<Integer> nnzclasses,
                                                             Map<Integer, List<Integer>> nnzclasses_eprio,
                                                             Map<Integer, List<Integer>> nnzclasses_hprio,
                                                             Map<Integer, List<Integer>> nnzclasses_ehprio,
                                                             int M,
                                                             int K,
                                                             Matrix nservers,
                                                             Matrix schedparam,
                                                             Matrix lldscaling_in,
                                                             Map<Station, SerializableFunction<Matrix, Matrix>> cdscaling,
                                                             Map<Station, SerializableFunction<Matrix, Matrix>> jdscaling,
                                                             List<SchedStrategy> sched,
                                                             List<Station> stations,
                                                             SolverOptions options) {
        Matrix lldscaling;

        if (gamma.size() == 0) {
            gamma = new ArrayList<Matrix>();
            gamma.add(new Matrix(K, M));
        }

        /* Evaluate lld and cd correction factors */
        Matrix totArvlQlenSeenByOpen = new Matrix(K, M);
        Matrix interpTotArvlQlen = new Matrix(M, 1);
        Matrix totArvlQlenSeenByClosed = new Matrix(M, K);
        Matrix stationaryQlen = new Matrix(M, K);
        Matrix selfArvlQlenSeenByClosed = new Matrix(M, K);
        for (int k = 0; k < M; k++) {
            double sum_qchain_k_nnz = 0.0;
            for (Integer r : nnzclasses) sum_qchain_k_nnz += Qchain_in.get(k, r);

            interpTotArvlQlen.set(k, 0, delta * sum_qchain_k_nnz);
            for (Integer r : nnzclasses) {
                selfArvlQlenSeenByClosed.set(k, r, deltaclass.get(r) * Qchain_in.get(k, r));
                if (sched.get(k) == SchedStrategy.HOL || sched.get(k) == SchedStrategy.FCFSPRIO) {
                    double sum_qchain_k_ehprio = 0.0;
                    for (Integer i : nnzclasses_ehprio.get(r)) sum_qchain_k_ehprio += Qchain_in.get(k, i);
                    totArvlQlenSeenByOpen.set(r, k, sum_qchain_k_ehprio);
                    totArvlQlenSeenByClosed.set(k, r,
                            deltaclass.get(r) * Qchain_in.get(k, r) + sum_qchain_k_ehprio - Qchain_in.get(k, r));
                } else {
                    totArvlQlenSeenByOpen.set(r, k, sum_qchain_k_nnz);
                    totArvlQlenSeenByClosed.set(k, r,
                            deltaclass.get(r) * Qchain_in.get(k, r) + sum_qchain_k_nnz - Qchain_in.get(k, r));
                }
                stationaryQlen.set(k, r, Qchain_in.get(k, r));
            }
        }

        if (lldscaling_in == null || lldscaling_in.isEmpty()) {
            lldscaling = new Matrix(M, 1);
            lldscaling.fill(1.0);
        } else {
            lldscaling = lldscaling_in.copy();
        }
        boolean isLin = "lin".equals(options.method) || "amva_lin".equals(options.method)
                || "qdlin".equals(options.method) || "amva_qdlin".equals(options.method);

        Matrix lldterm;
        if (isLin && !nnzclasses.isEmpty() && !ccl.isEmpty()) {
            lldterm = new Matrix(M, K);
            lldterm.fill(1.0);
            for (Integer r : nnzclasses) {
                // Linearizer fraction correction seen by an arriving class-r job,
                // evaluated at every station: entry k uses gamma(r,k,:).
                Matrix arg = new Matrix(M, 1);
                for (int k = 0; k < M; k++) {
                    double gcorr = 0.0;
                    for (Integer c : ccl) gcorr += Nchain_in.get(0, c) * gamma.get(c).get(r, k);
                    gcorr -= gamma.get(r).get(r, k);
                    arg.set(k, 0, 1.0 + interpTotArvlQlen.get(k, 0) + gcorr);
                }
                Matrix.extract(Pfqn_lldfun.pfqn_lldfun(arg, lldscaling, null), 0, M, 0, 1, lldterm, 0, r);
            }
        } else {
            // see _kb/06-solver-catalog.md for rationale
            lldterm = Pfqn_lldfun.pfqn_lldfun(interpTotArvlQlen.elementIncrease(1.0), lldscaling, null);
        }

        Matrix cdterm = new Matrix(M, K);
        cdterm.fill(1.0);
        for (Integer r : nnzclasses) {
            if (!(cdscaling == null || cdscaling.size() == 0)) {
                List<SerializableFunction<Matrix, Matrix>> cdscalingList = new ArrayList<SerializableFunction<Matrix, Matrix>>();
                for (int i = 0; i < M; i++) {
                    cdscalingList.add(cdscaling != null ? cdscaling.get(stations.get(i)) : null);
                }
                // qd-amva class-dependence term beta_{i,r}
                if (Double.isFinite(Nchain_in.get(0, r))) {
                    Matrix nvec = selfArvlQlenSeenByClosed.elementIncrease(1.0);
                    if (isLin) {
                        // self-fraction correction per station: entry k uses gamma(r,k,r)
                        double Nr = Nchain_in.get(0, r);
                        for (int k = 0; k < M; k++) {
                            double gself = (Nr - 1.0) * gamma.get(r).get(r, k);
                            for (int c = 0; c < K; c++) nvec.set(k, c, nvec.get(k, c) + gself);
                        }
                    }
                    Matrix.extract(Pfqn_cdfun.pfqn_cdfun(nvec, cdscalingList, M, r), 0, M, 0, 1, cdterm, 0, r);
                } else {
                    // gamma is zero for open classes, so the correction vanishes
                    Matrix.extract(Pfqn_cdfun.pfqn_cdfun(stationaryQlen.elementIncrease(1.0), cdscalingList, M, r),
                            0, M, 0, 1, cdterm, 0, r);
                }
            }
        }

        // joint-dependence term eta_i (non-product-form). Evaluated identically
        // to cdterm; kept separate to preserve the product-form vs joint
        // distinction (see Pfqn_jdfun). cd and jd are mutually exclusive per
        // station, so at most one of cdterm/jdterm differs from 1 at any station.
        Matrix jdterm = new Matrix(M, K);
        jdterm.fill(1.0);
        for (Integer r : nnzclasses) {
            if (!(jdscaling == null || jdscaling.size() == 0)) {
                List<SerializableFunction<Matrix, Matrix>> jdscalingList = new ArrayList<SerializableFunction<Matrix, Matrix>>();
                for (int i = 0; i < M; i++) {
                    jdscalingList.add(jdscaling != null ? jdscaling.get(stations.get(i)) : null);
                }
                if (Double.isFinite(Nchain_in.get(0, r))) {
                    Matrix nvec = selfArvlQlenSeenByClosed.elementIncrease(1.0);
                    if (isLin) {
                        double Nr = Nchain_in.get(0, r);
                        for (int k = 0; k < M; k++) {
                            double gself = (Nr - 1.0) * gamma.get(r).get(r, k);
                            for (int c = 0; c < K; c++) nvec.set(k, c, nvec.get(k, c) + gself);
                        }
                    }
                    Matrix.extract(Pfqn_jdfun.pfqn_jdfun(nvec, jdscalingList, M, r), 0, M, 0, 1, jdterm, 0, r);
                } else {
                    Matrix.extract(Pfqn_jdfun.pfqn_jdfun(stationaryQlen.elementIncrease(1.0), jdscalingList, M, r),
                            0, M, 0, 1, jdterm, 0, r);
                }
            }
        }

        Matrix msterm = new Matrix(0, 0);
        if ("softmin".equals(options.config.multiserver)) {
            if ("default".equals(options.method) || "amva_lin".equals(options.method) || "lin".equals(options.method)
                    || "amva_qdlin".equals(options.method) || "qdlin".equals(options.method)) {
                Matrix g = new Matrix(ccl.size(), M);
                for (Integer r : ccl) {
                    double param = ((Nt - 1.0) / Nt) * Nchain_in.get(0, r);
                    Matrix gamma_r = gamma.get(r);
                    int i = 0;
                    while (i < ccl.size()) {
                        int idx = ccl.get(i);
                        int j = 0;
                        while (j < M) {
                            g.set(i, j, g.get(i, j) + param * gamma_r.get(idx, j));
                            j++;
                        }
                        i++;
                    }
                }
                Matrix input = interpTotArvlQlen.add(1.0, meanColOrZeros(g, M));
                msterm = Pfqn_lldfun.pfqn_lldfun(input.elementIncrease(1.0), new Matrix(0, 0), nservers);
            } else {
                Matrix g = new Matrix(ccl.size(), M);
                int i = 0;
                while (i < ccl.size()) {
                    int idx = ccl.get(i);
                    int j = 0;
                    while (j < M) {
                        g.set(i, j, g.get(i, j) + (Nt - 1.0) * gamma.get(0).get(idx, j));
                        j++;
                    }
                    i++;
                }
                // see _kb/06-solver-catalog.md for rationale
                Matrix input = interpTotArvlQlen.add(1.0, meanColOrZeros(g, M));
                msterm = Pfqn_lldfun.pfqn_lldfun(input.elementIncrease(1.0), new Matrix(0, 0), nservers);
            }
        } else if ("seidmann".equals(options.config.multiserver)) {
            nservers.divide(1.0, msterm, false);
            msterm.apply(0.0, 1.0, "equal");
        } else if ("default".equals(options.config.multiserver) || "rolia".equals(options.config.multiserver)) {
            if ("default".equals(options.method) || "amva_lin".equals(options.method) || "lin".equals(options.method)
                    || "amva_qdlin".equals(options.method) || "qdlin".equals(options.method)) {
                Matrix g = new Matrix(ccl.size(), M);
                for (Integer r : ccl) {
                    double param = ((Nt - 1.0) / Nt) * Nchain_in.get(0, r);
                    Matrix gamma_r = gamma.get(r);
                    int i = 0;
                    while (i < ccl.size()) {
                        int idx = ccl.get(i);
                        int j = 0;
                        while (j < M) {
                            g.set(i, j, g.get(i, j) + param * gamma_r.get(idx, j));
                            j++;
                        }
                        i++;
                    }
                }
                Matrix input = interpTotArvlQlen.add(1.0, meanColOrZeros(g, M));
                msterm = Pfqn_lldfun.pfqn_lldfun(input.elementIncrease(1.0), new Matrix(0, 0), nservers);
            } else {
                Matrix g = new Matrix(1, M);
                for (Integer r : ccl) {
                    int j = 0;
                    while (j < M) {
                        g.set(0, j, g.get(0, j) + (Nt - 1.0) * gamma.get(0).get(r, j));
                        j++;
                    }
                }
                msterm = Pfqn_lldfun.pfqn_lldfun(interpTotArvlQlen.elementIncrease(1 + g.meanRow().value()), new Matrix(0, 0), nservers);
            }
            int i = 0;
            while (i < M) {
                SchedStrategy schedStrategy = sched.get(i);
                if (schedStrategy == SchedStrategy.FCFS || schedStrategy == SchedStrategy.SIRO || schedStrategy == SchedStrategy.LCFSPR) {
                    msterm.set(i, 0, 1.0 / nservers.get(i, 0));
                }
                i++;
            }
        } else if ("suri".equals(options.config.multiserver)) {
            msterm = new Matrix(M, 1);
            msterm.fill(1.0); // No demand scaling; multiserver handled via suriFactor in Wchain
        } else {
            throw new RuntimeException("Unrecognized multiserver approximation method");
        }

        // Compute Suri correction factor per station
        Matrix suriFactor = new Matrix(M, 1);
        suriFactor.fill(1.0);
        if ("suri".equals(options.config.multiserver)) {
            double suriAlpha = 4.464;
            double suriBeta = 0.676;
            for (int k = 0; k < M; k++) {
                double c = nservers.get(k, 0);
                if (c > 1 && Double.isFinite(c)) {
                    double rhoK = Math.min(Uchain_in.sumRows(k) / c, 1.0 - options.tol);
                    if (rhoK > 0) {
                        suriFactor.set(k, 0, Math.pow(rhoK, suriAlpha * (Math.pow(c, suriBeta) - 1.0)) / c);
                    } else {
                        suriFactor.set(k, 0, 1.0 / c);
                    }
                } else if (Double.isInfinite(c) || c == 0.0) {
                    suriFactor.set(k, 0, 0.0);
                }
            }
        }

        Matrix Wchain = new Matrix(M, K);
        Matrix STeff = new Matrix(STchain_in.getNumRows(), STchain_in.getNumCols());
        lldterm = lldterm.repmat(1, K);
        for (Integer r : nnzclasses) {
            for (int k = 0; k < M; k++) {
                STeff.set(k, r, STchain_in.get(k, r) * lldterm.get(k, r) * msterm.get(k, 0) * cdterm.get(k, r) * jdterm.get(k, r));
            }
        }

        /* if amva.qli or amva.fli, update now totArvlQlenSeenByClosed with STeff */
        if ("amva_qli".equals(options.method) || "qli".equals(options.method)) {
            List<Integer> infset = new ArrayList<Integer>();
            int i = 0;
            while (i < M) {
                if (sched.get(i) == SchedStrategy.INF) infset.add(i);
                i++;
            }

            int k = 0;
            while (k < M) {
                if (sched.get(k) == SchedStrategy.HOL || sched.get(k) == SchedStrategy.FCFSPRIO) {
                    for (Integer r : nnzclasses) {
                        double sum_Qchain_k_r_ephrio = 0.0;
                        for (Integer ii : nnzclasses_ehprio.get(r)) sum_Qchain_k_r_ephrio += Qchain_in.get(k, ii);

                        if (Math.abs(Nchain_in.get(0, r) - 1) < 1e-20) {
                            totArvlQlenSeenByClosed.set(k, r, sum_Qchain_k_r_ephrio - Qchain_in.get(k, r));
                        } else {
                            double qlinum = STeff.get(k, r) * (1 + sum_Qchain_k_r_ephrio - Qchain_in.get(k, r));
                            double qliden = 0.0;
                            for (Integer ii : infset) qliden += STeff.get(ii, r);
                            int m = 0;
                            while (m < M) {
                                double sum_Qchain_m_r_ephrio = 0.0;
                                for (Integer ii : nnzclasses_ehprio.get(r)) sum_Qchain_m_r_ephrio += Qchain_in.get(m, ii);
                                qliden += STeff.get(m, r) * (1 + sum_Qchain_m_r_ephrio - Qchain_in.get(m, r));
                                m++;
                            }
                            totArvlQlenSeenByClosed.set(k, r,
                                    sum_Qchain_k_r_ephrio - (1 / (Nchain_in.get(0, r) - 1)) * (Qchain_in.get(k, r) - qlinum / qliden));
                        }
                    }
                } else {
                    for (Integer r : nnzclasses) {
                        double sum_Qchain_k_r_nnzclasses = 0.0;
                        for (Integer ii : nnzclasses) sum_Qchain_k_r_nnzclasses += Qchain_in.get(k, ii);

                        if (Math.abs(Nchain_in.get(0, r) - 1) < 1e-20) {
                            totArvlQlenSeenByClosed.set(k, r, sum_Qchain_k_r_nnzclasses - Qchain_in.get(k, r));
                        } else {
                            double qlinum = STeff.get(k, r) * (1 + sum_Qchain_k_r_nnzclasses - Qchain_in.get(k, r));
                            double qliden = 0.0;
                            for (Integer ii : infset) qliden += STeff.get(ii, r);
                            int m = 0;
                            while (m < M) {
                                double sum_Qchain_m_r_nnzclasses = 0.0;
                                for (Integer ii : nnzclasses) sum_Qchain_m_r_nnzclasses += Qchain_in.get(m, ii);
                                qliden += STeff.get(m, r) * (1 + sum_Qchain_m_r_nnzclasses - Qchain_in.get(m, r));
                                m++;
                            }
                            totArvlQlenSeenByClosed.set(k, r,
                                    sum_Qchain_k_r_nnzclasses - (1 / (Nchain_in.get(0, r) - 1)) * (Qchain_in.get(k, r) - qlinum / qliden));
                        }
                    }
                }
                k++;
            }
        } else if ("amva_fli".equals(options.method) || "fli".equals(options.method)) {
            List<Integer> infset = new ArrayList<Integer>();
            int i = 0;
            while (i < M) {
                if (sched.get(i) == SchedStrategy.INF) infset.add(i);
                i++;
            }

            int k = 0;
            while (k < M) {
                if (sched.get(k) == SchedStrategy.HOL || sched.get(k) == SchedStrategy.FCFSPRIO) {
                    for (Integer r : nnzclasses) {
                        double sum_Qchain_k_r_ephrio = 0.0;
                        for (Integer ii : nnzclasses_ehprio.get(r)) sum_Qchain_k_r_ephrio += Qchain_in.get(k, ii);

                        if (Math.abs(Nchain_in.get(0, r) - 1) < 1e-20) {
                            totArvlQlenSeenByClosed.set(k, r, sum_Qchain_k_r_ephrio - Qchain_in.get(k, r));
                        } else {
                            double qlinum = STeff.get(k, r) * (1 + sum_Qchain_k_r_ephrio - Qchain_in.get(k, r));
                            double qliden = 0.0;
                            for (Integer ii : infset) qliden += STeff.get(ii, r);
                            int m = 0;
                            while (m < M) {
                                double sum_Qchain_m_r_ephrio = 0.0;
                                for (Integer ii : nnzclasses_ehprio.get(r)) sum_Qchain_m_r_ephrio += Qchain_in.get(m, ii);
                                qliden += STeff.get(m, r) * (1 + sum_Qchain_m_r_ephrio - Qchain_in.get(m, r));
                                m++;
                            }
                            totArvlQlenSeenByClosed.set(k, r,
                                    sum_Qchain_k_r_ephrio - (2 / Nchain_in.get(0, r)) * (Qchain_in.get(k, r) + qlinum / qliden));
                        }
                    }
                } else {
                    for (Integer r : nnzclasses) {
                        double sum_Qchain_k_r_nnzclasses = 0.0;
                        for (Integer ii : nnzclasses) sum_Qchain_k_r_nnzclasses += Qchain_in.get(k, ii);

                        if (Math.abs(Nchain_in.get(0, r) - 1) < 1e-20) {
                            totArvlQlenSeenByClosed.set(k, r, sum_Qchain_k_r_nnzclasses - Qchain_in.get(k, r));
                        } else {
                            double qlinum = STeff.get(k, r) * (1 + sum_Qchain_k_r_nnzclasses - Qchain_in.get(k, r));
                            double qliden = 0.0;
                            for (Integer ii : infset) qliden += STeff.get(ii, r);
                            int m = 0;
                            while (m < M) {
                                double sum_Qchain_m_r_nnzclasses = 0.0;
                                for (Integer ii : nnzclasses) sum_Qchain_m_r_nnzclasses += Qchain_in.get(m, ii);
                                qliden += STeff.get(m, r) * (1 + sum_Qchain_m_r_nnzclasses - Qchain_in.get(m, r));
                                m++;
                            }
                            totArvlQlenSeenByClosed.set(k, r,
                                    sum_Qchain_k_r_nnzclasses - (2 / Nchain_in.get(0, r)) * (Qchain_in.get(k, r) + qlinum / qliden));
                        }
                    }
                }
                k++;
            }
        }

        /* Compute response time */
        for (Integer r : nnzclasses) {
            List<Integer> sd = new ArrayList<Integer>(nnzclasses);
            List<Integer> sdprio = new ArrayList<Integer>(nnzclasses_ehprio.get(r));
            sd.remove((Integer) r);
            sdprio.remove((Integer) r);

            for (int k = 0; k < M; k++) {
                SchedStrategy schedK = sched.get(k);
                if (schedK == SchedStrategy.INF) {
                    Wchain.set(k, r, STeff.get(k, r));
                } else if (schedK == SchedStrategy.PS) {
                    String method = options.method;
                    boolean isPSDefaultMethod = "def".equals(method) || "amva".equals(method) || "amva_qd".equals(method)
                            || "amva_qdamva".equals(method) || "qd".equals(method) || "qdamva".equals(method)
                            || "lin".equals(method) || "qdlin".equals(method);
                    if (isPSDefaultMethod) {
                        if ("seidmann".equals(options.config.multiserver)) {
                            double multiServerTerm = STeff.get(k, r) * (nservers.get(k, 0) - 1);
                            if (ocl.contains(r)) {
                                Wchain.set(k, r, multiServerTerm + STeff.get(k, r) * (1 + totArvlQlenSeenByOpen.get(r, k)));
                            } else {
                                if ("default".equals(method) || "amva_lin".equals(method) || "lin".equals(method)
                                        || "amva_qdlin".equals(method) || "qdlin".equals(method)) {
                                    double tmp = 0.0;
                                    for (Integer c : ccl) tmp += Nchain_in.get(0, c) * gamma.get(c).get(r, k);
                                    Wchain.set(k, r,
                                            multiServerTerm + STeff.get(k, r) * (1 + interpTotArvlQlen.get(k, 0) + tmp - gamma.get(r).get(r, k)));
                                } else {
                                    Wchain.set(k, r,
                                            multiServerTerm + STeff.get(k, r) * (1 + interpTotArvlQlen.get(k, 0) + (Nt - 1) * gamma.get(0).get(r, k)));
                                }
                            }
                        } else if ("suri".equals(options.config.multiserver)) {
                            if (ocl.contains(r)) {
                                Wchain.set(k, r, STeff.get(k, r) * (1 + totArvlQlenSeenByOpen.get(r, k) * suriFactor.get(k, 0)));
                            } else {
                                if ("default".equals(method) || "amva_lin".equals(method) || "lin".equals(method)
                                        || "amva_qdlin".equals(method) || "qdlin".equals(method)) {
                                    double tmp = 0.0;
                                    for (Integer c : ccl) tmp += Nchain_in.get(0, c) * gamma.get(c).get(r, k);
                                    Wchain.set(k, r,
                                            STeff.get(k, r) * (1 + (interpTotArvlQlen.get(k, 0) + tmp - gamma.get(r).get(r, k)) * suriFactor.get(k, 0)));
                                } else {
                                    Wchain.set(k, r,
                                            STeff.get(k, r) * (1 + (interpTotArvlQlen.get(k, 0) + (Nt - 1) * gamma.get(0).get(r, k)) * suriFactor.get(k, 0)));
                                }
                            }
                        } else {
                            if (ocl.contains(r)) {
                                Wchain.set(k, r, STeff.get(k, r) * (1 + totArvlQlenSeenByOpen.get(r, k)));
                            } else {
                                if ("default".equals(method) || "amva_lin".equals(method) || "lin".equals(method)
                                        || "amva_qdlin".equals(method) || "qdlin".equals(method)) {
                                    double tmp = 0.0;
                                    for (Integer c : ccl) tmp += Nchain_in.get(0, c) * gamma.get(c).get(r, k);
                                    Wchain.set(k, r,
                                            Wchain.get(k, r) + STeff.get(k, r) * (1 + interpTotArvlQlen.get(k, 0) + tmp - gamma.get(r).get(r, k)));
                                } else {
                                    Wchain.set(k, r,
                                            STeff.get(k, r) * (1 + interpTotArvlQlen.get(k, 0) + (Nt - 1) * gamma.get(0).get(r, k)));
                                }
                            }
                        }
                    } else {
                        if ("seidmann".equals(options.config.multiserver)) {
                            double multiServerTerm = STeff.get(k, r) * (nservers.get(k, 0) - 1);
                            if (ocl.contains(r)) {
                                Wchain.set(k, r, multiServerTerm + STeff.get(k, r) * (1 + totArvlQlenSeenByOpen.get(r, k)));
                            } else {
                                if ("default".equals(method) || "amva_lin".equals(method) || "lin".equals(method)
                                        || "amva_qdlin".equals(method) || "qdlin".equals(method)) {
                                    double tmp = 0.0;
                                    for (Integer c : ccl) tmp += Nchain_in.get(0, c) * gamma.get(c).get(r, k);
                                    Wchain.set(k, r,
                                            multiServerTerm + STeff.get(k, r) * (1 + totArvlQlenSeenByClosed.get(k, 0) + tmp - gamma.get(r).get(r, k)));
                                } else {
                                    Wchain.set(k, r,
                                            multiServerTerm + STeff.get(k, r) * (1 + totArvlQlenSeenByClosed.get(k, 0) + (Nt - 1) * gamma.get(0).get(r, k)));
                                }
                            }
                        } else if ("suri".equals(options.config.multiserver)) {
                            if (ocl.contains(r)) {
                                Wchain.set(k, r, STeff.get(k, r) * (1 + totArvlQlenSeenByOpen.get(r, k) * suriFactor.get(k, 0)));
                            } else {
                                Wchain.set(k, r,
                                        STeff.get(k, r) * (1 + (totArvlQlenSeenByClosed.get(k, 0) + (Nt - 1) * gamma.get(0).get(r, k)) * suriFactor.get(k, 0)));
                            }
                        } else {
                            double currentWchain = Wchain.get(k, r);
                            if (ocl.contains(r)) {
                                Wchain.set(k, r, currentWchain + STeff.get(k, r) * (1 + totArvlQlenSeenByOpen.get(r, k)));
                            } else {
                                if ("default".equals(method) || "amva_lin".equals(method) || "lin".equals(method)
                                        || "amva_qdlin".equals(method) || "qdlin".equals(method)) {
                                    double tmp = 0.0;
                                    for (Integer c : ccl) tmp += Nchain_in.get(0, c) * gamma.get(c).get(r, k);
                                    Wchain.set(k, r,
                                            currentWchain + STeff.get(k, r) * (1 + totArvlQlenSeenByClosed.get(k, 0) + tmp - gamma.get(r).get(r, k)));
                                } else {
                                    Wchain.set(k, r,
                                            currentWchain + STeff.get(k, r) * (1 + totArvlQlenSeenByClosed.get(k, 0) + (Nt - 1) * gamma.get(0).get(r, k)));
                                }
                            }
                        }
                    }
                } else if (schedK == SchedStrategy.DPS) {
                    double wDpsBase = 0.0;
                    if (nservers.get(k, 0) > 1) {
                        // see _kb/06-solver-catalog.md for rationale
                        wDpsBase = STeff.get(k, r) * (nservers.get(k, 0) - 1);
                    }

                    Wchain.set(k, r, wDpsBase + STeff.get(k, r) * (1 + selfArvlQlenSeenByClosed.get(k, r)));
                    for (Integer s : sd) {
                        if (schedparam.get(k, s) == schedparam.get(k, r)) {
                            Wchain.set(k, r, Wchain.get(k, r) + STeff.get(k, r) * stationaryQlen.get(k, s));
                        } else if (schedparam.get(k, s) / schedparam.get(k, r) <= Double.POSITIVE_INFINITY) {
                            Wchain.set(k, r, Wchain.get(k, r) + STeff.get(k, r) * stationaryQlen.get(k, s) * schedparam.get(k, s) / schedparam.get(k, r));
                        }
                    }
                } else if (schedK == SchedStrategy.FCFS || schedK == SchedStrategy.SIRO || schedK == SchedStrategy.LCFSPR) {
                    if (STeff.get(k, r) <= 0) continue;

                    Matrix Uchain_r = new Matrix(Uchain_in.getNumRows(), Uchain_in.getNumCols());
                    int ii = 0;
                    while (ii < Uchain_in.getNumRows()) {
                        int j = 0;
                        while (j < Uchain_in.getNumCols()) {
                            Uchain_r.set(ii, j,
                                    (Uchain_in.get(ii, j) / Xchain_in.get(0, j)) * (Xchain_in.get(0, j) + tau.get(r, j)));
                            j++;
                        }
                        ii++;
                    }

                    Matrix Bk = new Matrix(1, K);
                    if (nservers.get(k, 0) > 1) {
                        Matrix deltaclass_r = new Matrix(Xchain_in.getNumRows(), Xchain_in.getNumCols());
                        deltaclass_r.fill(1.0);
                        deltaclass_r.set(0, r, deltaclass.get(0, r));
                        Matrix BK_tmp = new Matrix(1, K);
                        int i2 = 0;
                        while (i2 < K) {
                            BK_tmp.set(0, i2,
                                    deltaclass_r.get(0, i2) * Xchain_in.get(0, i2) * Vchain_in.get(k, i2) * STeff.get(k, i2));
                            i2++;
                        }
                        if (BK_tmp.elementSum() < 0.75) {
                            Bk = BK_tmp.copy();
                        } else {
                            Bk = BK_tmp.elementPower(nservers.get(k, 0) - 1);
                        }
                    } else {
                        Bk.fill(1.0);
                    }


                    if (nservers.get(k, 0) == 1.0 && (((lldscaling != null) && !lldscaling.isEmpty())
                            || ((cdscaling != null) && (cdscaling.size() != 0))
                            || ((jdscaling != null) && (jdscaling.size() != 0)))) {
                        if ("hvmva".equals(options.config.highvar)) {
                            double sum_uchain_r = 0.0;
                            double weightedSum = 0.0;
                            for (Integer s : ccl) {
                                sum_uchain_r += Uchain_r.get(k, s);
                                weightedSum += STeff.get(k, s) * Uchain_r.get(k, s) * (1.0 + SCVchain_in.get(k, s)) / 2.0;
                            }
                            Wchain.set(k, r, weightedSum + STeff.get(k, r) * (1 - sum_uchain_r));
                        } else {
                            Wchain.set(k, r, STeff.get(k, r));
                        }

                        double steff_mult_stationaryQlen = 0.0;
                        for (Integer s : sd) steff_mult_stationaryQlen += STeff.get(k, s) * stationaryQlen.get(k, s);

                        if (ocl.contains(r)) {
                            Wchain.set(k, r,
                                    Wchain.get(k, r) + (STeff.get(k, r) * stationaryQlen.get(k, r) + steff_mult_stationaryQlen));
                        } else {
                            Wchain.set(k, r,
                                    Wchain.get(k, r) + (STeff.get(k, r) * selfArvlQlenSeenByClosed.get(k, r) + steff_mult_stationaryQlen));
                        }
                    } else {
                        double steff_mult_stationarQlen_mult_Bk = 0.0;
                        for (Integer s : sd) steff_mult_stationarQlen_mult_Bk += STeff.get(k, s) * stationaryQlen.get(k, s) * Bk.get(0, s);

                        if ("softmin".equals(options.config.multiserver)) {
                            if (ocl.contains(r)) {
                                Wchain.set(k, r,
                                        STeff.get(k, r) + STeff.get(k, r) * stationaryQlen.get(k, r) * Bk.get(0, r) + steff_mult_stationarQlen_mult_Bk);
                            } else {
                                Wchain.set(k, r,
                                        STeff.get(k, r) + STeff.get(k, r) * selfArvlQlenSeenByClosed.get(k, r) * Bk.get(0, r) + steff_mult_stationarQlen_mult_Bk);
                            }
                        } else if ("suri".equals(options.config.multiserver)) {
                            if (ocl.contains(r)) {
                                double L_m = deltaclass.get(0, r) * stationaryQlen.get(k, r);
                                for (Integer s : sd) L_m += stationaryQlen.get(k, s);
                                Wchain.set(k, r, STeff.get(k, r) + STeff.get(k, r) * L_m * suriFactor.get(k, 0));
                            } else {
                                double L_m = selfArvlQlenSeenByClosed.get(k, r);
                                for (Integer s : sd) L_m += stationaryQlen.get(k, s);
                                Wchain.set(k, r, STeff.get(k, r) + STeff.get(k, r) * L_m * suriFactor.get(k, 0));
                            }
                        } else {
                            if (ocl.contains(r)) {
                                Wchain.set(k, r,
                                        STeff.get(k, r) * (nservers.get(k, 0) - 1) + STeff.get(k, r)
                                                + (STeff.get(k, r) * deltaclass.get(0, r) * stationaryQlen.get(k, r) * Bk.get(0, r) + steff_mult_stationarQlen_mult_Bk));
                            } else {
                                Wchain.set(k, r,
                                        STeff.get(k, r) * (nservers.get(k, 0) - 1) + STeff.get(k, r)
                                                + (STeff.get(k, r) * selfArvlQlenSeenByClosed.get(k, r) * Bk.get(0, r) + steff_mult_stationarQlen_mult_Bk));
                            }
                        }
                    }
                } else if (schedK == SchedStrategy.HOL || schedK == SchedStrategy.FCFSPRIO) {
                    if (STeff.get(k, r) <= 0) continue;

                    Matrix Uchain_r = new Matrix(Uchain_in.getNumRows(), Uchain_in.getNumCols());
                    int ii = 0;
                    while (ii < Uchain_in.getNumRows()) {
                        int j = 0;
                        while (j < Uchain_in.getNumCols()) {
                            Uchain_r.set(ii, j,
                                    (Uchain_in.get(ii, j) / Xchain_in.get(0, j)) * (Xchain_in.get(0, j) + tau.get(r, j)));
                            j++;
                        }
                        ii++;
                    }

                    double prioScaling = 0.0;
                    if ("default".equals(options.config.np_priority) || "cl".equals(options.config.np_priority)) {
                        double UHigherPrio = 0.0;
                        for (Integer h : nnzclasses_hprio.get(r)) {
                            UHigherPrio += Vchain_in.get(k, h) * STeff.get(k, h) * (Xchain_in.get(0, h) - Qchain_in.get(k, h) * tau.get(h));
                        }
                        prioScaling = FastMath.min(Math.max(options.tol, 1 - UHigherPrio), 1 - options.tol);
                    } else if ("shadow".equals(options.config.np_priority)) {
                        double UHigherPrio = 0.0;
                        for (Integer h : nnzclasses_hprio.get(r)) {
                            UHigherPrio += Vchain_in.get(k, h) * STeff.get(k, h) * Xchain_in.get(0, h);
                        }
                        prioScaling = FastMath.min(Math.max(options.tol, 1 - UHigherPrio), 1 - options.tol);
                    }

                    Matrix Bk = new Matrix(1, K);
                    if (nservers.get(k, 0) > 1) {
                        Matrix BK_tmp = new Matrix(1, K);
                        int i2 = 0;
                        while (i2 < K) {
                            BK_tmp.set(0, i2,
                                    deltaclass.get(0, i2) * Xchain_in.get(0, i2) * Vchain_in.get(k, i2) * STeff.get(k, i2));
                            i2++;
                        }

                        if (BK_tmp.elementSum() < 0.75) {
                            if ("softmin".equals(options.config.multiserver)) {
                                Bk = BK_tmp.copy();
                            } else if ("default".equals(options.config.multiserver) || "seidmann".equals(options.config.multiserver)) {
                                BK_tmp.divide(nservers.get(k, 0), Bk, true);
                            } else if ("suri".equals(options.config.multiserver)) {
                                Bk.fill(1.0);
                            }
                        } else {
                            if ("softmin".equals(options.config.multiserver)) {
                                Bk = BK_tmp.elementPower(nservers.get(k, 0));
                            } else if ("default".equals(options.config.multiserver) || "seidmann".equals(options.config.multiserver)) {
                                BK_tmp.divide(nservers.get(k, 0), BK_tmp, true);
                                Bk = BK_tmp.elementPower(nservers.get(k, 0));
                            } else if ("suri".equals(options.config.multiserver)) {
                                Bk.fill(1.0);
                            }
                        }
                    } else {
                        Bk.fill(1.0);
                    }

                    if (nservers.get(k, 0) == 1.0 && (((lldscaling != null) && !lldscaling.isEmpty())
                            || ((cdscaling != null) && (cdscaling.size() != 0))
                            || ((jdscaling != null) && (jdscaling.size() != 0)))) {
                        if ("hvmva".equals(options.config.highvar)) {
                            double sum_uchain_r = 0.0;
                            for (Integer s : ccl) sum_uchain_r += Uchain_r.get(k, s);
                            Wchain.set(k, r, (STeff.get(k, r) / prioScaling) * (1 - sum_uchain_r));
                            for (Integer s : ccl) {
                                double UHigherPrio_s = 0.0;
                                for (Integer h : nnzclasses_hprio.get(s)) {
                                    UHigherPrio_s += Vchain_in.get(k, h) * STeff.get(k, h) * (Xchain_in.get(0, h) - Qchain_in.get(k, h) * tau.get(h));
                                }
                                double prioScaling_s = FastMath.min(Math.max(options.tol, 1 - UHigherPrio_s), 1 - options.tol);
                                Wchain.set(k, r,
                                        Wchain.get(k, r) + (STeff.get(k, s) / prioScaling_s) * Uchain_r.get(k, s) * (1 + SCVchain_in.get(k, s)) / 2);
                            }
                        } else {
                            Wchain.set(k, r, STeff.get(k, r) / prioScaling);
                        }

                        if (ocl.contains(r)) {
                            Wchain.set(k, r,
                                    Wchain.get(k, r) + (STeff.get(k, r) * stationaryQlen.get(k, r)) / prioScaling);
                        } else {
                            Wchain.set(k, r,
                                    Wchain.get(k, r) + (STeff.get(k, r) * selfArvlQlenSeenByClosed.get(k, r)) / prioScaling);
                        }
                    } else {
                        if ("softmin".equals(options.config.multiserver)) {
                            if (ocl.contains(r)) {
                                Wchain.set(k, r,
                                        (STeff.get(k, r) / prioScaling) + STeff.get(k, r) * stationaryQlen.get(k, r) * Bk.get(0, r) / prioScaling);
                            } else {
                                Wchain.set(k, r,
                                        (STeff.get(k, r) / prioScaling) + STeff.get(k, r) * selfArvlQlenSeenByClosed.get(k, r) * Bk.get(0, r) / prioScaling);
                            }
                        } else if ("seidmann".equals(options.config.multiserver) || "default".equals(options.config.multiserver)) {
                            if (ocl.contains(r)) {
                                Wchain.set(k, r,
                                        (STeff.get(k, r) * (nservers.get(k, 0) - 1) / prioScaling) + (STeff.get(k, r) / prioScaling)
                                                + (STeff.get(k, r) * stationaryQlen.get(k, r) * Bk.get(0, r)) / prioScaling);
                            } else {
                                Wchain.set(k, r,
                                        (STeff.get(k, r) * (nservers.get(k, 0) - 1) / prioScaling) + (STeff.get(k, r) / prioScaling)
                                                + (STeff.get(k, r) * selfArvlQlenSeenByClosed.get(k, r) * Bk.get(0, r)) / prioScaling);
                            }
                        } else if ("suri".equals(options.config.multiserver)) {
                            if (ocl.contains(r)) {
                                Wchain.set(k, r,
                                        (STeff.get(k, r) / prioScaling) + (STeff.get(k, r) * stationaryQlen.get(k, r) * suriFactor.get(k, 0)) / prioScaling);
                            } else {
                                Wchain.set(k, r,
                                        (STeff.get(k, r) / prioScaling) + (STeff.get(k, r) * selfArvlQlenSeenByClosed.get(k, r) * suriFactor.get(k, 0)) / prioScaling);
                            }
                        }
                    }
                }
            }
        }
        return new Pair<Matrix, Matrix>(Wchain, STeff);
    }
}
