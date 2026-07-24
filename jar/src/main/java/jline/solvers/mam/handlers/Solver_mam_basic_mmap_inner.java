/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mam.handlers;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.GlobalConstants;
import jline.api.da.Da_fpi;
import jline.api.mam.Map_acf;
import jline.api.mam.Map_exponential;
import jline.api.mam.Map_mean;
import jline.api.mam.Map_normalize;
import jline.api.mam.Map_pie;
import jline.api.mam.Map_scale;
import jline.api.mam.Mmap_hide;
import jline.api.mam.Mmap_lambda;
import jline.api.mam.Qbd_depproc_etaqa;
import jline.api.mam.Qbd_depproc_etaqa_ps;
import jline.api.qsys.Qsys_mmck;
import jline.api.sn.SnBuildFjSyncMap;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.lib.butools.MMAPPH1FCFS;
import jline.lib.qmam.MAPMAP1Options;
import jline.lib.qmam.MAPMAP1Result;
import jline.lib.qmam.Q_CT_MAP_MAP_1;
import jline.solvers.SolverOptions;
import jline.solvers.mam.MAMResult;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * MAM/MMAP fork-join decomposition algorithm parameterised by per-class arrival
 * rates LAMBDA.
 *
 * <p>Port of matlab/src/solvers/MAM/solver_mam_basic_mmap_inner.m. Performs the
 * departure-process refinement loop only; population enforcement for closed
 * networks is the wrapper's responsibility (see
 * {@link Solver_mam_basic_mmap_closed}).</p>
 */
public final class Solver_mam_basic_mmap_inner {
    private Solver_mam_basic_mmap_inner() {}

    public static MAMResult solver_mam_basic_mmap_inner(final NetworkStruct sn,
                                                        final SolverOptions options,
                                                        final double[] lambda) {
        final SolverOptions.Config config = options.config;
        if (!config.containsKey("fj_sync_q_len")) {
            config.put("fj_sync_q_len", Integer.valueOf(2));
        }
        Object etaqaTrunc = config.get("etaqa_trunc");
        int etaqaTruncVal = (etaqaTrunc instanceof Integer) ? ((Integer) etaqaTrunc).intValue() : 0;
        if (!config.containsKey("etaqa_trunc") || etaqaTruncVal == 0) {
            config.put("etaqa_trunc", Integer.valueOf(8));
            etaqaTruncVal = 8;
        }
        final int etaqaN = etaqaTruncVal;

        final int I = sn.nnodes;
        final int M = sn.nstations;
        final int K = sn.nclasses;
        final Matrix V = Matrix.cellsum(sn.visits);

        final Matrix QN = new Matrix(M, K);
        final Matrix UN = new Matrix(M, K);
        final Matrix RN = new Matrix(M, K);
        final Matrix TN = new Matrix(M, K);
        final Matrix CN = new Matrix(1, K);
        final Matrix XN = new Matrix(1, K);

        final SnBuildFjSyncMap.FjSyncMap fjSyncMap = SnBuildFjSyncMap.sn_build_fj_sync_map(sn);

        // Local PH copy: MATLAB rescales PH{ist}{k} in place; we must not mutate sn.proc
        final Map<Integer, Map<Integer, MatrixCell>> PH = new HashMap<Integer, Map<Integer, MatrixCell>>();
        final Map<Integer, Map<Integer, Matrix>> pie = new HashMap<Integer, Map<Integer, Matrix>>();
        final Map<Integer, Map<Integer, Matrix>> D0 = new HashMap<Integer, Map<Integer, Matrix>>();
        for (int ist = 0; ist < M; ist++) {
            Map<Integer, MatrixCell> row = new HashMap<Integer, MatrixCell>();
            for (int k = 0; k < K; k++) {
                MatrixCell p = null;
                if (sn.proc != null && sn.proc.get(sn.stations.get(ist)) != null) {
                    p = sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(k));
                }
                row.put(k, p);
            }
            PH.put(ist, row);
            pie.put(ist, new HashMap<Integer, Matrix>());
            D0.put(ist, new HashMap<Integer, Matrix>());
        }

        // Prepare PH service distributions
        for (int ist = 0; ist < M; ist++) {
            SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
            if (sched == SchedStrategy.EXT) {
                for (int k = 0; k < K; k++) {
                    double rate = sn.rates.get(ist, k);
                    TN.set(ist, k, Double.isNaN(rate) ? 0.0 : rate);
                }
            } else if (sched == SchedStrategy.FCFS || sched == SchedStrategy.HOL
                    || sched == SchedStrategy.FCFSPRPRIO || sched == SchedStrategy.PS) {
                for (int k = 0; k < K; k++) {
                    MatrixCell p = PH.get(ist).get(k);
                    if (p != null && !p.isEmpty()) {
                        p = Map_scale.map_scale(p, Map_mean.map_mean(p.get(0), p.get(1))
                                / sn.nservers.get(ist));
                        PH.get(ist).put(k, p);
                    }
                    preparePie(PH, pie, D0, ist, k);
                }
            } else if (sched == SchedStrategy.INF) {
                for (int k = 0; k < K; k++) {
                    preparePie(PH, pie, D0, ist, k);
                }
            }
        }

        // departure-process fixed point (FJ parametric decomposition), driven on
        // the station queue lengths by the generic DA driver
        final Map<Integer, Map<Integer, MatrixCell>> DEP = new HashMap<Integer, Map<Integer, MatrixCell>>();

        Da_fpi.Sweep<Matrix> sweep = new Da_fpi.Sweep<Matrix>() {
            @Override
            public Da_fpi.SweepResult<Matrix> sweep(Matrix x, int itnum) {
                if (itnum == 1) {
                    initDep(sn, DEP, PH, V, lambda, I, K);
                }

                Map<Integer, MatrixCell> ARV = Solver_mam_traffic_mmap.solver_mam_traffic_mmap(
                        sn, DEP, config, fjSyncMap);

                Matrix xref = QN.copy();

                for (int ist = 0; ist < M; ist++) {
                    int ind = (int) sn.stationToNode.get(ist);
                    NodeType nt = sn.nodetype.get(ind);
                    SchedStrategy sched = sn.sched.get(sn.stations.get(ist));

                    if (nt == NodeType.Join) {
                        for (int k = 0; k < K; k++) {
                            TN.set(ist, k, lambda[k]);
                            UN.set(ist, k, 0.0);
                            QN.set(ist, k, 0.0);
                            RN.set(ist, k, 0.0);
                        }
                    } else if (nt == NodeType.Queue) {
                        MatrixCell arv = ARV.get(ind);
                        if (arv == null || arv.isEmpty()) {
                            continue;
                        }
                        if (arv.get(0).getNumRows() > config.space_max) {
                            arv = Solver_mam_traffic_mmap.compress(arv, config.compress);
                            ARV.put(ind, arv);
                        }

                        boolean finiteCapUsed = false;
                        if (sched == SchedStrategy.FCFS || sched == SchedStrategy.HOL
                                || sched == SchedStrategy.FCFSPRPRIO) {
                            // see _kb/06-solver-catalog.md for rationale
                            double capVal = sn.cap.get(ist);
                            boolean isFiniteCap = Double.isFinite(capVal)
                                    && capVal < (double) GlobalConstants.MaxInt;
                            if (isFiniteCap) {
                                int capK = (int) capVal;
                                Mam_detect_mmck.Result det = Mam_detect_mmck.mam_detect_mmck(sn, ist, K, arv);
                                double meanQ_fc;
                                double lossProb_fc;
                                if (det.isMmck) {
                                    Matrix lam = Mmap_lambda.mmap_lambda(arv);
                                    double aggrLambda = 0.0;
                                    for (int k = 0; k < lam.length(); k++) {
                                        double v = lam.get(k);
                                        if (!Double.isNaN(v)) aggrLambda += v;
                                    }
                                    Qsys_mmck.Result ex = Qsys_mmck.qsys_mmck(aggrLambda, det.muRate,
                                            (int) sn.nservers.get(ist), capK);
                                    meanQ_fc = ex.meanQueueLength;
                                    lossProb_fc = ex.lossProbability;
                                } else {
                                    Mam_truncate_renorm.Result tr = Mam_truncate_renorm.mam_truncate_renorm(
                                            arrivalMarks(arv, K), pie.get(ist), D0.get(ist), capK);
                                    meanQ_fc = tr.meanQ;
                                    lossProb_fc = tr.lossProb;
                                }
                                Matrix lambdaInflow = Mmap_lambda.mmap_lambda(arv);
                                double sumTN = 0.0;
                                double[] TN_eff = new double[K];
                                for (int k = 0; k < K; k++) {
                                    double v = lambdaInflow.get(k);
                                    if (Double.isNaN(v)) v = 0.0;
                                    TN_eff[k] = v * (1 - lossProb_fc);
                                    sumTN += TN_eff[k];
                                }
                                double[] S_actual = new double[K];
                                for (int k = 0; k < K; k++) {
                                    S_actual[k] = phMean(PH, ist, k) * sn.nservers.get(ist);
                                }
                                double Wq;
                                if (sumTN > 0) {
                                    double acc = 0.0;
                                    for (int k = 0; k < K; k++) {
                                        double t = TN_eff[k] * S_actual[k];
                                        if (!Double.isNaN(t)) acc += t;
                                    }
                                    double Savg_eff = acc / sumTN;
                                    Wq = Math.max(0.0, meanQ_fc / sumTN - Savg_eff);
                                } else {
                                    Wq = 0.0;
                                }
                                for (int k = 0; k < K; k++) {
                                    TN.set(ist, k, TN_eff[k]);
                                    UN.set(ist, k, TN.get(ist, k) * phMean(PH, ist, k));
                                    if (TN.get(ist, k) > 0) {
                                        RN.set(ist, k, Wq + S_actual[k]);
                                        QN.set(ist, k, TN.get(ist, k) * RN.get(ist, k));
                                    } else {
                                        RN.set(ist, k, 0.0);
                                        QN.set(ist, k, 0.0);
                                    }
                                }
                                finiteCapUsed = true;
                            } else {
                                Matrix lam = Mmap_lambda.mmap_lambda(arv);
                                double rho_ist = 0.0;
                                for (int k = 0; k < K; k++) {
                                    double v = lam.get(k) * phMean(PH, ist, k);
                                    if (Double.isNaN(v)) v = 0.0;
                                    // Exclude self-looping classes: pinned by the
                                    // closed wrapper, must not drive FCFS to saturation
                                    if (sn.isslc.get(k) == 0) {
                                        rho_ist += v;
                                    }
                                }
                                if (rho_ist < 1 - GlobalConstants.FineTol) {
                                    // see _kb/06-solver-catalog.md for rationale
                                    boolean useMapMap1 = (K == 1) && (sn.nservers.get(ist) == 1)
                                            && Math.abs(acfLag1(PH, ist, 0)) > GlobalConstants.CoarseTol;
                                    if (useMapMap1) {
                                        Matrix Carv0 = arv.get(0);
                                        Matrix Carv1 = arv.get(2);
                                        MatrixCell srv = PH.get(ist).get(0);
                                        MAPMAP1Result r = Q_CT_MAP_MAP_1.qCtMapMap1(Carv0, Carv1,
                                                srv.get(0), srv.get(1),
                                                new MAPMAP1Options("MaxNumComp", 100000, 0));
                                        Matrix ql = r.getQueueLength();
                                        double acc = 0.0;
                                        for (int n = 0; n < ql.length(); n++) {
                                            acc += n * ql.get(n);
                                        }
                                        QN.set(ist, 0, acc);
                                    } else {
                                        Map<String, Map<Integer, Matrix>> res = MMAPPH1FCFS.MMAPPH1FCFS(
                                                arrivalMarks(arv, K), pie.get(ist), D0.get(ist),
                                                Integer.valueOf(1), Integer.valueOf(2),
                                                null, null, false, false, null, null);
                                        for (int k = 0; k < K; k++) {
                                            double q = 0.0;
                                            if (res.containsKey("ncMoms") && res.get("ncMoms").containsKey(k)) {
                                                q = res.get("ncMoms").get(k).elementSum();
                                            }
                                            QN.set(ist, k, q);
                                        }
                                    }
                                } else {
                                    // Saturation: bound queue lengths so the wrapper's
                                    // bisection can recognise overload without NaNs
                                    for (int k = 0; k < K; k++) {
                                        if (Double.isFinite(sn.njobs.get(k))) {
                                            QN.set(ist, k, sn.njobs.get(k));
                                        } else {
                                            QN.set(ist, k, 1.0 / GlobalConstants.FineTol);
                                        }
                                    }
                                }
                                for (int k = 0; k < K; k++) {
                                    TN.set(ist, k, lam.get(k));
                                }
                            }
                        } else if (sched == SchedStrategy.PS) {
                            Matrix lam = Mmap_lambda.mmap_lambda(arv);
                            for (int k = 0; k < K; k++) {
                                TN.set(ist, k, lam.get(k));
                                UN.set(ist, k, TN.get(ist, k) * (1.0 / sn.rates.get(ist, k)));
                            }
                            // Self-looping classes must not count toward the PS
                            // sharing denominator
                            double sumU = 0.0;
                            for (int k = 0; k < K; k++) {
                                if (sn.isslc.get(k) == 0) {
                                    sumU += UN.get(ist, k);
                                }
                            }
                            double Uden = Math.min(1 - GlobalConstants.FineTol, sumU);
                            for (int k = 0; k < K; k++) {
                                QN.set(ist, k, UN.get(ist, k) / (1 - Uden));
                            }
                        }

                        if (!finiteCapUsed) {
                            double c = sn.nservers.get(ist);
                            for (int k = 0; k < K; k++) {
                                UN.set(ist, k, TN.get(ist, k) * phMean(PH, ist, k));
                                QN.set(ist, k, QN.get(ist, k)
                                        + TN.get(ist, k) * (phMean(PH, ist, k) * c) * (c - 1) / c);
                                RN.set(ist, k, QN.get(ist, k) / TN.get(ist, k));
                            }
                        }
                    } else {
                        if (sched == SchedStrategy.INF) {
                            MatrixCell arv = ARV.get(ind);
                            if (arv != null && !arv.isEmpty()) {
                                Matrix lam = Mmap_lambda.mmap_lambda(arv);
                                for (int k = 0; k < K; k++) {
                                    TN.set(ist, k, lam.get(k));
                                }
                            }
                            for (int k = 0; k < K; k++) {
                                if (TN.get(ist, k) > 0) {
                                    double S = 1.0 / sn.rates.get(ist, k);
                                    UN.set(ist, k, S * TN.get(ist, k));
                                    QN.set(ist, k, TN.get(ist, k) * S);
                                    RN.set(ist, k, S);
                                }
                            }
                        }
                        // SchedStrategy.EXT: Source, TN already set above
                    }
                }

                // Update departure processes
                updateDep(sn, DEP, ARV, PH, UN, TN, V, lambda, config, etaqaN, M, K);

                return new Da_fpi.SweepResult<Matrix>(QN.copy(), xref);
            }
        };

        Da_fpi.Options<Matrix> fpopts = new Da_fpi.Options<Matrix>(options.iter_max, options.iter_tol,
                new Da_fpi.Norm<Matrix>() {
                    @Override
                    public double eval(Matrix xn, Matrix xr) {
                        // MATLAB: max(abs(xn(:)-xr(:))./(xr(:)+GlobalConstants.FineTol))
                        double d = Double.NEGATIVE_INFINITY;
                        for (int i = 0; i < xn.length(); i++) {
                            double v = Math.abs(xn.get(i) - xr.get(i)) / (xr.get(i) + GlobalConstants.FineTol);
                            if (v > d) d = v;
                        }
                        return d;
                    }
                });
        // legacy loop tested convergence only from the third sweep
        fpopts.miniter = 3;

        Da_fpi.Result<Matrix> fpres = Da_fpi.run(sweep, QN.copy(), fpopts);
        int totiter = fpres.it;

        // Join: derive QN/RN from parallel branch means
        for (int joinIdx = 0; joinIdx < I; joinIdx++) {
            if (sn.nodetype.get(joinIdx) != NodeType.Join) {
                continue;
            }
            int joinStat = (int) sn.nodeToStation.get(joinIdx);
            if (joinStat < 0) {
                continue;
            }
            List<Integer> syncGroups = new ArrayList<Integer>();
            for (int j = 0; j < I; j++) {
                int g = (int) fjSyncMap.nodeSync.get(joinIdx, j);
                if (g > 0 && !syncGroups.contains(Integer.valueOf(g))) {
                    syncGroups.add(Integer.valueOf(g));
                }
            }
            java.util.Collections.sort(syncGroups);

            for (int r = 0; r < K; r++) {
                if (TN.get(joinStat, r) <= 0) {
                    continue;
                }
                double syncDelay = 0.0;
                double joinArrivalRate = 0.0;
                for (int gi = 0; gi < syncGroups.size(); gi++) {
                    int gid = syncGroups.get(gi).intValue();
                    List<Double> branchRt = new ArrayList<Double>();
                    double branchTputSum = 0.0;
                    for (int b = 0; b < I; b++) {
                        if ((int) fjSyncMap.nodeSync.get(joinIdx, b) != gid) {
                            continue;
                        }
                        int branchStat = (int) sn.nodeToStation.get(b);
                        if (branchStat < 0 || RN.get(branchStat, r) <= 0) {
                            continue;
                        }
                        branchRt.add(Double.valueOf(RN.get(branchStat, r)));
                        branchTputSum += TN.get(branchStat, r);
                    }
                    if (branchRt.size() < 2) {
                        continue;
                    }
                    double[] lambdai = new double[branchRt.size()];
                    for (int i = 0; i < branchRt.size(); i++) {
                        lambdai[i] = 1.0 / branchRt.get(i).doubleValue();
                    }
                    double maxBranchRt = inclusionExclusionMaxMean(lambdai);
                    double avg = 0.0;
                    for (int i = 0; i < branchRt.size(); i++) {
                        avg += branchRt.get(i).doubleValue();
                    }
                    avg /= branchRt.size();
                    syncDelay += Math.max(maxBranchRt - avg, 0.0);
                    joinArrivalRate += branchTputSum;
                }
                RN.set(joinStat, r, syncDelay);
                QN.set(joinStat, r, joinArrivalRate * syncDelay);
                UN.set(joinStat, r, 0.0);
            }
        }

        for (int r = 0; r < K; r++) {
            double acc = 0.0;
            for (int ist = 0; ist < M; ist++) {
                acc += RN.get(ist, r);
            }
            CN.set(0, r, acc);
        }
        QN.removeNaN();
        RN.removeNaN();
        UN.removeNaN();
        TN.removeNaN();

        MAMResult result = new MAMResult();
        result.QN = QN;
        result.UN = UN;
        result.RN = RN;
        result.TN = TN;
        result.CN = CN;
        result.XN = XN;
        result.iter = totiter;
        return result;
    }

    private static void preparePie(Map<Integer, Map<Integer, MatrixCell>> PH,
                                   Map<Integer, Map<Integer, Matrix>> pie,
                                   Map<Integer, Map<Integer, Matrix>> D0,
                                   int ist, int k) {
        MatrixCell p = PH.get(ist).get(k);
        if (p == null || p.isEmpty()) {
            Matrix d = new Matrix(1, 1, 1);
            d.set(0, 0, -GlobalConstants.Immediate);
            D0.get(ist).put(k, d);
            Matrix pi = new Matrix(1, 1, 1);
            pi.set(0, 0, 1.0);
            pie.get(ist).put(k, pi);
            PH.get(ist).put(k, Map_exponential.map_exponential(GlobalConstants.Immediate));
            return;
        }
        pie.get(ist).put(k, Map_pie.map_pie(p));
        Matrix d0 = p.get(0);
        if (d0.hasNaN()) {
            Matrix d = new Matrix(1, 1, 1);
            d.set(0, 0, -GlobalConstants.Immediate);
            D0.get(ist).put(k, d);
            Matrix pi = new Matrix(1, 1, 1);
            pi.set(0, 0, 1.0);
            pie.get(ist).put(k, pi);
            PH.get(ist).put(k, Map_exponential.map_exponential(GlobalConstants.Immediate));
        } else {
            D0.get(ist).put(k, d0);
        }
    }

    private static double phMean(Map<Integer, Map<Integer, MatrixCell>> PH, int ist, int k) {
        MatrixCell p = PH.get(ist).get(k);
        if (p == null || p.isEmpty()) {
            return 0.0;
        }
        return Map_mean.map_mean(p.get(0), p.get(1));
    }

    private static double acfLag1(Map<Integer, Map<Integer, MatrixCell>> PH, int ist, int k) {
        MatrixCell p = PH.get(ist).get(k);
        if (p == null || p.isEmpty()) {
            return 0.0;
        }
        Matrix lags = new Matrix(1, 1, 1);
        lags.set(0, 0, 1.0);
        Matrix a = Map_acf.map_acf(p, lags);
        if (a == null || a.length() == 0) {
            return 0.0;
        }
        return a.get(0);
    }

    /**
     * MATLAB {ARV{ind}{[1,3:end]}}: D0 plus the per-class marks, dropping the
     * aggregate D1 at position 2.
     */
    private static MatrixCell arrivalMarks(MatrixCell arv, int K) {
        MatrixCell d = new MatrixCell(K + 1);
        d.set(0, arv.get(0));
        for (int k = 0; k < K; k++) {
            d.set(k + 1, arv.get(2 + k));
        }
        return d;
    }

    private static void initDep(NetworkStruct sn,
                                Map<Integer, Map<Integer, MatrixCell>> DEP,
                                Map<Integer, Map<Integer, MatrixCell>> PH,
                                Matrix V, double[] lambda, int I, int K) {
        DEP.clear();
        for (int ind = 0; ind < I; ind++) {
            Map<Integer, MatrixCell> row = new HashMap<Integer, MatrixCell>();
            NodeType nt = sn.nodetype.get(ind);
            boolean isForkJoin = (nt == NodeType.Fork || nt == NodeType.Join);
            if (sn.isstation.get(ind) > 0 && !isForkJoin) {
                int ist = (int) sn.nodeToStation.get(ind);
                for (int r = 0; r < K; r++) {
                    MatrixCell p = PH.get(ist).get(r);
                    if (V.get(ist, r) > 0 && lambda[r] > 0) {
                        row.put(r, Map_scale.map_scale(p, 1.0 / (lambda[r] * V.get(ist, r))));
                    } else if (sn.isslc.get(r) > 0) {
                        // see _kb/06-solver-catalog.md for rationale
                        row.put(r, Map_exponential.map_exponential(1.0 / GlobalConstants.Zero));
                    } else {
                        row.put(r, p);
                    }
                }
            } else {
                for (int r = 0; r < K; r++) {
                    if (lambda[r] > 0) {
                        row.put(r, Map_exponential.map_exponential(1.0 / lambda[r]));
                    } else {
                        row.put(r, Map_exponential.map_exponential(1.0 / GlobalConstants.Immediate));
                    }
                }
            }
            DEP.put(ind, row);
        }
    }

    private static void updateDep(NetworkStruct sn,
                                  Map<Integer, Map<Integer, MatrixCell>> DEP,
                                  Map<Integer, MatrixCell> ARV,
                                  Map<Integer, Map<Integer, MatrixCell>> PH,
                                  Matrix UN, Matrix TN, Matrix V, double[] lambda,
                                  SolverOptions.Config config, int etaqaN, int M, int K) {
        for (int ist = 0; ist < M; ist++) {
            int ind = (int) sn.stationToNode.get(ist);
            NodeType nt = sn.nodetype.get(ind);
            SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
            if (nt == NodeType.Queue) {
                MatrixCell arv = ARV.get(ind);
                if (arv == null || arv.isEmpty()) {
                    continue;
                }
                boolean isFcfsLike = (sched == SchedStrategy.FCFS || sched == SchedStrategy.HOL
                        || sched == SchedStrategy.FCFSPRPRIO);
                boolean isPs = (sched == SchedStrategy.PS);
                if (!isFcfsLike && !isPs) {
                    continue;
                }
                double rho = 0.0;
                for (int k = 0; k < K; k++) {
                    rho += UN.get(ist, k);
                }
                for (int r = 0; r < K; r++) {
                    // see _kb/06-solver-catalog.md for rationale
                    Matrix types = (K > 1) ? new Matrix(1, K - 1) : new Matrix(0, 0);
                    int t = 0;
                    for (int q = 0; q < K; q++) {
                        if (q != r) {
                            types.set(0, t, q);
                            t++;
                        }
                    }
                    MatrixCell A = Mmap_hide.mmap_hide(arv, types);
                    MatrixCell Srv = PH.get(ist).get(r);
                    int na = A.get(0).getNumRows();
                    int ns = Srv.get(0).getNumRows();
                    int etaqa_sz = (etaqaN + 1) * na * ns;

                    if (isFcfsLike) {
                        MatrixCell dep;
                        if (etaqa_sz <= config.space_max && rho < 1 - GlobalConstants.FineTol) {
                            try {
                                dep = Qbd_depproc_etaqa.qbd_depproc_etaqa(A, Srv, etaqaN);
                                dep = Map_normalize.map_normalize(dep);
                            } catch (Exception e) {
                                dep = Srv;
                            }
                        } else {
                            dep = Srv;
                        }
                        if (V.get(ist, r) > 0 && lambda[r] > 0) {
                            dep = Map_scale.map_scale(dep, 1.0 / (lambda[r] * V.get(ist, r)));
                        }
                        DEP.get(ind).put(r, dep);
                    } else {
                        if (V.get(ist, r) > 0 && lambda[r] > 0) {
                            MatrixCell dep;
                            if (etaqa_sz <= config.space_max && rho < 1 - GlobalConstants.FineTol) {
                                try {
                                    dep = Qbd_depproc_etaqa_ps.qbd_depproc_etaqa_ps(A, Srv, etaqaN);
                                    dep = Map_normalize.map_normalize(dep);
                                } catch (Exception e) {
                                    dep = Srv;
                                }
                            } else {
                                dep = Srv;
                            }
                            dep = Map_scale.map_scale(dep, 1.0 / (lambda[r] * V.get(ist, r)));
                            DEP.get(ind).put(r, dep);
                        }
                    }
                }
            } else if (nt == NodeType.Join) {
                for (int r = 0; r < K; r++) {
                    if (TN.get(ist, r) > 0) {
                        DEP.get(ind).put(r, Map_exponential.map_exponential(1.0 / TN.get(ist, r)));
                    }
                }
            }
        }
    }

    private static double inclusionExclusionMaxMean(double[] lambdai) {
        // MATLAB: sum_{pow} (-1)^pow * sum(1./sum(nchoosek(lambdai,pow+1),2))
        double total = 0.0;
        int n = lambdai.length;
        for (int mask = 1; mask < (1 << n); mask++) {
            double subsetSum = 0.0;
            int bits = 0;
            for (int i = 0; i < n; i++) {
                if ((mask & (1 << i)) != 0) {
                    subsetSum += lambdai[i];
                    bits++;
                }
            }
            if (subsetSum <= 0) {
                continue;
            }
            double term = 1.0 / subsetSum;
            total += (bits % 2 == 1) ? term : -term;
        }
        return total;
    }
}
