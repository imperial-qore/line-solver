/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mam.handlers;

import java.util.ArrayList;
import java.util.List;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.mam.MAMResult;
import jline.util.matrix.Matrix;

/**
 * Closed-network wrapper around {@link Solver_mam_basic_mmap_inner}.
 *
 * <p>Port of matlab/src/solvers/MAM/solver_mam_basic_mmap_closed.m. Drives a
 * per-class bisection on the surrogate arrival rate LAMBDA so that the inner
 * solver's queue lengths match the closed population sn.njobs. Mirrors the
 * outer-loop structure of solver_mna_closed.</p>
 */
public final class Solver_mam_basic_mmap_closed {
    private Solver_mam_basic_mmap_closed() {}

    public static MAMResult solver_mam_basic_mmap_closed(NetworkStruct sn, SolverOptions options) {
        int K = sn.nclasses;
        int M = sn.nstations;

        // Per-class bisection bounds: upper = slowest non-INF station rate for that class
        List<Integer> nonInfStations = new ArrayList<Integer>();
        List<Integer> infStations = new ArrayList<Integer>();
        for (int i = 0; i < M; i++) {
            if (Double.isFinite(sn.nservers.get(i))) {
                nonInfStations.add(Integer.valueOf(i));
            } else {
                infStations.add(Integer.valueOf(i));
            }
        }

        double[] lambda_lb = new double[K];
        double[] lambda_ub = new double[K];
        for (int k = 0; k < K; k++) {
            List<Double> rates_k = new ArrayList<Double>();
            for (int idx = 0; idx < nonInfStations.size(); idx++) {
                double r = sn.rates.get(nonInfStations.get(idx).intValue(), k);
                if (Double.isFinite(r) && r > 0) {
                    rates_k.add(Double.valueOf(r));
                }
            }
            if (rates_k.isEmpty()) {
                List<Double> rates_inf = new ArrayList<Double>();
                for (int idx = 0; idx < infStations.size(); idx++) {
                    double r = sn.rates.get(infStations.get(idx).intValue(), k);
                    if (Double.isFinite(r) && r > 0) {
                        rates_inf.add(Double.valueOf(r));
                    }
                }
                if (rates_inf.isEmpty()) {
                    lambda_ub[k] = 1.0;
                } else {
                    double mx = rates_inf.get(0).doubleValue();
                    for (int i = 1; i < rates_inf.size(); i++) {
                        mx = Math.max(mx, rates_inf.get(i).doubleValue());
                    }
                    lambda_ub[k] = mx;
                }
            } else {
                double mn = rates_k.get(0).doubleValue();
                for (int i = 1; i < rates_k.size(); i++) {
                    mn = Math.min(mn, rates_k.get(i).doubleValue());
                }
                lambda_ub[k] = mn;
            }
        }

        // open classes contribute 0; only closed populations gate convergence
        double[] QNc = new double[K];
        for (int k = 0; k < K; k++) {
            double n = sn.njobs.get(k);
            QNc[k] = Double.isFinite(n) ? n : 0.0;
        }
        double[] QN_chain = new double[K];

        int it_out = 0;
        double[] lambda = new double[K];
        for (int k = 0; k < K; k++) {
            lambda[k] = lambda_ub[k];
            // Self-looping classes are pinned by the SLC clamp below; they must not
            // contribute a (saturating) surrogate arrival stream to the inner algorithm.
            if (sn.isslc.get(k) > 0) {
                lambda[k] = 0.0;
            }
        }

        // MATLAB assigns a struct copy; mutating iter_max/verbose on the caller's
        // options object would leak the inner budget back out to the dispatcher.
        SolverOptions inner_options = options.copy();
        inner_options.iter_max = Math.max(20, (int) Math.ceil(options.iter_max / 10.0));
        inner_options.verbose = VerboseLevel.SILENT;
        // see _kb/06-solver-catalog.md for rationale
        if (inner_options.config.space_max > 16) {
            inner_options.config.space_max = 16;
        }

        Matrix QN = new Matrix(M, K);
        Matrix UN = new Matrix(M, K);
        Matrix RN = new Matrix(M, K);
        Matrix TN = new Matrix(M, K);
        Matrix CN = new Matrix(1, K);
        Matrix XN = new Matrix(1, K);

        // Last successful inner-algorithm outputs (fallback if final trial diverges)
        Matrix QN_last = QN.copy();
        Matrix UN_last = UN.copy();
        Matrix RN_last = RN.copy();
        Matrix TN_last = TN.copy();
        Matrix CN_last = CN.copy();
        Matrix XN_last = XN.copy();
        boolean have_good = false;
        boolean algorithm_ok = true;

        double bisect_tol = Math.max(options.iter_tol, 1e-3);

        while (maxAbsDiff(QN_chain, QNc) > bisect_tol && it_out < options.iter_max) {
            it_out = it_out + 1;
            if (it_out > 1) {
                boolean bracket_collapsed = true;
                for (int k = 0; k < K; k++) {
                    if (!Double.isFinite(QNc[k]) || QNc[k] == 0 || sn.isslc.get(k) > 0) {
                        continue;
                    }
                    if (QN_chain[k] < QNc[k]) {
                        lambda_lb[k] = lambda[k];
                    } else {
                        lambda_ub[k] = lambda[k];
                    }
                    lambda[k] = 0.5 * (lambda_lb[k] + lambda_ub[k]);
                    // see _kb/06-solver-catalog.md for rationale
                    if ((lambda_ub[k] - lambda_lb[k])
                            > GlobalConstants.FineTol * Math.max(1.0, Math.abs(lambda_ub[k]))) {
                        bracket_collapsed = false;
                    }
                }
                // see _kb/06-solver-catalog.md for rationale
                if (bracket_collapsed) {
                    it_out = it_out - 1;
                    break;
                }
            }

            try {
                MAMResult r = Solver_mam_basic_mmap_inner.solver_mam_basic_mmap_inner(
                        sn, inner_options, lambda);
                QN = r.QN; UN = r.UN; RN = r.RN;
                TN = r.TN; CN = r.CN; XN = r.XN;
                algorithm_ok = true;
            } catch (Exception e) {
                // see _kb/06-solver-catalog.md for rationale
                algorithm_ok = false;
            }

            if (algorithm_ok) {
                // SLC clamp: all jobs at refstat for self-looping classes
                for (int k = 0; k < K; k++) {
                    if (sn.isslc.get(k) > 0) {
                        for (int i = 0; i < M; i++) {
                            QN.set(i, k, 0.0);
                        }
                        QN.set((int) sn.refstat.get(k), k, sn.njobs.get(k));
                    }
                }
                for (int k = 0; k < K; k++) {
                    double acc = 0.0;
                    for (int i = 0; i < M; i++) {
                        acc += QN.get(i, k);
                    }
                    if (Double.isNaN(acc) || Double.isInfinite(acc)) {
                        acc = 1.0 / GlobalConstants.FineTol;
                    }
                    QN_chain[k] = acc;
                }
                QN_last = QN.copy(); UN_last = UN.copy(); RN_last = RN.copy();
                TN_last = TN.copy(); CN_last = CN.copy(); XN_last = XN.copy();
                have_good = true;
            } else {
                for (int k = 0; k < K; k++) {
                    QN_chain[k] = 1.0 / GlobalConstants.FineTol;
                }
            }
        }

        // If the last trial diverged, fall back to the most recent successful one
        if (!algorithm_ok && have_good) {
            QN = QN_last; UN = UN_last; RN = RN_last;
            TN = TN_last; CN = CN_last; XN = XN_last;
        }

        // Final SLC pass: pin throughput/utilisation at refstat (mirrors solver_mna_closed)
        for (int k = 0; k < K; k++) {
            if (sn.isslc.get(k) > 0) {
                for (int i = 0; i < M; i++) {
                    QN.set(i, k, 0.0);
                }
                int ist = (int) sn.refstat.get(k);
                QN.set(ist, k, sn.njobs.get(k));
                TN.set(ist, k, sn.njobs.get(k) * sn.rates.get(ist, k));
                if (TN.get(ist, k) > 0) {
                    RN.set(ist, k, QN.get(ist, k) / TN.get(ist, k));
                } else {
                    RN.set(ist, k, 0.0);
                }
                UN.set(ist, k, (1.0 / sn.rates.get(ist, k)) * TN.get(ist, k));
            }
        }

        // Population redistribution within chain (matches solver_mna_closed)
        for (int c = 0; c < sn.nchains; c++) {
            Matrix inchain = sn.inchain.get(c);
            if (inchain == null || inchain.length() == 0) {
                continue;
            }
            if (Double.isFinite(sn.njobs.get(c))) {
                double sumQ = 0.0;
                for (int i = 0; i < M; i++) {
                    for (int j = 0; j < inchain.length(); j++) {
                        sumQ += QN.get(i, (int) inchain.get(j));
                    }
                }
                if (sumQ > 0) {
                    for (int i = 0; i < M; i++) {
                        for (int j = 0; j < inchain.length(); j++) {
                            int r = (int) inchain.get(j);
                            QN.set(i, r, sn.njobs.get(c) * QN.get(i, r) / sumQ);
                        }
                    }
                }
            }
        }

        // Delay/INF utilisation = mean number of jobs (matches solver_mna_closed)
        for (int ist = 0; ist < M; ist++) {
            if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.INF) {
                for (int k = 0; k < K; k++) {
                    UN.set(ist, k, QN.get(ist, k));
                }
            }
        }

        for (int k = 0; k < K; k++) {
            double acc = 0.0;
            for (int ist = 0; ist < M; ist++) {
                acc += RN.get(ist, k);
            }
            CN.set(0, k, acc);
        }
        QN.removeNaN();
        UN.removeNaN();
        RN.removeNaN();
        TN.removeNaN();
        CN.removeNaN();

        MAMResult result = new MAMResult();
        result.QN = QN;
        result.UN = UN;
        result.RN = RN;
        result.TN = TN;
        result.CN = CN;
        result.XN = XN;
        result.iter = it_out;
        return result;
    }

    private static double maxAbsDiff(double[] a, double[] b) {
        double d = 0.0;
        for (int i = 0; i < a.length; i++) {
            d = Math.max(d, Math.abs(a[i] - b[i]));
        }
        return d;
    }
}
