/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.analyzers;

import org.apache.commons.math3.util.FastMath;

import jline.api.qsys.Qsys_mg1_fb;
import jline.api.qsys.Qsys_mg1_lrpt;
import jline.api.qsys.Qsys_mg1_psjf;
import jline.api.qsys.Qsys_mg1_setf;
import jline.api.qsys.Qsys_mg1_srpt;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.util.matrix.Matrix;

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.line_warning;

/**
 * M/G/1 queueing systems with SIZE-BASED scheduling: SRPT, PSJF, FB (LAS),
 * LRPT and SETF.
 *
 * <p>Port of {@code matlab/src/solvers/MVA/solver_mva_qsys_sizebased_analyzer.m}.
 * The analytical response times come from A. Wierman and M. Harchol-Balter,
 * "Classifying scheduling policies with respect to unfairness in an M/GI/1",
 * SIGMETRICS 2003, evaluated by the {@code qsys_mg1_*} leaves that the JAR
 * already carried; only the analyzer wiring them to a model was missing, so
 * SolverMVA refused every model these disciplines exist for.
 *
 * <p>Applies to a multiclass OPEN Source-Queue-Sink system with Poisson
 * arrivals. The generic AMVA path has no size-based term at all: it would solve
 * the station as if the discipline were size-blind, which is not what SRPT
 * means.
 *
 * <p>NOTE ON THE INDEX SPACE. {@code sn.visits} is indexed by CHAIN, not by
 * station. Reading {@code visits.get(source_ist)} takes the chain whose number
 * happens to equal the Source's station index -- chain 1 -- so on the
 * multiclass models this analyzer exists for, every class beyond the first
 * would get the visit of a chain it does not belong to, which is zero. The
 * chain of each class is looked up explicitly below; the MATLAB reference
 * carries the same correction and the same comment.
 */
public final class Solver_mva_qsys_sizebased_analyzer {

    private Solver_mva_qsys_sizebased_analyzer() {
    }

    /** True for the five disciplines this analyzer serves. */
    public static boolean isSizeBasedPolicy(SchedStrategy s) {
        return s == SchedStrategy.SRPT || s == SchedStrategy.PSJF || s == SchedStrategy.FB
                || s == SchedStrategy.LRPT || s == SchedStrategy.SETF;
    }

    /**
     * @param sn        the network struct
     * @param options   solver options
     * @param schedType the queue's size-based discipline
     */
    public static MVAResult solver_mva_qsys_sizebased_analyzer(NetworkStruct sn,
                                                               SolverOptions options,
                                                               SchedStrategy schedType) {
        long startTime = System.nanoTime();
        int K = sn.nclasses;
        int M = sn.nstations;

        int source_ist = -1;
        int queue_ist = -1;
        int queue_ind = -1;
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Source) {
                source_ist = (int) sn.nodeToStation.get(i);
            } else if (sn.nodetype.get(i) == NodeType.Queue) {
                queue_ist = (int) sn.nodeToStation.get(i);
                queue_ind = i;
            }
        }
        if (source_ist < 0 || queue_ist < 0) {
            line_error("solver_mva_qsys_sizebased_analyzer",
                    "the size-based analyzer needs a Source-Queue-Sink system.");
        }
        int queue_isf = (int) sn.stationToStateful.get(queue_ist);

        // sn.visits is CHAIN-indexed: take the chain of each class
        int[] chainOf = new int[K];
        for (int k = 0; k < K; k++) {
            chainOf[k] = 0;
            for (int c = 0; c < sn.chains.getNumRows(); c++) {
                if (sn.chains.get(c, k) != 0.0) {
                    chainOf[k] = c;
                    break;
                }
            }
        }

        Matrix lambda = new Matrix(1, K);
        Matrix mu = new Matrix(1, K);
        Matrix cs = new Matrix(1, K);
        double[] visits = new double[K];
        for (int k = 0; k < K; k++) {
            visits[k] = sn.visits.get(chainOf[k]).get(queue_isf, k);
            lambda.set(0, k, sn.rates.get(source_ist, k) * visits[k]);
            mu.set(0, k, sn.rates.get(queue_ist, k));
            double scv = sn.scv.get(queue_ist, k);
            cs.set(0, k, (Double.isFinite(scv) && scv > 0) ? FastMath.sqrt(scv) : 1.0);
        }

        for (int k = 0; k < K; k++) {
            if (!(lambda.get(0, k) > 0) || !(mu.get(0, k) > 0)) {
                line_error("solver_mva_qsys_sizebased_analyzer",
                        "invalid arrival or service rates (must be positive).");
            }
        }

        double rho = 0.0;
        for (int k = 0; k < K; k++) {
            rho += lambda.get(0, k) / mu.get(0, k);
        }
        if (rho >= 1.0) {
            line_warning("solver_mva_qsys_sizebased_analyzer",
                    "System is unstable (rho = %.4f >= 1).", Double.valueOf(rho));
        }

        String actualmethod;
        // Ret.qsys_prio carries W in a STATIC field, so it is read immediately
        // after the call and never held across another one.
        if (schedType == SchedStrategy.SRPT) {
            Qsys_mg1_srpt.qsys_mg1_srpt(lambda, mu, cs);
            actualmethod = "mg1.srpt";
        } else if (schedType == SchedStrategy.PSJF) {
            Qsys_mg1_psjf.qsys_mg1_psjf(lambda, mu, cs);
            actualmethod = "mg1.psjf";
        } else if (schedType == SchedStrategy.FB) {
            Qsys_mg1_fb.qsys_mg1_fb(lambda, mu, cs);
            actualmethod = "mg1.fb";
        } else if (schedType == SchedStrategy.LRPT) {
            Qsys_mg1_lrpt.qsys_mg1_lrpt(lambda, mu, cs);
            actualmethod = "mg1.lrpt";
        } else if (schedType == SchedStrategy.SETF) {
            Qsys_mg1_setf.qsys_mg1_setf(lambda, mu, cs);
            actualmethod = "mg1.setf";
        } else {
            actualmethod = null;
            line_error("solver_mva_qsys_sizebased_analyzer",
                    "unsupported scheduling type for the size-based analyzer.");
        }
        Matrix W = Ret.qsys_prio.W;

        Matrix QN = new Matrix(M, K);
        Matrix UN = new Matrix(M, K);
        Matrix RN = new Matrix(M, K);
        Matrix TN = new Matrix(M, K);
        Matrix CN = new Matrix(M, K);
        Matrix AN = new Matrix(M, K);
        Matrix WN = new Matrix(M, K);
        Matrix XN = new Matrix(M, K);
        for (int k = 0; k < K; k++) {
            double lam = lambda.get(0, k);
            double wk = W.get(k);
            RN.set(queue_ist, k, wk * visits[k]);
            CN.set(queue_ist, k, RN.get(queue_ist, k));
            TN.set(source_ist, k, lam);
            TN.set(queue_ist, k, lam);
            AN.set(source_ist, k, lam);
            AN.set(queue_ist, k, lam);
            XN.set(queue_ist, k, lam);
            UN.set(queue_ist, k, lam / mu.get(0, k));
            QN.set(queue_ist, k, lam * wk);
        }
        if (queue_ind < 0) {
            line_error("solver_mva_qsys_sizebased_analyzer", "no Queue node found.");
        }

        MVAResult res = new MVAResult();
        res.QN = QN;
        res.UN = UN;
        res.RN = RN;
        res.TN = TN;
        res.CN = CN;
        res.XN = XN;
        res.AN = AN;
        res.WN = WN;
        res.logNormConstAggr = 0.0;
        res.iter = 1;
        res.method = actualmethod;
        res.runtime = (System.nanoTime() - startTime) / 1000000000.0;
        return res;
    }
}
