/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.handlers;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.api.npfqn.Npfqn_traffic_split_rr;
import jline.api.sn.SnIsOpenModel;
import jline.api.sn.SnRtStations;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.solvers.mva.SolverMVA;
import jline.util.Utils;
import jline.api.da.Da_fpi;
import jline.util.matrix.Matrix;

public final class Solver_qna {
    private Solver_qna() {}

    public static MVAResult solver_qna(NetworkStruct sn, SolverOptions options) {
        SolverOptions.Config config = options.config;
        config.space_max = 1;

        int K = sn.nclasses;
        // sn.rt and sn.visits are indexed by stateful node: project them onto stations
        jline.util.Pair<Matrix, Matrix> rtst = SnRtStations.snRtStations(sn);
        Matrix rt = rtst.getLeft();
        Matrix S = sn.rates.elementPow(-1.0);
        Matrix scv = sn.scv.copy();
        scv.removeNaN();

        int I = sn.nnodes;
        int M = sn.nstations;
        int C = sn.nchains;
        Matrix V = rtst.getRight();
        Matrix Q = new Matrix(M, K, M * K);

        Matrix U = new Matrix(M, K, M * K);
        Matrix R = new Matrix(M, K, M * K);
        Matrix T = new Matrix(M, K, M * K);
        Matrix X = new Matrix(1, K, K);

        Matrix lambda = new Matrix(1, C, C);

        // One predicate for the gate and the run: SolverMVA.methodFeatureSet
        // withholds "qna" for a discipline the station update below has no arm for,
        // and the update refuses it by name rather than leaving that station's row
        // of Q, U, R and T at zero and returning the table as a solution.
        String qnaReason = SolverMVA.qnaSchedulingReason(sn);
        if (!qnaReason.isEmpty()) {
            throw new RuntimeException(qnaReason);
        }

        int it = 0;

        for (int i = 0; i < sn.njobs.length(); i++) {
            if (Double.isFinite(sn.njobs.get(i))) {
                throw new RuntimeException("QNA does not support closed classes.");
            }
        }

        Matrix a1 = new Matrix(M, K, M * K);
        Matrix a2 = new Matrix(M, K, M * K);
        Matrix d2 = new Matrix(M, 1, M);
        Matrix f2 = new Matrix(M * K, M * K, (int) Math.pow((double) (M * K), 2.0));
        // deterministic (round-robin) split degrees, k=1 where the split is Markovian
        Matrix kRR = Npfqn_traffic_split_rr.npfqn_traffic_split_rr(sn);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                if (sn.nodetype.get((int) sn.stationToNode.get(j)) != NodeType.Source) {
                    for (int r = 0; r < K; r++) {
                        for (int s = 0; s < K; s++) {
                            if (rt.get(i * K + r, j * K + s) > 0) {
                                f2.set(i * K + r, j * K + s, 1 + rt.get(i * K + r, j * K + s) * (1 - kRR.get(i, r)));
                            }
                        }
                    }
                }
            }
        }
        Map<Integer, Matrix> lambdas_inchain = new HashMap<Integer, Matrix>();
        Map<Integer, Matrix> scvs_inchain = new HashMap<Integer, Matrix>();
        Matrix d2c = new Matrix(1, C, C);
        int last_source_idx = 0;
        if (SnIsOpenModel.snIsOpenModel(sn)) {
            for (int c = 0; c < C; c++) {
                Matrix inchain = sn.inchain.get(c);
                int sourceIdx = (int) sn.refstat.get((int) inchain.get(0));
                last_source_idx = sourceIdx;
                Matrix lic = new Matrix(1, inchain.length(), inchain.length());
                lambdas_inchain.put(Integer.valueOf(c), lic);
                for (int i = 0; i < inchain.length(); i++) {
                    lic.set(0, i, sn.rates.get(sourceIdx, (int) inchain.get(i)));
                }
                Matrix sic = new Matrix(1, inchain.length(), inchain.length());
                scvs_inchain.put(Integer.valueOf(c), sic);
                for (int i = 0; i < inchain.length(); i++) {
                    sic.set(0, i, scv.get(sourceIdx, (int) inchain.get(i)));
                }
                Matrix lambdas_inchain_C = lic.copy();
                lambdas_inchain_C.removeInfinite();
                lambda.set(c, lambdas_inchain_C.elementSum());
                d2c.set(c, qna_superpos(lic, sic));
                boolean openChain = false;
                for (int i = 0; i < inchain.length(); i++) {
                    if (Utils.isInf(sn.njobs.get((int) inchain.get(i)))) {
                        openChain = true;
                    }
                }
                if (openChain) {
                    for (int i = 0; i < inchain.length(); i++) {
                        T.set(sourceIdx, (int) inchain.get(i), lic.get(i));
                    }
                }
            }
            d2.set(last_source_idx,
                    Matrix.extractRows(d2c, last_source_idx, last_source_idx + 1, null).mult(lambda.transpose())
                            .get(0) / lambda.elementSum());
        }
        // flow fixed point on the queue lengths, driven by the generic DA
        // successive-substitution driver
        Da_fpi.Options<Matrix> fpopts = new Da_fpi.Options<Matrix>(
                options.iter_max + 1, options.iter_tol, new Da_fpi.Norm<Matrix>() {
            @Override
            public double eval(Matrix xnew, Matrix xref) {
                Matrix diff = xnew.add(-1.0, xref);
                diff.absEq();
                return diff.elementMax();
            }
        });
        fpopts.nanstop = true; // legacy while-loop exited on NaN convergence measure
        Da_fpi.Result<Matrix> fpres = Da_fpi.run(new Da_fpi.Sweep<Matrix>() {
            @Override
            public Da_fpi.SweepResult<Matrix> sweep(Matrix x, int it) {
            Matrix xref = Q.copy();

            if (it == 1) {
                for (int c = 0; c < C; c++) {
                    Matrix inchain = sn.inchain.get(c);
                    for (int m = 0; m < M; m++) {
                        for (int i = 0; i < inchain.length(); i++) {
                            T.set(m, (int) inchain.get(i), V.get(m, (int) inchain.get(i)) * lambda.get(c));
                        }
                    }
                }
            }
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < K; j++) {
                    a1.set(i, j, 0);
                    a2.set(i, j, 0);
                }
                double lambda_i = T.sumRows(i);
                for (int j = 0; j < M; j++) {
                    for (int r = 0; r < K; r++) {
                        for (int s = 0; s < K; s++) {
                            a1.set(i, r, a1.get(i, r) + T.get(j, s) * rt.get(j * K + s, i * K + r));
                            a2.set(i, r,
                                    a2.get(i, r) + 1 / lambda_i * f2.get(j * K + s, i * K + r) * T.get(j, s)
                                            * rt.get(j * K + s, i * K + r));
                        }
                    }
                }
            }

            for (int ind = 0; ind < I; ind++) {
                if (sn.isstation.get(ind) == 1.0) {
                    int ist = (int) sn.nodeToStation.get(ind);
                    if (sn.nodetype.get(ind) != NodeType.Join) {
                        if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.INF) {
                            // Departure SCV of a delay. solver_qna.m writes
                            // d2(ist,s) = a2(ist,s) over every class, which grows a
                            // vector declared as zeros(M,1) into an M-by-K matrix;
                            // every downstream read is then the SCALAR d2(ist),
                            // which linear indexing resolves to the FIRST class
                            // column. So the value that survives is a2(ist,1), and
                            // the loops over r and s write nothing else that is
                            // ever read. d2 is an (M,1) vector here, as the C++
                            // port also keeps it (solver_qna.h), so the same
                            // two-index write ran off the end of it the moment a
                            // second class existed -- which the fork-join
                            // transform creates on a model written with one.
                            d2.set(ist, a2.get(ist, 0));
                            for (int c = 0; c < C; c++) {
                                Matrix inchain = sn.inchain.get(c);
                                for (int k1 = 0; k1 < inchain.length(); k1++) {
                                    int k = (int) inchain.get(k1);
                                    T.set(ist, k, a1.get(ist, k));
                                    U.set(ist, k, S.get(ist, k) * T.get(ist, k));
                                    Q.set(ist, k, T.get(ist, k) * S.get(ist, k) * V.get(ist, k));
                                    R.set(ist, k, Q.get(ist, k) / T.get(ist, k));
                                }
                            }
                        } else if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.FCFS) {
                            Matrix mu_ist = new Matrix(1, K, K);
                            for (int i = 0; i < K; i++) {
                                mu_ist.set(i, sn.rates.get(ist, i));
                            }
                            mu_ist.removeNaN();
                            Matrix rho_ist_class = new Matrix(1, K, K);
                            for (int i = 0; i < K; i++) {
                                rho_ist_class.set(i, a1.get(ist, i) / (GlobalConstants.FineTol + sn.rates.get(ist, i)));
                            }
                            rho_ist_class.removeNaN();
                            double lambda_ist = a1.sumRows(ist);
                            int mi = (int) sn.nservers.get(ist);
                            double rho_ist = rho_ist_class.elementSum() / mi;
                            double c2 = 0.0;
                            if (rho_ist < 1 - options.tol) {
                                for (int k = 0; k < K; k++) {
                                    double alpha_mi;
                                    if (rho_ist > 0.7) {
                                        alpha_mi = (Math.pow(rho_ist, (double) mi) + rho_ist) / 2;
                                    } else {
                                        alpha_mi = FastMath.pow(rho_ist, (mi + 1) / 2.0);
                                    }
                                    double mubar = lambda_ist / rho_ist;
                                    c2 = -1.0;
                                    for (int r = 0; r < K; r++) {
                                        if (mu_ist.get(r) > 0) {
                                            c2 = c2 + a1.get(ist, r) / lambda_ist
                                                    * Math.pow(mubar / mi / mu_ist.get(r), 2.0)
                                                    * (scv.get(ist, r) + 1);
                                        }
                                    }
                                    double Wiq = (alpha_mi / mubar) * 1 / (1 - rho_ist) * (a2.sumRows(ist) + c2) / 2;
                                    Q.set(ist, k, a1.get(ist, k) / mu_ist.get(k) + a1.get(ist, k) * Wiq);
                                }
                                d2.set(ist,
                                        1 + Math.pow(rho_ist, 2.0) * (c2 - 1) / Math.sqrt((double) mi)
                                                + (1 - Math.pow(rho_ist, 2.0)) * (a2.sumRows(ist) - 1));
                            } else {
                                for (int k = 0; k < K; k++) {
                                    Q.set(ist, k, sn.njobs.get(k));
                                }
                                d2.set(ist, 1.0);
                            }
                            for (int k = 0; k < K; k++) {
                                T.set(ist, k, a1.get(ist, k));
                                U.set(ist, k, T.get(ist, k) * S.get(ist, k) / sn.nservers.get(ist));
                                R.set(ist, k, Q.get(ist, k) / T.get(ist, k));
                            }
                        }
                    }
                } else {
                    if (sn.nodetype.get(ind) == NodeType.Fork) {
                        throw new RuntimeException("Fork nodes not supported yet by QNA solver.");
                    }
                }
            }

            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    if (sn.nodetype.get((int) sn.stationToNode.get(j)) != NodeType.Source) {
                        for (int r = 0; r < K; r++) {
                            for (int s = 0; s < K; s++) {
                                if (rt.get(i * K + r, j * K + s) > 0) {
                                    // k-fold convolution then Bernoulli thinning at q=k*p: C^2 = (q/k)*d2+1-q
                                    f2.set(i * K + r, j * K + s, 1 + rt.get(i * K + r, j * K + s) * (d2.get(i) - kRR.get(i, r)));
                                }
                            }
                        }
                    }
                }
            }
            return new Da_fpi.SweepResult<Matrix>(Q.copy(), xref);
            }
        }, Q.copy(), fpopts);
        it = fpres.it;
        MVAResult result = new MVAResult();
        Matrix CN = R.sumCols();
        Q.absEq();
        Q.removeNaN();
        U.removeNaN();
        R.removeNaN();
        CN.removeNaN();
        X.removeNaN();
        result.XN = X;
        result.QN = Q;
        result.CN = CN;
        result.UN = U;
        result.RN = R;
        result.TN = T;
        result.iter = it;
        result.logNormConstAggr = 0.0;
        return result;
    }

    public static double qna_superpos(Matrix lambda, Matrix a2) {
        return jline.api.da.Da_traffic_superpos.da_traffic_superpos(lambda, a2);
    }
}
