package jline.solvers.mam.handlers;

import java.util.HashMap;
import java.util.Map;

import jline.GlobalConstants;
import jline.api.mam.Map_exponential;
import jline.api.mam.Map_pie;
import jline.api.mam.Mmap_shorten;
import jline.api.mam.Mmap_super;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.APH;
import jline.lib.butools.MMAPPH1FCFS;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.solvers.mam.MAMResult;
import jline.util.Utils;
import jline.api.da.Da_fpi;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Solver_mna_open {
    private Solver_mna_open() {}

    public static SolverResult solver_mna_open(NetworkStruct sn, SolverOptions options) {
        SolverOptions.Config config = options.config;
        config.space_max = 32;
        Object depScvObj = config.get("dep_scv");
        String depScv = (depScvObj instanceof String) ? (String) depScvObj : "qna";

        int K = sn.nclasses;
        Matrix rt = sn.rt.copy();
        Matrix S = sn.rates.elementPow(-1.0);
        Matrix scv = sn.scv.copy();
        scv.removeNaN();
        Map<jline.lang.nodes.Station, Map<jline.lang.JobClass, MatrixCell>> PH = sn.proc;
        int I = sn.nnodes;
        int M = sn.nstations;
        int C = sn.nchains;
        Matrix V = Matrix.cellsum(sn.visits);
        Matrix Q = new Matrix(M, K, M * K);
        Map<Integer, MatrixCell> pie = new HashMap<Integer, MatrixCell>();
        Map<Integer, MatrixCell> DO = new HashMap<Integer, MatrixCell>();

        Matrix U = new Matrix(M, K, M * K);
        Matrix R = new Matrix(M, K, M * K);
        Matrix T = new Matrix(M, K, M * K);
        Matrix X = new Matrix(1, K, K);
        for (int ist = 0; ist < M; ist++) {
            SchedStrategy schedI = sn.sched.get(sn.stations.get(ist));
            if (schedI == SchedStrategy.FCFS || schedI == SchedStrategy.INF || schedI == SchedStrategy.PS) {
                pie.put(ist, new MatrixCell());
                DO.put(ist, new MatrixCell());
                for (int k = 0; k < K; k++) {
                    MatrixCell phProc = PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k));
                    pie.get(ist).set(k, Map_pie.map_pie(phProc.get(0), phProc.get(1)));
                    Matrix d0Val = phProc.get(0);
                    if (d0Val.hasNaN()) {
                        DO.get(ist).set(k, Matrix.singleton(-GlobalConstants.Immediate));
                        pie.get(ist).set(k, Matrix.singleton(1.0));
                        PH.get(sn.stations.get(ist)).put(sn.jobclasses.get(k), Map_exponential.map_exponential(GlobalConstants.Immediate));
                    } else {
                        DO.get(ist).set(k, d0Val);
                    }
                }
            }
        }

        Matrix lambda = new Matrix(1, C, C);
        int it = 0;
        Matrix a1 = new Matrix(M, K, M * K);
        Matrix a2 = new Matrix(M, K, M * K);

        Matrix d2 = new Matrix(M, 1, M);
        Matrix f2 = new Matrix(M * K, M * K, (int) Math.pow(M * K, 2));
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                if (sn.nodetype.get((int) sn.stationToNode.get(j)) != NodeType.Source) {
                    for (int r = 0; r < K; r++) {
                        for (int s = 0; s < K; s++) {
                            if (rt.get(i * K + r, j * K + s) > 0) f2.set(i * K + r, j * K + s, 1);
                        }
                    }
                }
            }
        }
        MatrixCell lambdas_inchain = new MatrixCell();
        MatrixCell scvs_inchain = new MatrixCell();
        Matrix d2c = new Matrix(1, C, C);
        int last_source_idx = 0;
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            int sourceIdx = (int) sn.refstat.get((int) inchain.get(0));
            last_source_idx = sourceIdx;
            lambdas_inchain.set(c, new Matrix(1, inchain.length(), inchain.length()));
            for (int i = 0; i < inchain.length(); i++) {
                lambdas_inchain.get(c).set(0, i, sn.rates.get(sourceIdx, (int) inchain.get(i)));
            }
            scvs_inchain.set(c, new Matrix(1, inchain.length(), inchain.length()));
            for (int i = 0; i < inchain.length(); i++) {
                scvs_inchain.get(c).set(0, i, scv.get(sourceIdx, (int) inchain.get(i)));
            }
            Matrix lic = lambdas_inchain.get(c).copy();
            lic.removeInfinite();
            lambda.set(c, lic.elementSum());
            d2c.set(c, Qna_superpos.qna_superpos(lambdas_inchain.get(c), scvs_inchain.get(c)));
            for (int i = 0; i < inchain.length(); i++) {
                T.set(sourceIdx, (int) inchain.get(i), lambdas_inchain.get(c).get(i));
            }
        }
        d2.set(last_source_idx,
                Matrix.extractRows(d2c, last_source_idx, last_source_idx + 1, null).mult(lambda.transpose()).get(0)
                        / lambda.elementSum());
        // flow fixed point on the arrival rates a1 and SCVs a2, driven by the
        // generic DA successive-substitution driver
        Da_fpi.Options<Matrix[]> fpopts = new Da_fpi.Options<Matrix[]>(
                options.iter_max + 1, options.iter_tol, new Da_fpi.Norm<Matrix[]>() {
            @Override
            public double eval(Matrix[] xnew, Matrix[] xref) {
                Matrix d1m = xnew[0].add(-1.0, xref[0]);
                d1m.absEq();
                Matrix d2m = xnew[1].add(-1.0, xref[1]);
                d2m.absEq();
                double e1 = d1m.elementMax();
                double e2 = d2m.elementMax();
                if (e1 > options.iter_tol || e2 > options.iter_tol) {
                    return Double.MAX_VALUE; // keep iterating (legacy OR test)
                } else if (Double.isNaN(e1) || Double.isNaN(e2)) {
                    return Double.NaN;
                }
                return Math.max(e1, e2);
            }
        });
        fpopts.nanstop = true; // legacy while-loop exited on NaN convergence measure
        Da_fpi.Result<Matrix[]> fpres = Da_fpi.run(new Da_fpi.Sweep<Matrix[]>() {
            @Override
            public Da_fpi.SweepResult<Matrix[]> sweep(Matrix[] x, int it) {
            Matrix[] xref = new Matrix[]{a1.copy(), a2.copy()};

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
                for (int j = 0; j < K; j++) { a1.set(i, j, 0); a2.set(i, j, 0); }
                double lambda_i = T.sumRows(i);
                for (int j = 0; j < M; j++) {
                    for (int r = 0; r < K; r++) {
                        for (int s = 0; s < K; s++) {
                            a1.set(i, r, a1.get(i, r) + T.get(j, s) * rt.get(j * K + s, i * K + r));
                            a2.set(i, r, a2.get(i, r) + (1 / lambda_i) * f2.get(j * K + s, i * K + r) * T.get(j, s) * rt.get(j * K + s, i * K + r));
                        }
                    }
                }
            }

            for (int ind = 0; ind < I; ind++) {
                if (sn.isstation.get(ind) == 1.0) {
                    int ist = (int) sn.nodeToStation.get(ind);
                    if (sn.nodetype.get(ind) != NodeType.Join) {
                        SchedStrategy schedI = sn.sched.get(sn.stations.get(ist));
                        if (schedI == SchedStrategy.INF) {
                            for (int s = 0; s < K; s++) d2.set(ist, s, a2.get(ist, s));
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
                        } else if (schedI == SchedStrategy.PS) {
                            for (int c = 0; c < C; c++) {
                                Matrix inchain = sn.inchain.get(c);
                                for (int k1 = 0; k1 < inchain.length(); k1++) {
                                    int k = (int) inchain.get(k1);
                                    T.set(ist, k, lambda.get(c) * V.get(ist, k));
                                    U.set(ist, k, S.get(ist, k) * T.get(ist, k));
                                }
                                double Uden = Math.min(1 - GlobalConstants.FineTol, U.sumRows(ist));
                                for (int k1 = 0; k1 < inchain.length(); k1++) {
                                    int k = (int) inchain.get(k1);
                                    Q.set(ist, k, U.get(ist, k) / (1 - Uden));
                                    R.set(ist, k, Q.get(ist, k) / T.get(ist, k));
                                }
                            }
                        } else if (schedI == SchedStrategy.FCFS) {
                            Matrix mu_ist = new Matrix(1, K, K);
                            for (int i = 0; i < K; i++) mu_ist.set(i, sn.rates.get(ist, i));
                            mu_ist.removeNaN();
                            Matrix rho_ist_class = new Matrix(1, K, K);
                            for (int i = 0; i < K; i++) {
                                rho_ist_class.set(i, a1.get(ist, i) / (GlobalConstants.FineTol + sn.rates.get(ist, i)));
                            }
                            rho_ist_class.removeNaN();
                            double lambda_ist = a1.sumRows(ist);
                            int mi = (int) sn.nservers.get(ist);
                            double rho_ist = rho_ist_class.elementSum() / mi;
                            double c2;
                            if (rho_ist < 1 - options.tol) {
                                double mubar = lambda_ist / rho_ist;
                                c2 = -1.0;
                                for (int r = 0; r < K; r++) {
                                    if (mu_ist.get(r) > 0) {
                                        c2 = c2 + a1.get(ist, r) / lambda_ist * Math.pow(mubar / mi / mu_ist.get(r), 2.0) * (scv.get(ist, r) + 1);
                                    }
                                }
                                d2.set(ist, 1 + Math.pow(rho_ist, 2.0) * (c2 - 1) / Math.sqrt((double) mi)
                                        + (1 - Math.pow(rho_ist, 2.0)) * (a2.sumRows(ist) - 1));
                            } else {
                                for (int k = 0; k < K; k++) Q.set(ist, k, sn.njobs.get(k));
                                d2.set(ist, 1.0);
                            }
                            for (int k = 0; k < K; k++) {
                                T.set(ist, k, a1.get(ist, k));
                                U.set(ist, k, T.get(ist, k) * S.get(ist, k) / sn.nservers.get(ist));
                            }
                        }
                    }
                } else {
                    if (sn.nodetype.get(ind) == NodeType.Fork) {
                        throw new IllegalArgumentException("Fork nodes not supported yet by MNA solver");
                    }
                }
            }

            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    if (sn.nodetype.get((int) sn.stationToNode.get(j)) != NodeType.Source) {
                        for (int r = 0; r < K; r++) {
                            for (int s = 0; s < K; s++) {
                                if (rt.get(i * K + r, j * K + s) > 0) {
                                    f2.set(i * K + r, j * K + s, 1 + rt.get(i * K + r, j * K + s) * (d2.get(i) - 1));
                                }
                            }
                        }
                    }
                }
            }
            return new Da_fpi.SweepResult<Matrix[]>(new Matrix[]{a1.copy(), a2.copy()}, xref);
            }
        }, new Matrix[]{a1.copy(), a2.copy()}, fpopts);
        it = fpres.it;

        for (int ind = 0; ind < I; ind++) {
            if (sn.isstation.get(ind) == 1.0) {
                int ist = (int) sn.nodeToStation.get(ind);
                if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.FCFS) {
                    Matrix mu_ist = new Matrix(1, K, K);
                    for (int i = 0; i < K; i++) mu_ist.set(i, sn.rates.get(ist, i));
                    mu_ist.removeNaN();
                    Matrix rho_ist_class = new Matrix(1, K, K);
                    for (int i = 0; i < K; i++) {
                        rho_ist_class.set(i, a1.get(ist, i) / (GlobalConstants.FineTol + sn.rates.get(ist, i)));
                    }
                    rho_ist_class.removeNaN();
                    int mi = (int) sn.nservers.get(ist);
                    double rho_ist = rho_ist_class.elementSum() / mi;
                    if (rho_ist < 1 - options.tol) {
                        MatrixCell arri_class_anx;
                        MatrixCell arri_class = new MatrixCell();
                        MatrixCell arri_node = new MatrixCell();
                        for (int k = 0; k < K; k++) {
                            if (a1.get(ist, k) == 0.0) {
                                arri_class_anx = Map_exponential.map_exponential(Double.POSITIVE_INFINITY);
                            } else {
                                arri_class_anx = APH.fitMeanAndSCV(1 / a1.get(ist, k), a2.get(ist, k)).getProcess();
                            }
                            arri_class.set(0, arri_class_anx.get(0));
                            arri_class.set(1, arri_class_anx.get(1));
                            arri_class.set(2, arri_class_anx.get(1));
                            if (k == 0) {
                                arri_node.set(0, arri_class.get(0));
                                arri_node.set(1, arri_class.get(1));
                                arri_node.set(2, arri_class.get(2));
                            } else {
                                arri_node = Mmap_super.mmap_super(arri_node, arri_class);
                            }
                        }
                        Map<Integer, Matrix> qretMap = new HashMap<Integer, Matrix>();
                        Object resultObj = MMAPPH1FCFS.solve(Mmap_shorten.mmap_shorten(arri_node),
                                pie.get(ist).toMap(), DO.get(ist).toMap(),
                                1, null, null, null, false, false, null, null).get("ncMoms");
                        if (resultObj instanceof Map) {
                            Map<?, ?> raw = (Map<?, ?>) resultObj;
                            for (Map.Entry<?, ?> e : raw.entrySet()) {
                                if (e.getKey() instanceof Integer && e.getValue() instanceof Matrix) {
                                    qretMap.put((Integer) e.getKey(), (Matrix) e.getValue());
                                }
                            }
                        }
                        for (int i = 0; i < qretMap.size(); i++) {
                            Q.set(ist, i, qretMap.get(i).get(0));
                        }
                    } else {
                        for (int k = 0; k < K; k++) Q.set(ist, k, sn.njobs.get(k));
                    }
                    for (int k = 0; k < K; k++) R.set(ist, k, Q.get(ist, k) / T.get(ist, k));
                }
            }
        }

        MAMResult result = new MAMResult();
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
        result.lG = 0.0;
        return result;
    }
}
