package jline.solvers.mam.handlers;

import java.util.HashMap;
import java.util.Map;

import jline.GlobalConstants;
import jline.api.mam.Map_exponential;
import jline.api.mam.Map_pie;
import jline.api.mam.Map_lambda;
import jline.api.mam.Mmap_lambda;
import jline.api.mam.Mmap_shorten;
import jline.api.mam.Mmap_super_safe;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.APH;
import jline.lib.butools.MMAPPH1FCFS;
import jline.solvers.SolverOptions;
import jline.solvers.mam.MAMResult;
import jline.api.da.Da_fpi;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

public final class Solver_mna_closed {
    private Solver_mna_closed() {}

    public static MAMResult solver_mna_closed(NetworkStruct sn, SolverOptions options) {
        SolverOptions.Config config = options.config;
        config.space_max = 16;

        int K = sn.nclasses;
        Matrix rt = sn.rt.copy();
        Matrix S = sn.rates.elementPow(-1.0);
        Matrix scv = sn.scv.copy();
        scv.removeNaN();
        Map<jline.lang.nodes.Station, Map<JobClass, MatrixCell>> PH = sn.proc;
        int I = sn.nnodes;
        int M = sn.nstations;
        int C = sn.nchains;
        Matrix N = sn.njobs.copy();
        Matrix V = Matrix.cellsum(sn.visits);

        Map<Integer, MatrixCell> pie = new HashMap<Integer, MatrixCell>();
        Map<Integer, MatrixCell> D0 = new HashMap<Integer, MatrixCell>();
        Matrix Q = new Matrix(M, K, M * K);
        Matrix U = new Matrix(M, K, M * K);
        Matrix R = new Matrix(M, K, M * K);
        Matrix T = new Matrix(M, K, M * K);
        Matrix X = new Matrix(1, K, K);
        Matrix a1 = new Matrix(M, K, M * K);

        for (int ist = 0; ist < M; ist++) {
            SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
            if (sched == SchedStrategy.FCFS || sched == SchedStrategy.INF || sched == SchedStrategy.PS) {
                pie.put(ist, new MatrixCell());
                D0.put(ist, new MatrixCell());
                for (int kk = 0; kk < K; kk++) {
                    MatrixCell phProc = PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(kk));
                    pie.get(ist).set(kk, Map_pie.map_pie(phProc.get(0), phProc.get(1)));
                    Matrix d0Val = phProc.get(0);
                    if (d0Val.hasNaN()) {
                        D0.get(ist).set(kk, Matrix.singleton(-GlobalConstants.Immediate));
                        pie.get(ist).set(kk, Matrix.singleton(1.0));
                        PH.get(sn.stations.get(ist)).put(sn.jobclasses.get(kk),
                                Map_exponential.map_exponential(GlobalConstants.Immediate));
                    } else {
                        D0.get(ist).set(kk, d0Val);
                    }
                }
            }
        }

        Matrix lambda = new Matrix(1, K, K);
        Matrix QNc = sn.njobs.copy();
        int it_out = 0;
        Matrix lambda_lb = new Matrix(1, K, K);
        Matrix lambda_ub = new Matrix(1, K, K);
        for (int kk = 0; kk < K; kk++) {
            // The reference bounds the throughput by min(rates(nservers<Inf, k)),
            // which ranges over STATIONS; ranging over nodes instead read past the
            // end of nservers as soon as link() inserted a ClassSwitch node.
            for (int i = 0; i < M; i++) {
                if (sn.nservers.get(i) != Double.POSITIVE_INFINITY) {
                    if (lambda_ub.get(kk) == 0.0) {
                        lambda_ub.set(kk, sn.rates.get(i, kk));
                    } else {
                        lambda_ub.set(kk, Math.min(lambda_ub.get(kk), sn.rates.get(i, kk)));
                    }
                }
            }
        }

        // outer bisection on the class throughputs, driven against the closed
        // population target QNc by the generic DA driver
        Da_fpi.Options<Matrix> outopts = new Da_fpi.Options<Matrix>(
                options.iter_max, options.iter_tol, new Da_fpi.Norm<Matrix>() {
            @Override
            public double eval(Matrix xnew, Matrix xref) {
                Matrix diff = xnew.add(-1.0, xref);
                diff.absEq();
                return diff.elementMax();
            }
        });
        outopts.nanstop = true; // legacy while-loop exited on NaN convergence measure
        Da_fpi.Result<Matrix> outres = Da_fpi.run(new Da_fpi.Sweep<Matrix>() {
            @Override
            public Da_fpi.SweepResult<Matrix> sweep(Matrix QNprev, int it_out) {
            if (it_out == 1) {
                for (int kk = 0; kk < K; kk++) {
                    lambda.set(kk, lambda_ub.get(kk));
                }
            } else {
                for (int kk = 0; kk < K; kk++) {
                    if (QNprev.get(kk) < QNc.get(kk)) {
                        lambda_lb.set(kk, lambda.get(kk));
                    } else {
                        lambda_ub.set(kk, lambda.get(kk));
                    }
                    lambda.set(kk, (lambda_ub.get(kk) + lambda_lb.get(kk)) / 2);
                }
            }
            Q.zero();
            U.zero();
            R.zero();
            T.zero();
            X.zero();
            a1.zero();
            final Matrix a2 = new Matrix(M, K, M * K);
            final Matrix d2 = new Matrix(M, 1, M);
            Matrix f2 = new Matrix(M * K, M * K, (int) Math.pow(M * K, 2));
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    if (sn.nodetype.get((int) sn.stationToNode.get(j)) != NodeType.Source) {
                        for (int r = 0; r < K; r++) {
                            for (int s = 0; s < K; s++) {
                                if (rt.get(i * K + r, j * K + s) > 0) {
                                    f2.set(i * K + r, j * K + s, 1);
                                }
                            }
                        }
                    }
                }
            }
            // inner flow fixed point on the arrival rates a1 and SCVs a2
            Da_fpi.Options<Matrix[]> inopts = new Da_fpi.Options<Matrix[]>(
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
            inopts.nanstop = true;
            Da_fpi.run(new Da_fpi.Sweep<Matrix[]>() {
                @Override
                public Da_fpi.SweepResult<Matrix[]> sweep(Matrix[] xin, int it) {
                Matrix[] xref = new Matrix[]{a1.copy(), a2.copy()};

                for (int c = 0; c < C; c++) {
                    Matrix inchain = sn.inchain.get(c);
                    double chainSum = 0.0;
                    for (int kk = 0; kk < inchain.length(); kk++) {
                        chainSum += Q.sumCols().get((int) inchain.get(kk));
                    }
                    if (chainSum > 0) {
                        for (int kk = 0; kk < inchain.length(); kk++) {
                            int classIdx = (int) inchain.get(kk);
                            for (int m = 0; m < M; m++) {
                                Q.set(m, classIdx, sn.njobs.get(c) * Q.get(m, classIdx) / chainSum);
                            }
                        }
                    }
                }

                for (int kk = 0; kk < K; kk++) {
                    if (sn.isslc.get(kk) == 1.0) {
                        for (int m = 0; m < M; m++) {
                            Q.set(m, kk, 0.0);
                        }
                        Q.set((int) sn.refstat.get(kk), kk, sn.njobs.get(kk));
                    }
                }


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
                                        a2.get(i, r) + (1.0 / lambda_i) * f2.get(j * K + s, i * K + r) * T.get(j, s) * rt.get(j * K + s, i * K + r));
                            }
                        }
                    }
                }

                for (int ind = 0; ind < I; ind++) {
                    if (sn.isstation.get(ind) == 1.0) {
                        int ist = (int) sn.nodeToStation.get(ind);
                        if (sn.nodetype.get(ind) == NodeType.Join) {
                            // no-op
                        } else {
                            SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
                            if (sched == SchedStrategy.INF) {
                                // The reference writes the whole row into what it declared as
                                // an (M,1) vector, so the value its splitting step later reads
                                // back as d2(ist) is a2(ist,1), the FIRST class's flow SCV.
                                // Store exactly that: the wider write is out of bounds here
                                // and threw "Outside of matrix bounds" on any multiclass model
                                // with a Delay, which is every closed model the default now
                                // routes to this analyzer.
                                d2.set(ist, a2.get(ist, 0));
                                for (int c = 0; c < C; c++) {
                                    Matrix inchain = sn.inchain.get(c);
                                    for (int k1 = 0; k1 < inchain.length(); k1++) {
                                        int kk = (int) inchain.get(k1);
                                        T.set(ist, kk, a1.get(ist, kk));
                                        U.set(ist, kk, S.get(ist, kk) * T.get(ist, kk));
                                        Q.set(ist, kk, T.get(ist, kk) * S.get(ist, kk) * V.get(ist, kk));
                                        R.set(ist, kk, Q.get(ist, kk) / T.get(ist, kk));
                                    }
                                }
                            } else if (sched == SchedStrategy.PS) {
                                for (int c = 0; c < C; c++) {
                                    Matrix inchain = sn.inchain.get(c);
                                    for (int k1 = 0; k1 < inchain.length(); k1++) {
                                        int kk = (int) inchain.get(k1);
                                        T.set(ist, kk, lambda.get(c) * V.get(ist, kk));
                                        U.set(ist, kk, S.get(ist, kk) * T.get(ist, kk));
                                    }
                                    double Uden = Math.min(1 - GlobalConstants.FineTol, U.sumRows(ist));
                                    for (int k1 = 0; k1 < inchain.length(); k1++) {
                                        int kk = (int) inchain.get(k1);
                                        double Nc = sn.njobs.get(c);
                                        Q.set(ist, kk, (U.get(ist, kk) - Math.pow(U.get(ist, kk), Nc + 1)) / (1 - Uden));
                                        R.set(ist, kk, Q.get(ist, kk) / T.get(ist, kk));
                                    }
                                }
                            } else if (sched == SchedStrategy.FCFS) {
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
                                    for (int kk = 0; kk < K; kk++) {
                                        double mubar = lambda_ist / rho_ist;
                                        c2 = -1.0;
                                        for (int r = 0; r < K; r++) {
                                            if (mu_ist.get(r) > 0) {
                                                c2 = c2 + a1.get(ist, r) / lambda_ist * FastMath.pow(mubar / mi / mu_ist.get(r), 2) * (scv.get(ist, r) + 1);
                                            }
                                        }
                                    }
                                    d2.set(ist,
                                            1 + FastMath.pow(rho_ist, 2) * (c2 - 1) / FastMath.sqrt((double) mi)
                                                    + (1 - FastMath.pow(rho_ist, 2)) * (a2.sumRows(ist) - 1));
                                } else {
                                    for (int kk = 0; kk < K; kk++) {
                                        Q.set(ist, kk, sn.njobs.get(kk));
                                    }
                                    d2.set(ist, 1.0);
                                }
                                for (int kk = 0; kk < K; kk++) {
                                    T.set(ist, kk, a1.get(ist, kk));
                                    U.set(ist, kk, T.get(ist, kk) * S.get(ist, kk) / sn.nservers.get(ist));
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
                                        f2.set(i * K + r, j * K + s,
                                                1 + rt.get(i * K + r, j * K + s) * (d2.get(i) - 1));
                                    }
                                }
                            }
                        }
                    }
                }
                return new Da_fpi.SweepResult<Matrix[]>(new Matrix[]{a1.copy(), a2.copy()}, xref);
                }
            }, new Matrix[]{a1.copy(), a2.copy()}, inopts);

            for (int ind = 0; ind < I; ind++) {
                if (sn.isstation.get(ind) == 1.0) {
                    int ist = (int) sn.nodeToStation.get(ind);
                    if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.FCFS) {
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
                            for (int kk = 0; kk < K; kk++) {
                                if (a1.get(ist, kk) == 0.0) {
                                    arri_class_anx = Map_exponential.map_exponential(Double.POSITIVE_INFINITY);
                                } else {
                                    arri_class_anx = APH.fitMeanAndSCV(1.0 / a1.get(ist, kk), a2.get(ist, kk)).getProcess();
                                }
                                arri_class.set(0, arri_class_anx.get(0));
                                arri_class.set(1, arri_class_anx.get(1));
                                arri_class.set(2, arri_class_anx.get(1));
                                if (kk == 0) {
                                    arri_node.set(0, arri_class.get(0));
                                    arri_node.set(1, arri_class.get(1));
                                    arri_node.set(2, arri_class.get(2));
                                } else {
                                    Map<Integer, MatrixCell> mmapMap = new HashMap<Integer, MatrixCell>();
                                    mmapMap.put(0, arri_node);
                                    mmapMap.put(1, arri_class);
                                    arri_node = Mmap_super_safe.mmap_super_safe(mmapMap, config.space_max, "default");
                                }
                            }
                            Matrix finite_N = N.copy();
                            finite_N.removeInfinite();
                            int maxLevel = (int) finite_N.elementMax() + 1;
                            MatrixCell D = Mmap_shorten.mmap_shorten(arri_node);
                            Map<Integer, Matrix> pdistr = new HashMap<Integer, Matrix>();
                            Map<Integer, Matrix> Qret = new HashMap<Integer, Matrix>();
                            if (Map_lambda.map_lambda(D.get(0), D.get(1)) < GlobalConstants.FineTol) {
                                for (int kk = 0; kk < K; kk++) {
                                    Matrix pdistrK = new Matrix(1, 2, 2);
                                    pdistrK.set(0, 1 - GlobalConstants.FineTol);
                                    pdistrK.set(1, GlobalConstants.FineTol);
                                    pdistr.put(kk, pdistrK);
                                    Qret.put(kk, Matrix.singleton(GlobalConstants.FineTol / sn.rates.get(ist)));
                                }
                            } else {
                                Map<String, Map<Integer, Matrix>> result = MMAPPH1FCFS.MMAPPH1FCFS(
                                        D, pie.get(ist).toMap(), D0.get(ist).toMap(),
                                        null, maxLevel, null, null, false, false, null, null);
                                Map<Integer, Matrix> ncDistr = result != null ? result.get("ncDistr") : null;
                                if (ncDistr != null) pdistr = new HashMap<Integer, Matrix>(ncDistr);
                                for (int kk = 0; kk < K; kk++) {
                                    pdistr.put(kk,
                                            Matrix.extractRows(pdistr.get(kk).transpose(), 0, (int) N.get(kk) + 1, null));
                                    pdistr.get(kk).absEq();
                                    double sum = 0.0;
                                    for (int i = 0; i < pdistr.get(kk).length() - 1; i++) {
                                        sum = sum + pdistr.get(kk).get(i);
                                    }
                                    pdistr.get(kk).set(pdistr.get(kk).length() - 1, Math.abs(1 - sum));
                                    pdistr.get(kk).scaleEq(1.0 / pdistr.get(kk).elementSum());
                                    Matrix aMatrix = new Matrix(1, (int) N.get(kk) + 1, (int) N.get(kk) + 1);
                                    for (int i = 0; i < aMatrix.length(); i++) {
                                        aMatrix.set(i, (double) i);
                                    }
                                    Matrix bMatrix = new Matrix(1, (int) N.get(kk) + 1, (int) N.get(kk) + 1);
                                    for (int i = 0; i < aMatrix.length(); i++) {
                                        bMatrix.set(i, pdistr.get(kk).get(i));
                                    }
                                    Qret.put(kk,
                                            Matrix.singleton(Math.max(0.0, Math.min(N.get(kk), aMatrix.mult(bMatrix.transpose()).get(0)))));
                                }
                            }
                            for (int i = 0; i < Qret.size(); i++) {
                                Q.set(ist, i, Qret.get(i).get(0));
                            }
                        } else {
                            for (int kk = 0; kk < K; kk++) {
                                Q.set(ist, kk, sn.njobs.get(kk));
                            }
                        }
                        for (int kk = 0; kk < K; kk++) {
                            R.set(ist, kk, Q.get(ist, kk) / T.get(ist, kk));
                        }
                    }
                }
            }
            return new Da_fpi.SweepResult<Matrix>(Q.sumCols(), QNc);
            }
        }, new Matrix(1, K, K), outopts);
        it_out = outres.it;

        for (int kk = 0; kk < K; kk++) {
            if (sn.isslc.get(kk) == 1.0) {
                for (int m = 0; m < M; m++) {
                    Q.set(m, kk, 0.0);
                }
                int ist = (int) sn.refstat.get(kk);
                Q.set(ist, kk, sn.njobs.get(kk));
                T.set(ist, kk, sn.njobs.get(kk) * sn.rates.get(ist, kk));
                R.set(ist, kk, Q.get(ist, kk) / T.get(ist, kk));
                U.set(ist, kk, S.get(ist, kk) * T.get(ist, kk));
            }
        }

        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            if (Double.isFinite(sn.njobs.get(c))) {
                double chainSum = 0.0;
                for (int kk = 0; kk < inchain.length(); kk++) {
                    chainSum += Q.sumCols().get((int) inchain.get(kk));
                }
                if (chainSum > 0) {
                    for (int kk = 0; kk < inchain.length(); kk++) {
                        int classIdx = (int) inchain.get(kk);
                        for (int m = 0; m < M; m++) {
                            Q.set(m, classIdx, sn.njobs.get(c) * Q.get(m, classIdx) / chainSum);
                        }
                    }
                }
            }
        }

        for (int ist = 0; ist < M; ist++) {
            if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.INF) {
                for (int kk = 0; kk < K; kk++) {
                    U.set(ist, kk, Q.get(ist, kk));
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
        result.iter = it_out;
        result.lG = 0.0;

        return result;
    }
}
