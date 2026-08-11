package jline.solvers.mam.handlers;

import java.util.HashMap;
import java.util.Map;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.api.mam.Map_exponential;
import jline.api.mam.Map_mean;
import jline.api.mam.Map_normalize;
import jline.api.mam.Map_pie;
import jline.api.mam.Map_scale;
import jline.api.mam.Mmap_compress;
import jline.api.mam.Mmap_hide;
import jline.api.mam.Mmap_lambda;
import jline.api.mam.Ph_reindex;
import jline.api.mam.Qbd_depproc_etaqa;
import jline.api.mam.Qbd_depproc_etaqa_ps;
import jline.io.InputOutput;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.lib.butools.MMAPPH1FCFS;
import jline.solvers.SolverOptions;
import jline.solvers.mam.MAMResult;
import jline.api.da.Da_fpi;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Solver_mam {
    private Solver_mam() {}

    public static MAMResult solver_mam(NetworkStruct sn, SolverOptions options) {
        long startTime = System.nanoTime();
        SolverOptions.Config config = options.config;
        @SuppressWarnings("rawtypes")
        Map PH = sn.proc;
        int M = sn.nstations;
        int K = sn.nclasses;
        int C = sn.nchains;
        Matrix V = new Matrix(sn.visits.get(0));
        if (sn.nchains > 1) {
            for (int i = 1; i < sn.nchains; i++) {
                // see _kb/06-solver-catalog.md for rationale
                V.addEq(1.0, sn.visits.get(i));
            }
        }

        Matrix QN = new Matrix(M, K, M * K);
        Matrix UN = new Matrix(M, K, M * K);
        Matrix RN = new Matrix(M, K, M * K);
        Matrix TN = new Matrix(M, K, M * K);
        Matrix AN = new Matrix(M, K, M * K);
        Matrix WN = new Matrix(M, K, M * K);
        Matrix CN = new Matrix(1, K, K);
        Matrix XN = new Matrix(1, K, K);

        Matrix lambda = new Matrix(1, K, K);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            Matrix lambdas_inchain = new Matrix(1, inchain.length(), inchain.length());
            for (int i = 0; i < inchain.length(); i++) {
                double rate = sn.rates.get((int) sn.refstat.get((int) inchain.get(0)), (int) inchain.get(i));
                lambdas_inchain.set(0, i, rate);
            }
            double sum = 0.0;
            for (int j = 0; j < lambdas_inchain.length(); j++) {
                if (Double.isFinite(lambdas_inchain.get(j))) {
                    sum = sum + lambdas_inchain.get(j);
                }
            }
            for (int i = 0; i < inchain.length(); i++) {
                lambda.set(0, (int) inchain.get(i), sum);
            }
        }

        Matrix chain = new Matrix(1, K, K);
        for (int k = 0; k < K; k++) {
            for (int i = 0; i < chain.getNumCols(); i++) {
                if (sn.chains.get(i, k) != 0.0) {
                    chain.set(0, k, i);
                    break;
                }
            }
        }

        // Check for non-FCFS queues
        for (int ist = 0; ist < sn.nstations; ist++) {
            SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
            if (sched != SchedStrategy.EXT && sched != SchedStrategy.FCFS && sched != SchedStrategy.HOL
                    && sched != SchedStrategy.FCFSPRIO && sched != SchedStrategy.PS) {
                if (options.verbose != VerboseLevel.SILENT) {
                    InputOutput.line_warning(InputOutput.mfilename(new Object()),
                            "The dec.mmap method does not support non-FCFS queues.");
                }
                MAMResult emptyResult = new MAMResult();
                emptyResult.QN = new Matrix(0, 0);
                emptyResult.UN = new Matrix(0, 0);
                emptyResult.RN = new Matrix(0, 0);
                emptyResult.TN = new Matrix(0, 0);
                emptyResult.CN = new Matrix(0, 0);
                emptyResult.XN = new Matrix(0, 0);
                emptyResult.iter = 0;
                emptyResult.method = "";
                emptyResult.runtime = (System.nanoTime() - startTime) / 1000000000.0;
                return emptyResult;
            }
        }

        boolean isopen = true;
        for (int i = 0; i < sn.njobs.length(); i++) {
            if (Double.isFinite(sn.njobs.get(i))) {
                isopen = false;
                break;
            }
        }

        MAMResult result = new MAMResult();
        if (isopen) {
            int last_it = 0;
            Map<Integer, MatrixCell> pie = new HashMap<Integer, MatrixCell>();
            Map<Integer, MatrixCell> D0 = new HashMap<Integer, MatrixCell>();
            for (int ist = 0; ist < M; ist++) {
                if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.EXT) {
                    for (int i = 0; i < TN.getNumCols(); i++) {
                        if (!Double.isNaN(sn.rates.get(ist, i))) {
                            TN.set(ist, i, sn.rates.get(ist, i));
                        }
                    }
                } else if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.FCFS
                        || sn.sched.get(sn.stations.get(ist)) == SchedStrategy.HOL
                        || sn.sched.get(sn.stations.get(ist)) == SchedStrategy.FCFSPRIO) {
                    for (int k = 0; k < K; k++) {
                        @SuppressWarnings("unchecked")
                        Map<Object, MatrixCell> stationMap = (Map<Object, MatrixCell>) PH.get(sn.stations.get(ist));
                        Matrix D0_ = stationMap.get(sn.jobclasses.get(k)).get(0).copy();
                        D0_.scaleEq(1.0 / sn.nservers.get(ist));
                        Matrix D1 = stationMap.get(sn.jobclasses.get(k)).get(1).copy();
                        D1.scaleEq(1.0 / sn.nservers.get(ist));
                        stationMap.put(sn.jobclasses.get(k),
                                Map_scale.map_scale(stationMap.get(sn.jobclasses.get(k)).get(0),
                                        stationMap.get(sn.jobclasses.get(k)).get(1),
                                        Map_mean.map_mean(D0_, D1)));
                        if (k == 0) {
                            pie.put(Integer.valueOf(ist), new MatrixCell());
                            D0.put(Integer.valueOf(ist), new MatrixCell());
                        }
                        pie.get(Integer.valueOf(ist)).set(k,
                                Map_pie.map_pie(stationMap.get(sn.jobclasses.get(k)).get(0),
                                        stationMap.get(sn.jobclasses.get(k)).get(1)));

                        D0.get(Integer.valueOf(ist)).set(k, stationMap.get(sn.jobclasses.get(k)).get(0));
                        if (D0.get(Integer.valueOf(ist)).get(k).hasNaN()) {
                            stationMap.put(sn.jobclasses.get(k),
                                    Map_exponential.map_exponential(GlobalConstants.Immediate));
                            pie.get(Integer.valueOf(ist)).set(k, Matrix.singleton(1.0));
                            D0.get(Integer.valueOf(ist)).set(k, Matrix.singleton(-GlobalConstants.Immediate));
                        }
                    }
                }
            }
            int it_max = options.iter_max;
            @SuppressWarnings("unchecked")
            final Map<Integer, Map<Integer, MatrixCell>> DEP = (Map<Integer, Map<Integer, MatrixCell>>) (Map) Ph_reindex.ph_reindex(sn);
            // departure-process fixed point (parametric decomposition), driven
            // on the station queue lengths by the generic DA driver
            Da_fpi.Options<Matrix> fpopts = new Da_fpi.Options<Matrix>(
                    it_max, options.iter_tol, new Da_fpi.Norm<Matrix>() {
                @Override
                public double eval(Matrix xnew, Matrix xref) {
                    Matrix diff = xnew.add(-1.0, xref);
                    diff.absEq();
                    return diff.elementDivide(xref).elementMax();
                }
            });
            fpopts.miniter = 4; // legacy tested convergence from it0 >= 3 (4th sweep)
            Da_fpi.Result<Matrix> fpres = Da_fpi.run(new Da_fpi.Sweep<Matrix>() {
                @Override
                public Da_fpi.SweepResult<Matrix> sweep(Matrix xda, int itda) {
                int it = itda - 1;
                if (it == 1) {
                    // see _kb/06-solver-catalog.md for rationale
                    for (int ist = 0; ist < M; ist++) {
                        for (int r = 0; r < K; r++) {
                            // see _kb/06-solver-catalog.md for rationale
                            if (sn.markidx != null && ist < sn.markidx.getNumRows()
                                    && sn.markidx.get(ist, r) > 0
                                    && sn.sched.get(sn.stations.get(ist)) == SchedStrategy.EXT) {
                                continue;
                            }
                            @SuppressWarnings("unchecked")
                            Map<Object, MatrixCell> stationMap = (Map<Object, MatrixCell>) PH.get(sn.stations.get(ist));
                            DEP.get(Integer.valueOf(ist)).put(Integer.valueOf(r),
                                    Map_scale.map_scale(stationMap.get(sn.jobclasses.get(r)).get(0),
                                            stationMap.get(sn.jobclasses.get(r)).get(1),
                                            1 / (lambda.get(r) * V.get(ist, r))));
                        }
                    }
                }

                Map<Integer, MatrixCell> ARV = Solver_mam_traffic.solver_mam_traffic(sn, DEP, config);
                Matrix QN_1 = QN.copy();
                for (int ist = 0; ist < M; ist++) {
                    int ind = (int) sn.stationToNode.get(ist);
                    if (sn.nodetype.get(ind) == NodeType.Queue) {
                        if (ARV.get(Integer.valueOf(ind)).get(0).length() > config.space_max) {
                            System.out.println("Arrival process at node " + ind + " is now at "
                                    + ARV.get(Integer.valueOf(ind)).get(0).length() + " states. Compressing");
                            MatrixCell mmapCell = ARV.get(Integer.valueOf(ind));
                            Matrix[] mmapArray = new Matrix[mmapCell.size()];
                            for (int idx = 0; idx < mmapCell.size(); idx++) {
                                mmapArray[idx] = mmapCell.get(idx);
                            }
                            Matrix[] compressedArray = Mmap_compress.mmap_compress(mmapArray);
                            MatrixCell compressedCell = new MatrixCell(compressedArray.length);
                            for (int i = 0; i < compressedArray.length; i++) {
                                compressedCell.set(i, compressedArray[i]);
                            }
                            ARV.put(Integer.valueOf(ind), compressedCell);
                        }
                        SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
                        if (sched == SchedStrategy.FCFS || sched == SchedStrategy.HOL || sched == SchedStrategy.FCFSPRIO) {
                            // see _kb/06-solver-catalog.md for rationale
                            MatrixCell D = new MatrixCell();
                            MatrixCell arvCell = ARV.get(Integer.valueOf(ind));
                            D.set(0, arvCell.get(0));
                            for (int i = 2; i < arvCell.size(); i++) {
                                D.set(i - 1, arvCell.get(i));
                            }
                            Map<Integer, Matrix> pieMap = new HashMap<Integer, Matrix>();
                            for (int kk = 0; kk < pie.get(Integer.valueOf(ist)).size(); kk++) {
                                pieMap.put(Integer.valueOf(kk), pie.get(Integer.valueOf(ist)).get(kk));
                            }
                            Map<Integer, Matrix> D0Map = new HashMap<Integer, Matrix>();
                            for (int kk = 0; kk < D0.get(Integer.valueOf(ist)).size(); kk++) {
                                D0Map.put(Integer.valueOf(kk), D0.get(Integer.valueOf(ist)).get(kk));
                            }
                            Map<String, Map<Integer, Matrix>> mmapResult = MMAPPH1FCFS.MMAPPH1FCFS(D,
                                    pieMap, D0Map, Integer.valueOf(1), Integer.valueOf(2),
                                    null, null, false, false, null, null);
                            Map<Integer, Matrix> Qret = mmapResult.get("ncMoms");
                            for (int k = 0; k < K; k++) {
                                QN.set(ist, k, Qret.get(Integer.valueOf(k)).elementSum());
                            }
                        } else if (sched == SchedStrategy.PS) {
                            for (int k = 0; k < K; k++) {
                                @SuppressWarnings("unchecked")
                                Map<Object, MatrixCell> stationMap = (Map<Object, MatrixCell>) PH.get(sn.stations.get(ist));
                                UN.set(ist, k, TN.get(ist, k) * Map_mean.map_mean(
                                        stationMap.get(sn.jobclasses.get(k)).get(0),
                                        stationMap.get(sn.jobclasses.get(k)).get(1)));
                            }
                            double Uden = Math.min(1 - GlobalConstants.FineTol, UN.sumRows(ist));
                            for (int k = 0; k < K; k++) {
                                QN.set(ist, k, UN.get(ist, k) / (1 - Uden));
                            }
                        }
                        TN.insertSubMatrix(ist, 0, ist + 1, TN.getNumCols(),
                                Mmap_lambda.mmap_lambda(ARV.get(Integer.valueOf(ind))));
                    }
                    for (int k = 0; k < K; k++) {
                        @SuppressWarnings("unchecked")
                        Map<Object, MatrixCell> stationMap = (Map<Object, MatrixCell>) PH.get(sn.stations.get(ist));
                        UN.set(ist, k, TN.get(ist, k) * Map_mean.map_mean(
                                stationMap.get(sn.jobclasses.get(k)).get(0),
                                stationMap.get(sn.jobclasses.get(k)).get(1)));
                        double meanServiceTime = Map_mean.map_mean(
                                stationMap.get(sn.jobclasses.get(k)).get(0),
                                stationMap.get(sn.jobclasses.get(k)).get(1));
                        QN.set(ist, k, QN.get(ist, k) + TN.get(ist, k) * meanServiceTime
                                * sn.nservers.get(ist) * (sn.nservers.get(ist) - 1) / sn.nservers.get(ist));

                        if (V.get(ist, k) <= GlobalConstants.Zero) {
                            RN.set(ist, k, 0.0);
                        } else {
                            RN.set(ist, k, QN.get(ist, k) / TN.get(ist, k));
                        }
                    }
                }

                for (int ist = 0; ist < M; ist++) {
                    int ind = (int) sn.stationToNode.get(ist);
                    if (sn.nodetype.get(ind) == NodeType.Queue) {
                        for (int r = 0; r < K; r++) {
                            Matrix types = new Matrix(1, K - 1, K - 1);
                            int idx = 0;
                            for (int i = 0; i < K; i++) {
                                if (i != r) {
                                    types.set(idx, (double) i);
                                    idx++;
                                }
                            }
                            MatrixCell A = Mmap_hide.mmap_hide(ARV.get(Integer.valueOf(ind)), types);
                            @SuppressWarnings("unchecked")
                            Map<Object, MatrixCell> stationMap = (Map<Object, MatrixCell>) PH.get(sn.stations.get(ist));
                            MatrixCell S_r = stationMap.get(sn.jobclasses.get(r));
                            int na = A.get(0).getNumRows();
                            int ns = S_r.get(0).getNumRows();
                            int etaqa_n = 8;
                            Object etaqaCfg = config.get("etaqa_trunc");
                            if (etaqaCfg instanceof Integer) {
                                etaqa_n = ((Integer) etaqaCfg).intValue();
                            }
                            int etaqa_sz = (etaqa_n + 1) * na * ns;
                            double rho = UN.sumRows(ist);

                            if (etaqa_sz <= config.space_max && rho < 1 - GlobalConstants.FineTol) {
                                try {
                                    SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
                                    if (sched == SchedStrategy.FCFS || sched == SchedStrategy.HOL
                                            || sched == SchedStrategy.FCFSPRIO) {
                                        MatrixCell dep = Qbd_depproc_etaqa.qbd_depproc_etaqa(A, S_r, etaqa_n);
                                        DEP.get(Integer.valueOf(ist)).put(Integer.valueOf(r),
                                                Map_normalize.map_normalize(dep.get(0), dep.get(1)));
                                    } else if (sched == SchedStrategy.PS) {
                                        MatrixCell dep = Qbd_depproc_etaqa_ps.qbd_depproc_etaqa_ps(A, S_r, etaqa_n);
                                        DEP.get(Integer.valueOf(ist)).put(Integer.valueOf(r),
                                                Map_normalize.map_normalize(dep.get(0), dep.get(1)));
                                    }
                                } catch (Exception e) {
                                    DEP.get(Integer.valueOf(ist)).put(Integer.valueOf(r), new MatrixCell());
                                    DEP.get(Integer.valueOf(ist)).get(Integer.valueOf(r)).set(0,
                                            stationMap.get(sn.jobclasses.get(r)).get(0));
                                    DEP.get(Integer.valueOf(ist)).get(Integer.valueOf(r)).set(1,
                                            stationMap.get(sn.jobclasses.get(r)).get(1));
                                }
                            } else {
                                DEP.get(Integer.valueOf(ist)).put(Integer.valueOf(r), new MatrixCell());
                                DEP.get(Integer.valueOf(ist)).get(Integer.valueOf(r)).set(0,
                                        stationMap.get(sn.jobclasses.get(r)).get(0));
                                DEP.get(Integer.valueOf(ist)).get(Integer.valueOf(r)).set(1,
                                        stationMap.get(sn.jobclasses.get(r)).get(1));
                            }
                            DEP.get(Integer.valueOf(ist)).put(Integer.valueOf(r),
                                    Map_scale.map_scale(
                                            DEP.get(Integer.valueOf(ist)).get(Integer.valueOf(r)).get(0),
                                            DEP.get(Integer.valueOf(ist)).get(Integer.valueOf(r)).get(1),
                                            1 / (lambda.get(r) * V.get(ist, r))));
                        }
                    }
                }
                return new Da_fpi.SweepResult<Matrix>(QN.copy(), QN_1);
                }
            }, QN.copy(), fpopts);
            last_it = fpres.it - 1;
            result.iter = last_it;
            if (options.verbose != VerboseLevel.SILENT) {
                System.out.println("MAM parametric decomposition completed in " + last_it + " iterations");
            }
        } else {
            if (options.verbose != VerboseLevel.SILENT) {
                InputOutput.line_warning(InputOutput.mfilename(new Object()),
                        "This model is not supported by SolverMAM yet. Returning with no result.");
            }
        }
        result.QN = QN;
        result.UN = UN;
        result.RN = RN;
        result.TN = TN;
        result.WN = WN;
        result.AN = AN;
        result.CN = CN;
        result.XN = XN;
        result.runtime = (System.nanoTime() - startTime) / 1000000000.0;
        return result;
    }
}
