/**
 * QRF Analyzer Adapter for SolverCTMC.
 *
 * Port of MATLAB solver_ctmc_qrf_analyzer.m.
 *
 * @since LINE 3.0
 */
package jline.solvers.ctmc.analyzers;

import jline.api.mam.Map_mean;
import jline.api.mapqn.Mapqn_parameters;
import jline.api.mapqn.Mapqn_qr_bounds_bas_parameters;
import jline.api.mapqn.Mapqn_qr_bounds_rsrd_parameters;
import jline.api.mapqn.Mapqn_qrf_noblo_mem;
import jline.api.mapqn.Mapqn_qrf_noblo_mmi;
import jline.api.mapqn.Mapqn_qrf_noblo_mmi_ld;
import jline.api.mapqn.Mapqn_qrf_noblo_mmi_linear;
import jline.api.mapqn.Mapqn_solution;
import jline.api.sn.SnRefreshVisits;
import jline.io.InputOutput;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Solver_ctmc_qrf_analyzer {
    private Solver_ctmc_qrf_analyzer() {}

    public static SolverCTMC.AnalyzerResult solver_ctmc_qrf_analyzer(NetworkStruct sn, SolverOptions options) {
        long T0 = System.nanoTime();

        int M = sn.nstations;
        int K = sn.nclasses;
        int N = (int) sn.njobs.elementSum();
        Matrix S = sn.nservers;

        if (K != 1) {
            InputOutput.line_error("Solver_ctmc_qrf_analyzer",
                    String.format("QRF methods only support single-class networks (found %d classes).", K));
        }
        for (int r = 0; r < K; r++) {
            if (Double.isInfinite(sn.njobs.get(0, r))) {
                InputOutput.line_error("Solver_ctmc_qrf_analyzer", "QRF methods only support closed networks.");
            }
        }

        int[] KPhases = new int[M];
        double[][][][] MAPs = new double[M][][][];
        for (int i = 0; i < M; i++) {
            Station station = sn.stations.get(i);
            JobClass jobClass = sn.jobclasses.get(0);
            MatrixCell procCell = sn.proc.get(station) != null ? sn.proc.get(station).get(jobClass) : null;
            if (procCell != null && procCell.size() > 0) {
                Matrix D0 = procCell.get(0);
                Matrix D1 = procCell.get(1);
                int nPhases = D0.getNumRows();
                KPhases[i] = nPhases;
                double[][] d0Arr = new double[nPhases][nPhases];
                double[][] d1Arr = new double[nPhases][nPhases];
                for (int h = 0; h < nPhases; h++) {
                    for (int kk = 0; kk < nPhases; kk++) {
                        d0Arr[h][kk] = D0.get(h, kk);
                        d1Arr[h][kk] = D1.get(h, kk);
                    }
                }
                MAPs[i] = new double[][][]{d0Arr, d1Arr};
            } else {
                KPhases[i] = 1;
                MAPs[i] = new double[][][]{
                        new double[][]{{-1.0}},
                        new double[][]{{1.0}}
                };
            }
        }

        double[][] rt = new double[M][M];
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                rt[i][j] = sn.rt.get(i, j);
            }
        }

        int Kmax = 0;
        for (int v : KPhases) if (v > Kmax) Kmax = v;
        double[][][] mu = new double[M][Kmax][Kmax];
        double[][][] v = new double[M][Kmax][Kmax];
        for (int i = 0; i < M; i++) {
            double[][][] map = MAPs[i];
            double[][] D0 = map[0];
            double[][] D1 = map[1];
            for (int h = 0; h < KPhases[i]; h++) {
                for (int kk = 0; kk < KPhases[i]; kk++) {
                    // see _kb/06-solver-catalog.md for rationale
                    mu[i][h][kk] = D1[h][kk];
                    v[i][h][kk] = (h == kk) ? 0.0 : D0[h][kk];
                }
            }
        }

        double[] UN_qrf;
        double[] QN_qrf;
        String method = options.method;

        if ("qrf.mmi".equals(method)) {
            Mapqn_solution sol = Mapqn_qrf_noblo_mmi.solve(M, 1, KPhases, N, mu, v, rt);
            UN_qrf = extractUN(sol, M);
            QN_qrf = extractQN(sol, M);
        } else if ("qrf.mem".equals(method)) {
            Mapqn_solution sol = Mapqn_qrf_noblo_mem.solve(MAPs, N, rt);
            UN_qrf = extractUN(sol, M);
            QN_qrf = extractQN(sol, M);
        } else if ("qrf.mmi.ld".equals(method)) {
            double[][] alpha = options.qrfAlpha;
            Mapqn_solution sol = Mapqn_qrf_noblo_mmi_ld.solve(MAPs, N, rt, alpha);
            UN_qrf = extractUN(sol, M);
            QN_qrf = extractQN(sol, M);
        } else if ("qrf.mmi.linear".equals(method)) {
            double[][] alpha = options.qrfAlpha;
            Mapqn_solution sol = Mapqn_qrf_noblo_mmi_linear.solve(MAPs, N, rt, alpha);
            UN_qrf = extractUN(sol, M);
            QN_qrf = extractQN(sol, M);
        } else if ("qrf.bas".equals(method) || "qrf.rsrd".equals(method)) {
            // see _kb/06-solver-catalog.md for rationale
            InputOutput.line_error("Solver_ctmc_qrf_analyzer",
                    String.format("The '%s' method requires blocking parameters (f, MR, BB, MM, "
                            + "MM1, ZZ, ZM), which this implementation cannot yet accept. There is "
                            + "no meaningful default: assuming no blocking (MR=1) yields a bound "
                            + "roughly 31x farther from exact. Use the MATLAB SolverBA with "
                            + "options.config.qrf_params, or choose another qrf.* method.", method));
            UN_qrf = new double[M];
            QN_qrf = new double[M];
        } else {
            InputOutput.line_error("Solver_ctmc_qrf_analyzer",
                    String.format("Unknown QRF method: %s", method));
            UN_qrf = new double[M];
            QN_qrf = new double[M];
        }

        double qnSum = 0.0;
        for (double q : QN_qrf) qnSum += q;
        if (qnSum > 0) {
            for (int i = 0; i < M; i++) {
                QN_qrf[i] = QN_qrf[i] / qnSum * N;
            }
        }

        Matrix QN = new Matrix(M, K);
        Matrix UN = new Matrix(M, K);
        Matrix TN = new Matrix(M, K);
        Matrix RN = new Matrix(M, K);
        Matrix XN = new Matrix(1, K);
        Matrix CN = new Matrix(1, K);

        for (int i = 0; i < M; i++) {
            QN.set(i, 0, QN_qrf[i]);
        }

        Matrix V;
        if (sn.visits != null && !sn.visits.isEmpty()) {
            V = Matrix.cellsum(sn.visits);
        } else {
            NetworkStruct snUpdated = SnRefreshVisits.snRefreshVisits(sn, sn.chains, sn.rt, sn.rtnodes);
            V = Matrix.cellsum(snUpdated.visits);
        }

        int refstat = (int) sn.refstat.get(0, 0);

        // see _kb/06-solver-catalog.md for rationale
        double[] stimes = new double[M];
        for (int i = 0; i < M; i++) {
            MatrixCell procI = sn.proc.get(sn.stations.get(i)) != null
                    ? sn.proc.get(sn.stations.get(i)).get(sn.jobclasses.get(0)) : null;
            if (procI != null && procI.size() > 0) {
                stimes[i] = Map_mean.map_mean(procI.get(0), procI.get(1));
            }
        }

        for (int i = 0; i < M; i++) {
            if (!Double.isInfinite(S.get(i, 0)) && stimes[i] > 0 && V.get(i, 0) > 0
                    && UN_qrf[i] > 0) {
                XN.set(0, 0, UN_qrf[i] * S.get(i, 0) / (V.get(i, 0) * stimes[i]));
                break;
            }
        }
        if (XN.get(0, 0) == 0) {
            // No finite-server station carries load: fall back to the reference
            // station, where QN = XN * V * stime holds exactly for a delay.
            if (stimes[refstat] > 0 && V.get(refstat, 0) > 0) {
                XN.set(0, 0, QN.get(refstat, 0) / (V.get(refstat, 0) * stimes[refstat]));
            }
        }

        for (int i = 0; i < M; i++) {
            double stime = stimes[i];
            if (stime > 0) {
                TN.set(i, 0, XN.get(0, 0) * V.get(i, 0));
                if (Double.isInfinite(S.get(i, 0))) {
                    // Delay (infinite server): UN = QN by LINE convention
                    UN.set(i, 0, QN.get(i, 0));
                } else {
                    // Finite server: the QRF utilization directly
                    UN.set(i, 0, UN_qrf[i]);
                }
                if (TN.get(i, 0) > 0) {
                    RN.set(i, 0, QN.get(i, 0) / TN.get(i, 0));
                }
            }
        }

        if (XN.get(0, 0) > 0) {
            CN.set(0, 0, (double) N / XN.get(0, 0));
        }

        cleanNaN(QN);
        cleanNaN(UN);
        cleanNaN(RN);
        cleanNaN(TN);
        cleanNaN(XN);
        cleanNaN(CN);

        double runtime = (System.nanoTime() - T0) / 1_000_000_000.0;

        return new SolverCTMC.AnalyzerResult(
                QN, UN, RN, TN, CN, XN,
                new Matrix(0, 0),
                new Matrix(0, 0),
                new Matrix(0, 0),
                new MatrixCell(),
                runtime,
                "solver_ctmc_qrf_analyzer",
                sn);
    }

    private static double[] extractUN(Mapqn_solution sol, int M) {
        double[] UN = new double[M];
        for (int i = 0; i < M; i++) {
            UN[i] = sol.getVariable("UN_" + (i + 1));
        }
        return UN;
    }

    private static double[] extractQN(Mapqn_solution sol, int M) {
        double[] QN = new double[M];
        for (int i = 0; i < M; i++) {
            QN[i] = sol.getVariable("QN_" + (i + 1));
        }
        return QN;
    }

    private static void cleanNaN(Matrix m) {
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                if (Double.isNaN(m.get(i, j))) {
                    m.set(i, j, 0.0);
                }
            }
        }
    }

    private static Mapqn_qr_bounds_bas_parameters snToQrfBasParams(
            NetworkStruct sn, int[] KPhases, int N,
            double[][][] mu, double[][][] v,
            double[][] rt, SolverOptions options) {
        int M = sn.nstations;

        Matrix[] muMat = new Matrix[M];
        Matrix[] vMat = new Matrix[M];
        for (int i = 0; i < M; i++) {
            int Ki = KPhases[i];
            Matrix matMu = new Matrix(Ki, Ki);
            Matrix matV = new Matrix(Ki, Ki);
            for (int h = 0; h < Ki; h++) {
                for (int kk = 0; kk < Ki; kk++) {
                    matMu.set(h, kk, mu[i][h][kk]);
                    matV.set(h, kk, v[i][h][kk]);
                }
            }
            muMat[i] = matMu;
            vMat[i] = matV;
        }

        int[] F = new int[M];
        for (int i = 0; i < M; i++) {
            if (sn.cap != null && sn.cap.getNumRows() > 0) {
                int c = (int) sn.cap.get(i, 0);
                F[i] = (c == Integer.MAX_VALUE || Double.isInfinite(sn.cap.get(i, 0))) ? N : c;
            } else {
                F[i] = N;
            }
        }

        Matrix rMat = new Matrix(M, M);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                rMat.set(i, j, rt[i][j]);
            }
        }

        int f = 1;
        int MR = 1;
        Matrix BB = new Matrix(1, M);
        Matrix MM = new Matrix(1, 2);
        Matrix MM1 = new Matrix(1, M);
        int[] ZZ = new int[]{0};

        return new Mapqn_qr_bounds_bas_parameters(M, N, MR, f, KPhases, F, MM, MM1, ZZ, BB, muMat, vMat, rMat);
    }

    private static Mapqn_qr_bounds_rsrd_parameters snToQrfRsrdParams(
            NetworkStruct sn, int[] KPhases, int N,
            double[][][] mu, double[][][] v,
            double[][] rt, SolverOptions options) {
        int M = sn.nstations;

        Matrix[] muMat = new Matrix[M];
        Matrix[] vMat = new Matrix[M];
        for (int i = 0; i < M; i++) {
            int Ki = KPhases[i];
            Matrix matMu = new Matrix(Ki, Ki);
            Matrix matV = new Matrix(Ki, Ki);
            for (int h = 0; h < Ki; h++) {
                for (int kk = 0; kk < Ki; kk++) {
                    matMu.set(h, kk, mu[i][h][kk]);
                    matV.set(h, kk, v[i][h][kk]);
                }
            }
            muMat[i] = matMu;
            vMat[i] = matV;
        }

        int[] F = new int[M];
        for (int i = 0; i < M; i++) {
            if (sn.cap != null && sn.cap.getNumRows() > 0) {
                int c = (int) sn.cap.get(i, 0);
                F[i] = (c == Integer.MAX_VALUE || Double.isInfinite(sn.cap.get(i, 0))) ? N : c;
            } else {
                F[i] = N;
            }
        }

        double[][] alpha;
        if (options.qrfAlpha != null) {
            alpha = new double[M][];
            for (int i = 0; i < M; i++) {
                alpha[i] = options.qrfAlpha[i].clone();
            }
        } else {
            alpha = new double[M][N];
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < N; j++) alpha[i][j] = 1.0;
            }
        }

        Matrix rMat = new Matrix(M, M);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                rMat.set(i, j, rt[i][j]);
            }
        }

        return new Mapqn_qr_bounds_rsrd_parameters(M, N, F, KPhases, muMat, vMat, alpha, rMat);
    }

    private static double[][] deriveQnFromBounds(
            Mapqn_parameters params, int M, int N, Matrix S,
            double[][][][] MAPs, NetworkStruct sn) {
        double[] UBounds = new double[M];
        for (int i = 0; i < M; i++) {
            int queueIdx = i + 1;
            Mapqn_solution minSol;
            Mapqn_solution maxSol;
            if (params instanceof Mapqn_qr_bounds_bas_parameters) {
                minSol = jline.api.mapqn.Mapqn_qr_bounds_bas.solve((Mapqn_qr_bounds_bas_parameters) params, queueIdx, "min");
                maxSol = jline.api.mapqn.Mapqn_qr_bounds_bas.solve((Mapqn_qr_bounds_bas_parameters) params, queueIdx, "max");
            } else {
                minSol = jline.api.mapqn.Mapqn_qr_bounds_rsrd.solve((Mapqn_qr_bounds_rsrd_parameters) params, queueIdx, "min");
                maxSol = jline.api.mapqn.Mapqn_qr_bounds_rsrd.solve((Mapqn_qr_bounds_rsrd_parameters) params, queueIdx, "max");
            }
            UBounds[i] = (minSol.getObjectiveValue() + maxSol.getObjectiveValue()) / 2.0;
        }

        Matrix V;
        if (sn.visits != null && !sn.visits.isEmpty()) {
            V = Matrix.cellsum(sn.visits);
        } else {
            NetworkStruct snUpdated = SnRefreshVisits.snRefreshVisits(sn, sn.chains, sn.rt, sn.rtnodes);
            V = Matrix.cellsum(snUpdated.visits);
        }

        double XNEst = 0.0;
        for (int i = 0; i < M; i++) {
            if (!Double.isInfinite(S.get(i, 0))) {
                MatrixCell procI = sn.proc.get(sn.stations.get(i)) != null
                        ? sn.proc.get(sn.stations.get(i)).get(sn.jobclasses.get(0)) : null;
                if (procI != null && procI.size() > 0 && UBounds[i] > 0 && V.get(i, 0) > 0) {
                    double stime = Map_mean.map_mean(procI.get(0), procI.get(1));
                    if (stime > 0) {
                        XNEst = UBounds[i] * S.get(i, 0) / (V.get(i, 0) * stime);
                        break;
                    }
                }
            }
        }

        double[] QN_qrf = new double[M];
        double[] UN_qrf = UBounds.clone();
        for (int i = 0; i < M; i++) {
            MatrixCell procI = sn.proc.get(sn.stations.get(i)) != null
                    ? sn.proc.get(sn.stations.get(i)).get(sn.jobclasses.get(0)) : null;
            if (procI != null && procI.size() > 0) {
                double stime = Map_mean.map_mean(procI.get(0), procI.get(1));
                double TNi = XNEst * V.get(i, 0);
                if (Double.isInfinite(S.get(i, 0))) {
                    QN_qrf[i] = TNi * stime;
                    UN_qrf[i] = QN_qrf[i];
                } else {
                    if (UBounds[i] < 1) {
                        QN_qrf[i] = TNi * stime / (1 - UBounds[i]);
                    } else {
                        QN_qrf[i] = (double) N;
                    }
                }
            }
        }

        double qnSum = 0.0;
        for (double q : QN_qrf) qnSum += q;
        if (qnSum > 0) {
            for (int i = 0; i < M; i++) {
                QN_qrf[i] = QN_qrf[i] / qnSum * N;
            }
        }

        return new double[][]{UN_qrf, QN_qrf};
    }
}
