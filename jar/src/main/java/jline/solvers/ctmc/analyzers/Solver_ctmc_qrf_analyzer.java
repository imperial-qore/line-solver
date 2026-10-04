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
import jline.api.mapqn.Mapqn_qr_bounds_bas;
import jline.api.mapqn.Mapqn_qr_bounds_bas_parameters;
import jline.api.mapqn.Mapqn_qr_bounds_rsrd;
import jline.api.mapqn.Mapqn_qr_bounds_rsrd_parameters;
import jline.api.mapqn.Mapqn_qrf_noblo_bethe;
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
import jline.api.sn.SnToQrfAlpha;
import jline.api.sn.SnToQrfBlocking;
import jline.api.sn.SnToQrfCapacity;
import jline.solvers.QrfParams;
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
        // THE LOAD-DEPENDENT ARMS SERVE DELAY, MULTISERVER AND LOAD-DEPENDENT
        // STATIONS; THE REST STILL CANNOT. alpha(i,n) multiplies every rate out
        // of station i at population n, which IS the rate law of a delay
        // (alpha = n), of a c-server station (alpha = min(n,c)) and of limited
        // load dependence, so 'qrf.mmi.ld' and 'qrf.mmi.linear' answer the
        // model's OWN chain on all three rather than an approximation of it.
        // SnToQrfAlpha derives alpha and owns the one restriction that
        // survives: a station serving several jobs at once must be exponential,
        // because the QRF local state carries one phase per station.
        //
        // Every other arm builds a population-free q, so it models each station
        // as one server, and a c>1 station solved as c=1 is not a bound in
        // either direction (measured +200% at c=3,N=1 and -10% at c=3,N=3).
        // They keep refusing, by naming the two arms that do serve the model.
        // Same gate as the MATLAB, python and C++ twins.
        String methodRaw = options.method == null ? "qrf.mmi" : options.method;
        SnToQrfAlpha.Result alphaRes = SnToQrfAlpha.snToQrfAlpha(sn);
        boolean methodIsLd = "qrf.mmi.ld".equals(methodRaw) || "qrf.mmi.linear".equals(methodRaw);
        if (alphaRes.ld && !methodIsLd) {
            InputOutput.line_error("Solver_ctmc_qrf_analyzer", String.format(
                    "the '%s' method models every station as a single server: its transition "
                            + "rates carry no population index, so it has nowhere to put the rate "
                            + "of a delay, a multiserver or a load-dependent station. Use "
                            + "'qrf.mmi.ld' or 'qrf.mmi.linear', which do.", methodRaw));
        }
        if (!alphaRes.msg.isEmpty()) {
            InputOutput.line_error("Solver_ctmc_qrf_analyzer", String.format(
                    "The '%s' method cannot be applied: %s", methodRaw, alphaRes.msg));
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
        double[] BN_qrf = null;
        String method = options.method;

        if ("qrf.mmi".equals(method)) {
            Mapqn_solution sol = Mapqn_qrf_noblo_mmi.solve(M, 1, KPhases, N, mu, v, rt);
            UN_qrf = extractUN(sol, M);
            QN_qrf = extractQN(sol, M);
        } else if ("qrf.mem".equals(method)) {
            Mapqn_solution sol = Mapqn_qrf_noblo_mem.solve(MAPs, N, rt);
            UN_qrf = extractUN(sol, M);
            QN_qrf = extractQN(sol, M);
        } else if ("qrf.bethe".equals(method)) {
            // Same polytope and same phase-1 start as qrf.mmi; the objective is
            // the tree-reweighted (Bethe) free entropy at the uniform
            // spanning-tree weight lambda = 1/M, the largest uniform weight at
            // which the program is convex.
            Mapqn_solution sol = Mapqn_qrf_noblo_bethe.solve(M, 1, KPhases, N, mu, v, rt);
            UN_qrf = extractUN(sol, M);
            QN_qrf = extractQN(sol, M);
        } else if ("qrf.mmi.ld".equals(method)) {
            // options.qrfAlpha overrides the derivation, as options.qrfParams
            // does for the blocking tables; absent it, alpha comes from the model.
            double[][] alpha = options.qrfAlpha != null ? options.qrfAlpha : alphaRes.alpha;
            Mapqn_solution sol = Mapqn_qrf_noblo_mmi_ld.solve(MAPs, N, rt, alpha);
            UN_qrf = extractUN(sol, M);
            QN_qrf = extractQN(sol, M);
            BN_qrf = extractBN(sol, M);
        } else if ("qrf.mmi.linear".equals(method)) {
            double[][] alpha = options.qrfAlpha != null ? options.qrfAlpha : alphaRes.alpha;
            Mapqn_solution sol = Mapqn_qrf_noblo_mmi_linear.solve(MAPs, N, rt, alpha);
            UN_qrf = extractUN(sol, M);
            QN_qrf = extractQN(sol, M);
            BN_qrf = extractBN(sol, M);
        } else if ("qrf.bas".equals(method) || "qrf.rsrd".equals(method)) {
            // The blocking tables are DERIVED from sn rather than demanded from
            // the caller: the model fixes every one of them. The refusal this
            // replaces was right only while the alternative was to INVENT them
            // -- substituting no blocking (MR=1) measured 4.16667 from exact on
            // sanity_CQN_rm_{fcfs,ps}_1class where real tables sit at 0.133333,
            // i.e. 31x closer. options.qrfParams stays an explicit override.
            //
            // 'qrf.rsrd' needs NO tables at all: Mapqn_qr_bounds_rsrd_parameters
            // has no f/MR/BB/MM/MM1/ZZ field, and its PBB constraint sums over
            // every queue that can be full, so it also admits SEVERAL finite
            // buffers where 'qrf.bas' admits one.
            Mapqn_parameters qparams;
            if ("qrf.bas".equals(method)) {
                qparams = snToQrfBasParams(sn, KPhases, N, mu, v, rt, options);
            } else {
                qparams = snToQrfRsrdParams(sn, KPhases, N, mu, v, rt, options);
            }
            double[][] du = deriveQnFromBounds(qparams, M, N, S, MAPs, sn);
            UN_qrf = du[0];
            QN_qrf = du[1];
        } else if ("qrf.bas.mmi".equals(method) || "qrf.bas.mem".equals(method)
                || "qrf.bas.bethe".equals(method)) {
            // Same BAS polytope as the LP method name, so the same derivation serves
            // it: snToQrfBasParams builds the tables from sn when the caller
            // supplied none, and there is nothing left to refuse here.
            Mapqn_qr_bounds_bas_parameters bp =
                    (Mapqn_qr_bounds_bas_parameters) snToQrfBasParams(sn, KPhases, N, mu, v, rt, options);
            Mapqn_solution sol;
            if ("qrf.bas.mmi".equals(method)) {
                sol = jline.api.mapqn.Mapqn_qrf_bas_nlp.solveMmi(bp);
            } else if ("qrf.bas.bethe".equals(method)) {
                sol = jline.api.mapqn.Mapqn_qrf_bas_nlp.solveBethe(bp);
            } else {
                sol = jline.api.mapqn.Mapqn_qrf_bas_nlp.solveMem(bp);
            }
            UN_qrf = extractUN(sol, M);
            QN_qrf = extractQN(sol, M);
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

        // The alpha-free arms leave BN null because their alpha is identically
        // 1, and there BN = P(n >= 1) = UN_qrf: a single server's departure rate
        // is proportional to the probability that it is busy. Setting it here
        // rather than in each arm keeps the readout below one formula.
        if (BN_qrf == null) {
            BN_qrf = UN_qrf;
        }

        // System throughput from the ALPHA-WEIGHTED marginal mean BN, the mean
        // number of jobs actually in service: E[min(n,c)] at a c-server station,
        // E[n] at a delay, P(n >= 1) at a single server. That is what the
        // departure rate is proportional to, so T_i = BN_i / stime_i holds
        // exactly at the relaxed point and XN = T_i / V_i. The single-server
        // case is the former UN_qrf[i]*S/(V*stime) unchanged, S being 1 and BN
        // being UN_qrf there; a delay no longer needs the refstat fallback,
        // since alpha = n makes BN = E[n] = QN, which is what it computed.
        for (int i = 0; i < M; i++) {
            if (stimes[i] > 0 && V.get(i, 0) > 0 && BN_qrf[i] > 0) {
                XN.set(0, 0, BN_qrf[i] / (V.get(i, 0) * stimes[i]));
                break;
            }
        }
        if (XN.get(0, 0) == 0) {
            // No station carries load: fall back to the reference station,
            // where QN = XN * V * stime holds exactly for a delay.
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
                    // Finite server: the busy fraction of the station's DECLARED
                    // peak capacity, BN being the mean number in service. The
                    // normalizer is nservers times the reachable lld peak,
                    // LINE's one U = T*S/peak convention (see SnToQrfAlpha); at
                    // peak 1 it is UN_qrf[i], what the alpha-free arms report.
                    UN.set(i, 0, BN_qrf[i] / alphaRes.peak[i]);
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

    private static double[] extractBN(Mapqn_solution sol, int M) {
        double[] BN = new double[M];
        for (int i = 0; i < M; i++) {
            BN[i] = sol.getVariable("BN_" + (i + 1));
        }
        return BN;
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

        // F is an OCCUPANCY BOUND, not a declared capacity: SnToQrfCapacity
        // decides binding through SnGetBufferSize, which folds classcap and the
        // reachable population in as raw sn.cap does not.
        SnToQrfCapacity.Result capResult = SnToQrfCapacity.snToQrfCapacity(sn);
        if (!capResult.msg.isEmpty()) {
            InputOutput.line_error("Solver_ctmc_qrf_analyzer",
                    "The QRF blocking bounds cannot be applied: " + capResult.msg);
        }
        int[] F = capResult.F;

        Matrix rMat = new Matrix(M, M);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                rMat.set(i, j, rt[i][j]);
            }
        }

        // Blocking data: options.qrfParams when the caller supplied it, else
        // derived from sn. There is no third branch -- a silent no-blocking
        // default answers a different model, which is what SnToQrfBlocking
        // exists to avoid.
        QrfParams qp = options.qrfParams;
        if (qp == null) {
            SnToQrfBlocking.Result derived =
                    SnToQrfBlocking.snToQrfBlocking(sn, SnToQrfBlocking.DEFAULT_MAXVARS);
            if (!derived.msg.isEmpty()) {
                InputOutput.line_error("Solver_ctmc_qrf_analyzer", String.format(
                        "The 'qrf.bas' method cannot be applied to this model: %s Supply "
                                + "options.qrfParams explicitly to override the derivation.",
                        derived.msg));
            }
            qp = derived.params;
        }
        int f = qp.f;
        int MR = qp.MR;
        Matrix BB = (qp.BB != null) ? qp.BB : new Matrix(1, M);
        Matrix MM = (qp.MM != null) ? qp.MM : new Matrix(1, 2);
        Matrix MM1 = (qp.MM1 != null) ? qp.MM1 : new Matrix(1, M);
        int[] ZZ = (qp.ZZ != null) ? qp.ZZ : new int[]{0};
        if (qp.F != null) {
            F = qp.F;
        }

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

        // F is an OCCUPANCY BOUND, not a declared capacity: SnToQrfCapacity
        // decides binding through SnGetBufferSize, which folds classcap and the
        // reachable population in as raw sn.cap does not.
        SnToQrfCapacity.Result capResult = SnToQrfCapacity.snToQrfCapacity(sn);
        if (!capResult.msg.isEmpty()) {
            InputOutput.line_error("Solver_ctmc_qrf_analyzer",
                    "The QRF blocking bounds cannot be applied: " + capResult.msg);
        }
        int[] F = capResult.F;

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
        // One LP, MAXIMISING the utilization of queue 1, with every station's
        // utilization read off THAT vertex -- the native-Python convention.
        // The relaxation contains the true joint distribution, so the max
        // direction is the one that yields an upper bound; the libraries
        // default to 'U1min' instead, whose vertex has U = 0 at the objective
        // station and reads out as zero throughput. Solving per station and
        // averaging the two senses, as this method used to, returns a point
        // that is neither bound and that no other codebase produces.
        Mapqn_solution sol;
        if (params instanceof Mapqn_qr_bounds_bas_parameters) {
            sol = jline.api.mapqn.Mapqn_qr_bounds_bas.solve((Mapqn_qr_bounds_bas_parameters) params, 1, "max");
        } else {
            sol = jline.api.mapqn.Mapqn_qr_bounds_rsrd.solve((Mapqn_qr_bounds_rsrd_parameters) params, 1, "max");
        }
        // An unsolved LP comes back as objective NaN with an EMPTY variable
        // map, which reads downstream as a utilization of exactly zero and
        // hence a throughput of zero -- a wrong answer that looks like a
        // computed one. Fail loudly instead.
        if (Double.isNaN(sol.getObjectiveValue())) {
            InputOutput.line_error("Solver_ctmc_qrf_analyzer",
                    "The QRF bound LP did not solve (OSQP returned no solution). Check "
                            + "options.qrfParams: BB (buffer capacity), F (per-station population "
                            + "ceiling) and ZZ must describe a feasible blocking network.");
        }
        double[] UBounds = new double[M];
        for (int i = 0; i < M; i++) {
            UBounds[i] = sol.getUtilization(i + 1);
        }

        // A closed network with N > 0 jobs cannot have every station idle: an
        // all-zero utilization vector means the LP settled on the trivial point
        // rather than on the bound, and must not be reported as a result.
        double uTotal = 0.0;
        for (int i = 0; i < M; i++) uTotal += UBounds[i];
        if (N > 0 && uTotal <= 0.0) {
            InputOutput.line_error("Solver_ctmc_qrf_analyzer",
                    "The QRF bound LP returned zero utilization at every station of a closed "
                            + "network, i.e. the trivial point rather than a bound. Check "
                            + "options.qrfParams: MM/MM1/ZZ must encode the population constraint "
                            + "of the blocking network.");
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
