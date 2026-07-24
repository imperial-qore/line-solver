package jline.solvers.mam.handlers;

import jline.api.mam.Mam_is_renewal_map;
import jline.api.mam.Map_acf;
import jline.api.mam.Map_exponential;
import jline.api.mam.Qbd_setupdelayoff;
import jline.lang.NodeParam;
import jline.lang.nodeparam.QueueNodeParam;
import jline.lang.processes.Distribution;
import jline.api.mam.Map_lambda;
import jline.api.mam.Map_mean;
import jline.api.mam.Map_pie;
import jline.api.mam.Map_scale;
import jline.api.mam.Mmap_exponential;
import jline.api.mam.Mmap_lambda;
import jline.api.mam.Mmap_mark;
import jline.api.mam.Mmap_normalize;
import jline.api.mam.Mmap_scale;
import jline.api.mam.Mmap_shorten;
import jline.api.mam.Mmap_super_safe;
import jline.api.mam.Qbd_raprap1;
import jline.api.mam.QbdRapRap1Result;
import jline.api.mam.Qbd_setupdelayoff;
import jline.api.qsys.Qsys_dmc;
import jline.api.qsys.Qsys_mapdc;
import jline.api.qsys.Qsys_mdc_crommelin;
import jline.api.qsys.Qsys_phmc;
import jline.api.qsys.Qsys_phm1;
import jline.api.sn.SnGetDemandsChain;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Queue;
import jline.lang.processes.APH;
import jline.lang.nodes.Station;
import jline.GlobalConstants;
import jline.lang.constant.NodeType;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SchedStrategy;
import jline.lib.butools.MMAPPH1FCFS;
import jline.lib.qmam.Q_CT_MAP_MAP_1;
import jline.lib.qmam.MAPMAP1Options;
import jline.lib.qmam.MAPMAP1Result;
import jline.lib.butools.MMAPPH1NPPR;
import jline.lib.butools.MMAPPH1PRPR;
import jline.solvers.SolverOptions;
import jline.solvers.mam.MAMResult;
import jline.util.Utils;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import static jline.io.InputOutput.line_warning;
import static jline.io.InputOutput.line_warning_always;
import static jline.io.InputOutput.mfilename;

public final class Solver_mam_basic {

    private Solver_mam_basic() {
    }

    /**
     * True if the service process of any class at this station is a
     * matrix-exponential (ME) or rational-arrival-process (RAP)
     * representation. Neither is phase-type, so MMAPPH1FCFS -- which reads the
     * service only through the (pie, D0) phase-type pair -- returns the wrong
     * numbers for them and the RAP/RAP/1 QBD must be used instead.
     */
    /**
     * Render a server count the way MATLAB and Python %g do, so that a warning
     * assembled here is byte-identical to the one the other two codebases
     * assemble. Java's %g always keeps six significant digits and would print
     * a single server as "1.00000".
     */
    private static String fmtServers(double nservers) {
        if (Double.isInfinite(nservers)) return nservers > 0 ? "Inf" : "-Inf";
        if (Double.isNaN(nservers)) return "NaN";
        if (nservers == Math.rint(nservers) && Math.abs(nservers) < 1e15) {
            return String.valueOf((long) nservers);
        }
        return String.format("%g", nservers);
    }

    /**
     * True if station {@code jst}'s class-0 process is renewal. Returns false
     * when no usable (D0, D1) pair is present, which is the conservative side:
     * callers use this to decide whether a marginal-only closed form may be
     * applied, and applying one to a process whose correlation structure could
     * not be established is the failure this guard exists to prevent. Mirrors
     * MATLAB {@code mam_srcproc_is_renewal.m}.
     */
    private static boolean srcProcIsRenewal(NetworkStruct sn, int jst) {
        if (jst < 0 || jst >= sn.nstations) return false;
        Map<JobClass, MatrixCell> procMap = sn.proc.get(sn.stations.get(jst));
        if (procMap == null) return false;
        MatrixCell cell = procMap.get(sn.jobclasses.get(0));
        if (cell == null || cell.size() < 2) return false;
        Matrix D0 = cell.get(0);
        Matrix D1 = cell.get(1);
        if (D0 == null || D1 == null || D0.getNumRows() != D1.getNumRows()) return false;
        if (D0.hasNaN() || D1.hasNaN()) return false;
        return Mam_is_renewal_map.mam_is_renewal_map(D0, D1);
    }

    /**
     * True unless chain {@code c} is fed by an ME or RAP source, in which case
     * the assembled system-arrival (D0, D1) pair is legitimately non-Markovian
     * and {@code mmap_normalize} must not be applied to it: that routine clips
     * negative entries and re-derives the diagonal, which repairs roundoff on a
     * genuine MAP but substitutes a DIFFERENT, Markovian process with different
     * autocorrelation on an ME or RAP. Mirrors MATLAB
     * {@code mam_chain_arrival_is_markovian.m}.
     */
    private static boolean chainArrivalIsMarkovian(NetworkStruct sn, int c) {
        Matrix inchain = sn.inchain.get(c);
        if (inchain == null || inchain.length() == 0) return true;
        double totJobs = 0.0;
        for (int i = 0; i < inchain.length(); i++) {
            totJobs += sn.njobs.get((int) inchain.get(i));
        }
        if (!Utils.isInf(totJobs)) return true;   // closed chain: Poisson surrogate
        int jst = (int) sn.refstat.get((int) inchain.get(0));
        if (jst < 0 || jst >= sn.nstations) return true;
        Map<JobClass, ProcessType> classMap = sn.procid.get(sn.stations.get(jst));
        if (classMap == null) return true;
        for (int i = 0; i < inchain.length(); i++) {
            ProcessType pt = classMap.get(sn.jobclasses.get((int) inchain.get(i)));
            if (pt == ProcessType.ME || pt == ProcessType.RAP) return false;
        }
        return true;
    }

    private static boolean isMEorRAPService(NetworkStruct sn, int stationIdx) {
        Station station = sn.stations.get(stationIdx);
        Map<JobClass, ProcessType> classMap = sn.procid.get(station);
        if (classMap == null) return false;
        for (JobClass jobClass : sn.jobclasses) {
            ProcessType procType = classMap.get(jobClass);
            if (procType == ProcessType.RAP || procType == ProcessType.ME) return true;
        }
        return false;
    }

    /**
     * Analyse a single-class RAP/RAP/1 queue.
     *
     * The arrival process is any (D0,D1) pair: a Poisson stream, a MAP and a
     * PH renewal stream are all RAPs (a MAP is a RAP whose matrices happen to
     * be nonnegative), so the algorithm applies whenever the SERVICE is ME or
     * RAP, irrespective of the arrival type.
     *
     * qbd_raprap1 returns the full stationary queue-length distribution
     * pqueue, so every queue-length moment and the queue-length distribution
     * are recovered here and returned under the same keys MMAPPH1FCFS
     * populates ("ncMoms", "ncDistr"), class index 0. The algorithm carries no
     * sojourn-time information, so "stMoms"/"stDistr" are deliberately absent:
     * the caller derives the mean response time from Little's law
     * (RN = QN/TN, see the post-loop rescale below), and higher sojourn-time
     * moments are not obtainable from this QBD without a separate
     * passage-time analysis.
     */
    private static Map<String, Map<Integer, Matrix>> solveRapRap1(
            MatrixCell arrivalProc, MatrixCell serviceRAP, int stationIdx, int numQLMoms) {
        QbdRapRap1Result rapResult;
        try {
            rapResult = Qbd_raprap1.qbd_raprap1(arrivalProc, serviceRAP);
        } catch (RuntimeException e) {
            throw new RuntimeException("RAP/RAP/1 analysis failed at station " + stationIdx
                    + ", which has a matrix-exponential or rational service process: "
                    + e.getMessage(), e);
        }

        Matrix pqueue = rapResult.getPqueue();
        int numLevels = pqueue.getNumRows();
        int numPhases = pqueue.getNumCols();
        Matrix levelProb = new Matrix(1, numLevels, numLevels);
        for (int i = 0; i < numLevels; i++) {
            double p = 0.0;
            for (int j = 0; j < numPhases; j++) {
                p += pqueue.get(i, j);
            }
            levelProb.set(0, i, p);
        }

        Map<String, Map<Integer, Matrix>> result = new HashMap<String, Map<Integer, Matrix>>();

        int nMoms = FastMath.max(1, numQLMoms);
        Matrix moms = new Matrix(1, nMoms, nMoms);
        for (int m = 1; m <= nMoms; m++) {
            double mom = 0.0;
            for (int i = 0; i < numLevels; i++) {
                mom += FastMath.pow((double) i, m) * levelProb.get(0, i);
            }
            moms.set(0, m - 1, mom);
        }
        // The first moment is reported by qbd_raprap1 itself; keep its value so
        // the mean queue length agrees exactly with QbdRapRap1Result.getQN().
        moms.set(0, 0, rapResult.getQN());
        Map<Integer, Matrix> ncMoms = new HashMap<Integer, Matrix>();
        ncMoms.put(0, moms);
        result.put("ncMoms", ncMoms);

        Map<Integer, Matrix> ncDistr = new HashMap<Integer, Matrix>();
        ncDistr.put(0, levelProb);
        result.put("ncDistr", ncDistr);

        return result;
    }

    private static int getSpaceMax(Object config) {
        try {
            java.lang.reflect.Field f = config.getClass().getField("space_max");
            return ((Number) f.get(config)).intValue();
        } catch (Exception e) {
            return 0;
        }
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    public static MAMResult solver_mam_basic(NetworkStruct sn, SolverOptions options) {
        Object config = options.config;
        double tol = options.tol;

        // see _kb/06-solver-catalog.md for rationale
        Set<Integer> meWarned = new HashSet<Integer>();

        Map<Station, Map<JobClass, MatrixCell>> PH = sn.proc;
        int I = sn.nnodes;
        int M = sn.nstations;
        int K = sn.nclasses;
        int C = sn.nchains;
        Matrix N = sn.njobs.transpose();
        Matrix V = Matrix.cellsum(sn.visits);
        Matrix S = Matrix.ones(sn.rates.getNumRows(), sn.rates.getNumCols());
        S = S.elementDivide(sn.rates);
        Matrix Strue = S.copy(); // service times as declared, used to report utilization below
        Matrix Lchain = SnGetDemandsChain.snGetDemandsChain(sn).Dchain;

        // see _kb/06-solver-catalog.md for rationale
        double[] slcjobs = new double[M];
        if (sn.isslc != null) {
            for (int k = 0; k < K; k++) {
                if (sn.isslc.get(k) == 1.0) {
                    int ist_k = (int) sn.refstat.get(k, 0);
                    if (!Utils.isInf(sn.nservers.get(ist_k))) {
                        slcjobs[ist_k] += sn.njobs.get(k);
                    }
                }
            }
        }
        for (int ist = 0; ist < M; ist++) {
            if (slcjobs[ist] > 0) {
                for (int k = 0; k < K; k++) {
                    if (sn.isslc.get(k) != 1.0) {
                        S.set(ist, k, S.get(ist, k) * (1 + slcjobs[ist]));
                    }
                }
                // The chain demands drive the throughput fixed point and must be
                // consistent with the inflated service times.
                for (int c = 0; c < Lchain.getNumCols(); c++) {
                    Lchain.set(ist, c, Lchain.get(ist, c) * (1 + slcjobs[ist]));
                }
            }
        }

        Matrix QN = new Matrix(M, K, M * K);
        Matrix UN = new Matrix(M, K, M * K);
        Matrix RN = new Matrix(M, K, M * K);
        Matrix TN = new Matrix(M, K, M * K);
        Matrix WN = new Matrix(M, K, M * K);
        Matrix AN = new Matrix(M, K, M * K);
        Matrix CN = new Matrix(1, K, K);
        Matrix XN = new Matrix(1, K, K);

        Map<Integer, MatrixCell> pie = new HashMap<Integer, MatrixCell>();
        Map<Integer, MatrixCell> D0 = new HashMap<Integer, MatrixCell>();

        Matrix lambda = new Matrix(1, C, C);
        Map<Integer, MatrixCell> chainSysArrivals = new HashMap<Integer, MatrixCell>();
        Matrix TN_1 = new Matrix(M, K, M * K);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < K; j++) {
                TN_1.set(i, j, Double.POSITIVE_INFINITY);
            }
        }

        boolean[] mapdcStations = new boolean[M];

        int it = 0;

        for (int ist = 0; ist < M; ist++) {
            SchedStrategy schd = sn.sched.get(sn.stations.get(ist));
            if (schd == SchedStrategy.FCFS || schd == SchedStrategy.HOL || schd == SchedStrategy.FCFSPRIO
                    || schd == SchedStrategy.FCFSPRPRIO || schd == SchedStrategy.PS) {
                pie.put(ist, new MatrixCell());
                D0.put(ist, new MatrixCell());
                for (int k = 0; k < K; k++) {
                    Map<JobClass, MatrixCell> stMap = PH.get(sn.stations.get(ist));
                    MatrixCell cell = stMap.get(sn.jobclasses.get(k));
                    // see _kb/06-solver-catalog.md for rationale
                    Map<JobClass, ProcessType> procIdMap = sn.procid.get(sn.stations.get(ist));
                    ProcessType procTypeK = procIdMap == null ? null : procIdMap.get(sn.jobclasses.get(k));
                    if (procTypeK == ProcessType.ME || procTypeK == ProcessType.RAP) {
                        double ratio = Map_mean.map_mean(cell.get(0), cell.get(1))
                                / (S.get(ist, k) / sn.nservers.get(ist));
                        MatrixCell scaled = new MatrixCell(2);
                        scaled.set(0, cell.get(0).scale(ratio));
                        scaled.set(1, cell.get(1).scale(ratio));
                        stMap.put(sn.jobclasses.get(k), scaled);
                    } else {
                        stMap.put(sn.jobclasses.get(k),
                                Map_scale.map_scale(cell.get(0), cell.get(1), S.get(ist, k) / sn.nservers.get(ist)));
                    }
                    cell = stMap.get(sn.jobclasses.get(k));
                    pie.get(ist).set(k, Map_pie.map_pie(cell.get(0), cell.get(1)));
                    D0.get(ist).set(k, cell.get(0));
                    if (D0.get(ist).get(k).hasNaN()) {
                        stMap.put(sn.jobclasses.get(k), Map_exponential.map_exponential(GlobalConstants.Immediate));
                        pie.get(ist).set(k, Matrix.singleton(1.0));
                        D0.get(ist).set(k, Matrix.singleton(-GlobalConstants.Immediate));
                    }
                }
            }
        }

        boolean isOpen = false;
        boolean isClosed = false;

        outer1:
        for (int i = 0; i < sn.njobs.getNumRows(); i++) {
            for (int j = 0; j < sn.njobs.getNumCols(); j++) {
                if (!Double.isFinite(sn.njobs.get(i, j))) { isOpen = true; break outer1; }
            }
        }

        outer2:
        for (int i = 0; i < sn.njobs.getNumRows(); i++) {
            for (int j = 0; j < sn.njobs.getNumCols(); j++) {
                if (Double.isFinite(sn.njobs.get(i, j))) { isClosed = true; break outer2; }
            }
        }

        boolean isMixed = isOpen && isClosed;

        // see _kb/06-solver-catalog.md for rationale
        boolean[] isslcchain = new boolean[C];
        for (int c = 0; c < C; c++) {
            Matrix inchain_c = sn.inchain.get(c);
            boolean allslc = inchain_c.length() > 0;
            for (int i = 0; i < inchain_c.length(); i++) {
                if (sn.isslc == null || sn.isslc.get((int) inchain_c.get(i)) != 1.0) {
                    allslc = false;
                    break;
                }
            }
            isslcchain[c] = allslc;
        }

        MatrixCell lambdas_inchain = new MatrixCell(C);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            lambdas_inchain.set(c, new Matrix(inchain.getNumCols(), 1, inchain.getNumCols()));
            int ist = (int) sn.refstat.get((int) inchain.value(), 0);
            if (!PH.containsKey(sn.stations.get(ist))) {
                PH.put(sn.stations.get(ist), new HashMap<JobClass, MatrixCell>());
            }
            for (int j = 0; j < inchain.getNumCols(); j++) {
                lambdas_inchain.get(c).set(j, 0, sn.rates.get(ist, (int) inchain.get(0, j)));
            }
            double sum = 0.0;
            for (int k = 0; k < lambdas_inchain.get(c).getNumRows(); k++) {
                if (Double.isFinite(lambdas_inchain.get(c).get(k, 0))) {
                    sum += lambdas_inchain.get(c).get(k, 0);
                }
            }
            lambda.set(0, c, sum);
            if (isslcchain[c]) {
                lambda.set(0, c, 0.0);
                for (int j = 0; j < inchain.getNumCols(); j++) {
                    lambdas_inchain.get(c).set(j, 0, 0.0);
                }
            }
            boolean openChain = false;
            for (int j = 0; j < inchain.getNumCols(); j++) {
                if (Utils.isInf(sn.njobs.get((int) inchain.get(0, j)))) { openChain = true; break; }
            }
            if (openChain) {
                for (int k = 0; k < K; k++) {
                    if (Double.isNaN(PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k)).get(0).get(0))) {
                        PH.get(sn.stations.get(ist)).put(sn.jobclasses.get(k),
                                Map_exponential.map_exponential(Double.POSITIVE_INFINITY));
                    }
                }
                double k = inchain.get(0);
                chainSysArrivals.put(c, new MatrixCell());
                chainSysArrivals.get(c).set(0, PH.get(sn.stations.get(ist)).get(sn.jobclasses.get((int) k)).get(0));
                chainSysArrivals.get(c).set(1, PH.get(sn.stations.get(ist)).get(sn.jobclasses.get((int) k)).get(1));
                chainSysArrivals.get(c).set(2, PH.get(sn.stations.get(ist)).get(sn.jobclasses.get((int) k)).get(1));
                for (int ki = 1; ki < inchain.length(); ki++) {
                    k = inchain.get(ki);
                    if (Double.isNaN(PH.get(sn.stations.get(ist)).get(sn.jobclasses.get((int) k)).get(0).get(0))) {
                        PH.get(sn.stations.get(ist)).put(sn.jobclasses.get((int) k),
                                Map_exponential.map_exponential(Double.POSITIVE_INFINITY));
                    }
                    Map<Integer, MatrixCell> MMAPS = new HashMap<Integer, MatrixCell>();
                    MMAPS.put(0, chainSysArrivals.get(c));
                    MMAPS.put(1, new MatrixCell());
                    MMAPS.get(1).set(0, PH.get(sn.stations.get(ist)).get(sn.jobclasses.get((int) k)).get(0));
                    MMAPS.get(1).set(1, PH.get(sn.stations.get(ist)).get(sn.jobclasses.get((int) k)).get(1));
                    MMAPS.get(1).set(2, PH.get(sn.stations.get(ist)).get(sn.jobclasses.get((int) k)).get(1));
                    chainSysArrivals.put(c, Mmap_super_safe.mmap_super_safe(MMAPS, getSpaceMax(config), "default"));
                }
                for (int i = 0; i < inchain.getNumCols(); i++) {
                    TN.set(ist, (int) inchain.get(0, i), lambdas_inchain.get(c).get(i, 0));
                }
            }
        }

        Matrix sd = new Matrix(sn.nservers.getNumRows(), sn.nservers.getNumCols(), sn.nservers.getNonZeroLength());
        for (int i = 0; i < sn.nservers.getNumRows(); i++) {
            if (!Utils.isInf(sn.nservers.get(i, 0))) sd.set(i, 0, 1);
        }

        boolean[] isclosedchain = new boolean[C];
        boolean anyClosedChain = false;
        boolean anyOpenChain = false;
        for (int c = 0; c < C; c++) {
            Matrix inchain_c = sn.inchain.get(c);
            boolean openChain = false;
            for (int j = 0; j < inchain_c.getNumCols(); j++) {
                if (Utils.isInf(sn.njobs.get((int) inchain_c.get(0, j)))) { openChain = true; break; }
            }
            isclosedchain[c] = !openChain && !isslcchain[c];
            anyClosedChain = anyClosedChain || isclosedchain[c];
            anyOpenChain = anyOpenChain || openChain;
        }
        // see _kb/06-solver-catalog.md for rationale
        boolean ismixed = anyClosedChain && anyOpenChain;
        // see _kb/06-solver-catalog.md for rationale
        double Ulim = 1 - GlobalConstants.CoarseTol;

        Matrix dif_matrix = TN.add(-1.0, TN_1);
        double dif = FastMath.abs(dif_matrix.elementMaxAbs());

        while (dif > tol && it <= options.iter_max) {
            it++;
            TN_1 = TN.copy();
            // see _kb/06-solver-catalog.md for rationale
            double Umax = -Double.POSITIVE_INFINITY;
            boolean sawFiniteServerRow = false;
            boolean allRowsNaN = true;
            for (int i = 0; i < M; i++) {
                if (sd.get(i, 0) == 0.0) {
                    continue;
                }
                sawFiniteServerRow = true;
                double sum = 0.0;
                for (int j = 0; j < K; j++) sum += UN.get(i, j);
                if (Double.isNaN(sum)) {
                    continue;
                }
                allRowsNaN = false;
                if (sum > Umax) {
                    Umax = sum;
                }
            }
            if (sawFiniteServerRow && allRowsNaN) {
                Umax = Double.NaN;
            }
            if (ismixed || Umax < 1) {
                for (int c = 0; c < C; c++) {
                    Matrix inchain = sn.inchain.get(c);
                    if (isclosedchain[c]) {
                        double Nc = 0.0;
                        for (int i = 0; i < inchain.length(); i++) {
                            Nc = Nc + sn.njobs.get((int) inchain.get(i));
                        }
                        Matrix QN_omitnan = new Matrix(QN);
                        QN_omitnan.removeNaN();
                        Matrix QN_col_sum = new Matrix(inchain.length(), 1, inchain.length());
                        for (int i = 0; i < inchain.length(); i++) {
                            QN_col_sum.set(i, 0, QN_omitnan.sumCols((int) inchain.get(i)));
                        }
                        double QNc = QN_col_sum.elementSum();
                        QNc = FastMath.max(options.tol, QNc);
                        double TNlb = Nc / Lchain.sumCols(c);
                        if (it == 1) {
                            lambda.set(0, c, TNlb);
                        } else {
                            lambda.set(0, c,
                                    lambda.get(0, c) * it / options.iter_max
                                            + (Nc / QNc) * lambda.get(0, c) * (options.iter_max - it) / options.iter_max);
                        }
                    }
                }
            }
            if (ismixed) {
                // see _kb/06-solver-catalog.md for rationale
                double theta = Double.POSITIVE_INFINITY;
                boolean anyBinding = false;
                for (int i = 0; i < M; i++) {
                    if (sd.get(i, 0) == 0.0) continue;
                    double uOpen = 0.0;
                    double uClosed = 0.0;
                    for (int c = 0; c < C; c++) {
                        double u = Lchain.get(i, c) * lambda.get(0, c);
                        if (isclosedchain[c]) {
                            uClosed += u;
                        } else {
                            uOpen += u;
                        }
                    }
                    if (uClosed > tol) {
                        anyBinding = true;
                        theta = FastMath.min(theta, (Ulim - uOpen) / uClosed);
                    }
                }
                if (anyBinding && theta < 1) {
                    double scale = FastMath.max(0.0, theta);
                    for (int c = 0; c < C; c++) {
                        if (isclosedchain[c]) {
                            lambda.set(0, c, lambda.get(0, c) * scale);
                        }
                    }
                }
            } else if (Umax >= 1) {
                lambda.divideEq(Umax);
            }

            for (int c = 0; c < C; c++) {
                Matrix inchain = sn.inchain.get(c);
                double njobsChain = 0.0;
                for (int i = 0; i < inchain.length(); i++) {
                    njobsChain += sn.njobs.get((int) inchain.get(i));
                }
                if (!Double.isInfinite(njobsChain)) {
                    // see _kb/06-solver-catalog.md for rationale
                    Matrix lambda_c = new Matrix(1, 1, 1);
                    lambda_c.set(0, 0, lambda.get(c));
                    chainSysArrivals.put(c, Mmap_exponential.mmap_exponential(lambda_c));
                }
                for (int m = 0; m < M; m++) {
                    for (int i = 0; i < inchain.length(); i++) {
                        TN.set(m, (int) inchain.get(i), V.get(m, (int) inchain.get(i)) * lambda.get(c));
                    }
                }
            }

            for (int ind = 0; ind < I; ind++) {
                if (sn.isstation.get(ind) == 1.0) {
                    int ist = (int) sn.nodeToStation.get(ind);
                    if (sn.nodetype.get(ind) == NodeType.Join) {
                        for (int c = 0; c < C; c++) {
                            Matrix inchain = sn.inchain.get(c);
                            for (int i = 0; i < inchain.length(); i++) {
                                int fanin = 0;
                                for (int j = 0; j < sn.rtnodes.getNumRows(); j++) {
                                    if (sn.rtnodes.get(j, (int) ((ind - 1) * K + inchain.get(i))) != 0.0) fanin++;
                                }
                                int k = (int) inchain.get(i);
                                TN.set(ist, k, lambda.get(c) * V.get(ist, k) / fanin);
                                UN.set(ist, k, 0.0);
                                QN.set(ist, k, 0.0);
                                RN.set(ist, k, 0.0);
                            }
                        }
                    } else if (sn.nodetype.get(ind) == NodeType.Fork) {
                        for (int c = 0; c < C; c++) {
                            Matrix inchain = sn.inchain.get(c);
                            for (int i = 0; i < inchain.length(); i++) {
                                int k = (int) inchain.get(i);
                                TN.set(ist, k, lambda.get(c) * V.get(ist, k));
                                UN.set(ist, k, 0.0);
                                QN.set(ist, k, 0.0);
                                RN.set(ist, k, 0.0);
                            }
                        }
                    } else {
                        SchedStrategy schd2 = sn.sched.get(sn.stations.get(ist));
                        if (schd2 == SchedStrategy.INF) {
                            for (int c = 0; c < C; c++) {
                                Matrix inchain = sn.inchain.get(c);
                                for (int i = 0; i < inchain.length(); i++) {
                                    int k = (int) inchain.get(i);
                                    if (V.get(ist, k) <= GlobalConstants.Zero) {
                                        // see _kb/06-solver-catalog.md for rationale
                                        TN.set(ist, k, 0.0);
                                        UN.set(ist, k, 0.0);
                                        QN.set(ist, k, 0.0);
                                        RN.set(ist, k, 0.0);
                                        continue;
                                    }
                                    TN.set(ist, k, lambda.get(c) * V.get(ist, k));
                                    // see _kb/06-solver-catalog.md for rationale
                                    UN.set(ist, k, S.get(ist, k) * TN.get(ist, k));
                                    // see _kb/06-solver-catalog.md for rationale
                                    QN.set(ist, k, TN.get(ist, k) * S.get(ist, k));
                                    RN.set(ist, k, QN.get(ist, k) / TN.get(ist, k));
                                }
                            }
                        } else if (schd2 == SchedStrategy.PS) {
                            for (int c = 0; c < C; c++) {
                                Matrix inchain = sn.inchain.get(c);
                                for (int i = 0; i < inchain.length(); i++) {
                                    int k = (int) inchain.get(i);
                                    if (V.get(ist, k) <= GlobalConstants.Zero) {
                                        // see _kb/06-solver-catalog.md for rationale
                                        TN.set(ist, k, 0.0);
                                        UN.set(ist, k, 0.0);
                                        continue;
                                    }
                                    TN.set(ist, k, lambda.get(c) * V.get(ist, k));
                                    // Utilization Law: a c-server station holds TN*S/c of its capacity.
                                    UN.set(ist, k, S.get(ist, k) * TN.get(ist, k) / sn.nservers.get(ist));
                                }
                                double Uden = FastMath.min(1 - GlobalConstants.FineTol, UN.sumRows(ist));
                                for (int i = 0; i < inchain.length(); i++) {
                                    int k = (int) inchain.get(i);
                                    if (V.get(ist, k) <= GlobalConstants.Zero) {
                                        QN.set(ist, k, 0.0);
                                        RN.set(ist, k, 0.0);
                                        continue;
                                    }
                                    QN.set(ist, k, UN.get(ist, k) / (1 - Uden));
                                    RN.set(ist, k, QN.get(ist, k) / TN.get(ist, k));
                                }
                            }
                        } else if (schd2 == SchedStrategy.HOL || schd2 == SchedStrategy.FCFSPRIO
                                || schd2 == SchedStrategy.FCFS || schd2 == SchedStrategy.FCFSPRPRIO) {
                            Map<Integer, MatrixCell> chainArrivalAtNode = new HashMap<Integer, MatrixCell>();
                            Map<Integer, MatrixCell> rates = new HashMap<Integer, MatrixCell>();
                            MatrixCell aggrArrivalAtNode = new MatrixCell();
                            for (int c = 0; c < C; c++) {
                                Matrix a = Matrix.extractRows(V, ist, ist + 1, null);
                                a.scaleEq(lambda.get(c));
                                if (c == 0) rates.put(ist, new MatrixCell());
                                rates.get(ist).set(c, a);
                                Matrix inchain = sn.inchain.get(c);
                                Matrix markProb = new Matrix(1, inchain.length(), inchain.length());
                                for (int i = 0; i < inchain.length(); i++) {
                                    markProb.set(0, i, rates.get(ist).get(c).get((int) inchain.get(i)));
                                }
                                markProb.scaleEq(1.0 / markProb.elementSum());
                                markProb.removeNaN();
                                boolean chainIsMarkovian = chainArrivalIsMarkovian(sn, c);
                                chainArrivalAtNode.put(c, Mmap_mark.mmap_mark(chainSysArrivals.get(c), markProb));
                                // see _kb/06-solver-catalog.md for rationale
                                if (chainIsMarkovian) {
                                    chainArrivalAtNode.put(c, Mmap_normalize.mmap_normalize(chainArrivalAtNode.get(c)));
                                }
                                // see _kb/06-solver-catalog.md for rationale
                                boolean anyPosRate = false;
                                for (int i = 0; i < inchain.length(); i++) {
                                    if (rates.get(ist).get(c).get((int) inchain.get(i)) > 0) {
                                        anyPosRate = true;
                                        break;
                                    }
                                }
                                if (anyPosRate) {
                                    Matrix b = new Matrix(1, inchain.length(), inchain.length());
                                    for (int i = 0; i < inchain.length(); i++) {
                                        b.set(0, i, 1 / rates.get(ist).get(c).get((int) inchain.get(i)));
                                    }
                                    chainArrivalAtNode.put(c, Mmap_scale.mmap_scale(chainArrivalAtNode.get(c), b));
                                }
                                if (c == 0) {
                                    if (chainIsMarkovian) {
                                        Map<Integer, MatrixCell> MMAPS = new HashMap<Integer, MatrixCell>();
                                        MMAPS.put(0, chainArrivalAtNode.get(c));
                                        MMAPS.put(1, Mmap_exponential.mmap_exponential(Matrix.singleton(0.0), 1));
                                        MatrixCell aggrArrivalAtNode_temp =
                                                Mmap_super_safe.mmap_super_safe(MMAPS, getSpaceMax(config), "default");
                                        aggrArrivalAtNode.set(0, aggrArrivalAtNode_temp.get(0));
                                        aggrArrivalAtNode.set(1, aggrArrivalAtNode_temp.get(1));
                                        aggrArrivalAtNode.set(2, aggrArrivalAtNode_temp.get(1));
                                    } else {
                                        // see _kb/06-solver-catalog.md for rationale
                                        aggrArrivalAtNode.set(0, chainArrivalAtNode.get(c).get(0));
                                        aggrArrivalAtNode.set(1, chainArrivalAtNode.get(c).get(1));
                                        aggrArrivalAtNode.set(2, chainArrivalAtNode.get(c).get(1));
                                    }
                                    double lc = Map_lambda.map_lambda(chainArrivalAtNode.get(c).get(0),
                                            chainArrivalAtNode.get(c).get(1));
                                    if (lc > 0) {
                                        aggrArrivalAtNode = Mmap_scale.mmap_scale(aggrArrivalAtNode,
                                                Matrix.singleton(1.0 / lc));
                                    }
                                } else {
                                    if (!chainIsMarkovian) {
                                        // see _kb/06-solver-catalog.md for rationale
                                        line_warning_always(mfilename(new Object() {}),
                                                "Chain %d has a matrix-exponential or rational arrival process, "
                                                        + "which the superposition of several chains at station %s "
                                                        + "cannot represent exactly. Its autocorrelation is "
                                                        + "approximated by the closest Markovian arrival process.",
                                                c, sn.stations.get(ist).getName());
                                    }
                                    Map<Integer, MatrixCell> MMAPS = new HashMap<Integer, MatrixCell>();
                                    MMAPS.put(0, aggrArrivalAtNode);
                                    MMAPS.put(1, chainArrivalAtNode.get(c));
                                    aggrArrivalAtNode =
                                            Mmap_super_safe.mmap_super_safe(MMAPS, getSpaceMax(config), "default");
                                }
                            }
                            Map<Integer, Matrix> Qret = new HashMap<Integer, Matrix>();
                            PriorityAnalysis priorityAnalysis = Solver_mam_passage_time.analyzePriorities(sn, "basic");

                            SchedStrategy schedStrat = sn.sched.get(sn.stations.get(ist));
                            if ((schedStrat == SchedStrategy.HOL || schedStrat == SchedStrategy.FCFSPRIO
                                    || schedStrat == SchedStrategy.FCFSPRPRIO) && !priorityAnalysis.isIdentical) {
                                if (priorityAnalysis.isAllDistinct) {
                                    Map<Double, Integer> priorityMap = new LinkedHashMap<Double, Integer>();
                                    for (int k = 0; k < K; k++) {
                                        priorityMap.put(sn.classprio.get(k), k);
                                    }
                                    List<Double> sortedPriorities = new ArrayList<Double>(priorityMap.keySet());
                                    Collections.sort(sortedPriorities, Collections.reverseOrder());
                                    int[] iK = new int[K];
                                    for (int i = 0; i < sortedPriorities.size(); i++) {
                                        iK[i] = priorityMap.get(sortedPriorities.get(i));
                                    }

                                    MatrixCell reorderedArrival = new MatrixCell();
                                    reorderedArrival.set(0, aggrArrivalAtNode.get(0));
                                    for (int i = 0; i < iK.length; i++) {
                                        // MMAP cell layout: 0=D0, 1=D1 total, 2+k=class-k marked
                                        // matrix; BUTools expects {D0, Dlowest..Dhighest}
                                        reorderedArrival.set(1 + i, aggrArrivalAtNode.get(2 + iK[i]));
                                    }

                                    Map<Integer, Matrix> reorderedPie = new LinkedHashMap<Integer, Matrix>();
                                    Map<Integer, Matrix> reorderedD0 = new LinkedHashMap<Integer, Matrix>();
                                    for (int i = 0; i < iK.length; i++) {
                                        reorderedPie.put(i, pie.get(ist).get(iK[i]));
                                        reorderedD0.put(i, D0.get(ist).get(iK[i]));
                                    }

                                    MatrixCell reorderedPieCell = new MatrixCell();
                                    MatrixCell reorderedD0Cell = new MatrixCell();
                                    for (Map.Entry<Integer, Matrix> e : reorderedPie.entrySet()) {
                                        reorderedPieCell.set(e.getKey(), e.getValue());
                                    }
                                    for (Map.Entry<Integer, Matrix> e : reorderedD0.entrySet()) {
                                        reorderedD0Cell.set(e.getKey(), e.getValue());
                                    }

                                    Map<String, Map<Integer, Matrix>> prioResult;
                                    if (schedStrat == SchedStrategy.FCFSPRPRIO) {
                                        prioResult = MMAPPH1PRPR.MMAPPH1PRPR(
                                                reorderedArrival, reorderedPieCell, reorderedD0Cell,
                                                1, null, null, null, null, null, null);
                                    } else {
                                        prioResult = MMAPPH1NPPR.MMAPPH1NPPR(
                                                reorderedArrival, reorderedPieCell, reorderedD0Cell,
                                                1, null, null, null, null, null, null);
                                    }

                                    Map<Integer, Matrix> prioResultMoms = new HashMap<Integer, Matrix>();
                                    if (prioResult.get("ncMoms") != null) {
                                        prioResultMoms.putAll(prioResult.get("ncMoms"));
                                    }

                                    Qret = new HashMap<Integer, Matrix>();
                                    for (int i = 0; i < iK.length; i++) {
                                        if (prioResultMoms.containsKey(i)) {
                                            Qret.put(iK[i], prioResultMoms.get(i));
                                        }
                                    }
                                } else {
                                    throw new RuntimeException(priorityAnalysis.message);
                                }
                            } else {
                                Matrix sn_rates_ist_k = new Matrix(1, K, K);
                                for (int i = 0; i < K; i++) {
                                    sn_rates_ist_k.set(i, sn.rates.get(ist, i) * sn.nservers.get(ist));
                                }
                                Matrix lambdaVec = Mmap_lambda.mmap_lambda(aggrArrivalAtNode);
                                Matrix ratesVec = sn_rates_ist_k.elementIncrease(GlobalConstants.FineTol);
                                double aggrUtil = 0.0;
                                if (lambdaVec.getNumCols() == ratesVec.getNumCols()) {
                                    for (int i = 0; i < lambdaVec.getNumCols(); i++) {
                                        double lam = lambdaVec.get(i);
                                        double rate = ratesVec.get(i);
                                        if (Double.isFinite(lam) && Double.isFinite(rate) && rate > 0) {
                                            aggrUtil += lam / rate;
                                        }
                                    }
                                } else if (lambdaVec.getNumCols() == 1) {
                                    double totalLambda = lambdaVec.get(0);
                                    if (Double.isFinite(totalLambda)) {
                                        for (int i = 0; i < K; i++) {
                                            double rate = ratesVec.get(i);
                                            if (Double.isFinite(rate) && rate > 0) {
                                                aggrUtil += totalLambda / rate;
                                            }
                                        }
                                    }
                                } else {
                                    int minCols = Math.min(lambdaVec.getNumCols(), ratesVec.getNumCols());
                                    for (int i = 0; i < minCols; i++) {
                                        double lam = lambdaVec.get(i);
                                        double rate = ratesVec.get(i);
                                        if (Double.isFinite(lam) && Double.isFinite(rate) && rate > 0) {
                                            aggrUtil += lam / rate;
                                        }
                                    }
                                }
                                if (aggrUtil < 1 - GlobalConstants.FineTol) {
                                    boolean closed = true;
                                    for (int i = 0; i < N.length(); i++) {
                                        if (Utils.isInf(N.get(i))) { closed = false; break; }
                                    }

                                    // see _kb/06-solver-catalog.md for rationale
                                    boolean isFunctionOpen = !closed && sn.isfunction != null
                                            && ist < sn.isfunction.getNumRows()
                                            && sn.isfunction.get(ist, 0) == 1.0;

                                    boolean isMapDc = false;
                                    boolean isDMc = false;
                                    int dmcSourceIdx = -1;
                                    boolean isPhM1 = false;
                                    int phM1SourceIdx = -1;
                                    if (!closed && K == 1 && !isFunctionOpen) {
                                        Station station = sn.stations.get(ist);
                                        JobClass jobClass = sn.jobclasses.get(0);
                                        Map<JobClass, ProcessType> classMap = sn.procid.get(station);
                                        ProcessType procType = classMap == null ? null : classMap.get(jobClass);
                                        isMapDc = (procType == ProcessType.DET);
                                        if (!isMapDc && procType == ProcessType.EXP) {
                                            for (int jst = 0; jst < sn.nstations; jst++) {
                                                if (jst == ist) continue;
                                                Station srcStation = sn.stations.get(jst);
                                                Map<JobClass, ProcessType> srcMap = sn.procid.get(srcStation);
                                                ProcessType srcProc = srcMap == null ? null : srcMap.get(jobClass);
                                                if (srcProc == ProcessType.DET) {
                                                    isDMc = true; dmcSourceIdx = jst; break;
                                                }
                                            }
                                        }
                                        if (!isMapDc && !isDMc && procType == ProcessType.EXP
                                                && !Utils.isInf(sn.nservers.get(ist)) && sn.nservers.get(ist) >= 1.0) {
                                            for (int jst = 0; jst < sn.nstations; jst++) {
                                                if (jst == ist) continue;
                                                Station srcStation = sn.stations.get(jst);
                                                Map<JobClass, ProcessType> srcMap = sn.procid.get(srcStation);
                                                ProcessType srcProc = srcMap == null ? null : srcMap.get(jobClass);
                                                if (srcProc == null) continue;
                                                if (srcProc != ProcessType.EXP && srcProc != ProcessType.DET
                                                        && srcProc != ProcessType.IMMEDIATE
                                                        && srcProc != ProcessType.DISABLED) {
                                                    // see _kb/06-solver-catalog.md for rationale
                                                    if (!srcProcIsRenewal(sn, jst)) continue;
                                                    isPhM1 = true; phM1SourceIdx = jst; break;
                                                } else if (srcProc == ProcessType.EXP && sn.nservers.get(ist) > 1.0
                                                        && sn.nstations == 2
                                                        && sn.nodetype.get((int) sn.stationToNode.get(jst)) == NodeType.Source) {
                                                    // see _kb/06-solver-catalog.md for rationale
                                                    isPhM1 = true; phM1SourceIdx = jst; break;
                                                }
                                            }
                                        }
                                    }

                                    if (isPhM1) {
                                        double muQ = S.get(ist, 0) > 0.0 ? 1.0 / S.get(ist, 0) : Double.POSITIVE_INFINITY;
                                        jline.util.Pair<double[], double[][]> phPair =
                                                Qsys_phm1.extractPhPairForPhm1(sn, phM1SourceIdx, 0);
                                        if (phPair != null) {
                                            try {
                                                double[] alphaVec = phPair.getFirst();
                                                double[][] Tmat = phPair.getSecond();
                                                int cServ = (int) sn.nservers.get(ist);
                                                jline.api.qsys.PhMcResult res = Qsys_phmc.qsys_phmc(
                                                        alphaVec, Tmat, muQ, cServ);
                                                Qret.put(0, Matrix.singleton(res.getMeanQueueLength()));
                                                mapdcStations[ist] = true;
                                            } catch (Exception e) {
                                                isPhM1 = false;
                                            }
                                        } else {
                                            isPhM1 = false;
                                        }
                                    }
                                    if (isFunctionOpen) {
                                        // see _kb/06-solver-catalog.md for rationale
                                        Distribution setupDist = null;
                                        Distribution delayOffDist = null;
                                        NodeParam np = sn.nodeparam.get(sn.stations.get(ist));
                                        if (np instanceof QueueNodeParam) {
                                            QueueNodeParam qnp = (QueueNodeParam) np;
                                            for (int k = 0; k < K; k++) {
                                                JobClass jobClass = sn.jobclasses.get(k);
                                                if (qnp.setupTime.get(jobClass) != null
                                                        && qnp.delayoffTime.get(jobClass) != null) {
                                                    setupDist = qnp.setupTime.get(jobClass);
                                                    delayOffDist = qnp.delayoffTime.get(jobClass);
                                                }
                                            }
                                        }
                                        if (setupDist == null || delayOffDist == null) {
                                            throw new RuntimeException(
                                                    "Station " + ist + " is marked as a setup/delay-off station "
                                                            + "but carries no setup or delay-off distribution.");
                                        }
                                        double alpharate = 1.0 / setupDist.getMean();
                                        double alphascv = setupDist.getSCV();
                                        double betarate = 1.0 / delayOffDist.getMean();
                                        double betascv = delayOffDist.getSCV();

                                        // rho_k = lambda_k / mu_k is the per-class term of
                                        // aggrUtil, so aggrUtil is the aggregate load.
                                        double rhoTotal = aggrUtil;
                                        double aggrLambdaTotal = 0.0;
                                        for (int i = 0; i < lambdaVec.getNumCols(); i++) {
                                            double lam = lambdaVec.get(i);
                                            if (Double.isFinite(lam)) aggrLambdaTotal += lam;
                                        }
                                        Qret = new HashMap<Integer, Matrix>();
                                        if (rhoTotal > 0) {
                                            // see _kb/06-solver-catalog.md for rationale
                                            double aggrRate = aggrLambdaTotal / rhoTotal;
                                            double qTotal = Qbd_setupdelayoff.qbd_setupdelayoff(
                                                    aggrLambdaTotal, aggrRate, alpharate, alphascv, betarate, betascv);
                                            for (int k = 0; k < K; k++) {
                                                double lam = k < lambdaVec.getNumCols() ? lambdaVec.get(k) : 0.0;
                                                double rate = k < ratesVec.getNumCols() ? ratesVec.get(k) : 0.0;
                                                double rhoK = (Double.isFinite(lam) && Double.isFinite(rate) && rate > 0)
                                                        ? lam / rate : 0.0;
                                                Qret.put(k, Matrix.singleton(qTotal * rhoK / rhoTotal));
                                            }
                                        } else {
                                            for (int k = 0; k < K; k++) {
                                                Qret.put(k, Matrix.singleton(0.0));
                                            }
                                        }
                                    } else if (isPhM1) {
                                        // already handled
                                    } else if (isDMc) {
                                        int numServers = (int) sn.nservers.get(ist);
                                        double muQ = S.get(ist, 0) > 0.0 ? 1.0 / S.get(ist, 0) : Double.POSITIVE_INFINITY;
                                        double lamD = sn.rates.get(dmcSourceIdx, 0);
                                        try {
                                            jline.api.qsys.DmcResult res = Qsys_dmc.qsys_dmc(lamD, muQ, numServers);
                                            Qret.put(0, Matrix.singleton(res.getMeanQueueLength()));
                                            mapdcStations[ist] = true;
                                        } catch (Exception e) {
                                            Map<Integer, Matrix> pieMap = new HashMap<Integer, Matrix>();
                                            Map<Integer, Matrix> d0Map = new HashMap<Integer, Matrix>();
                                            MatrixCell pieCell = pie.get(ist);
                                            MatrixCell d0Cell = D0.get(ist);
                                            for (int kk = 0; kk < pieCell.size(); kk++) pieMap.put(kk, pieCell.get(kk));
                                            for (int kk = 0; kk < d0Cell.size(); kk++) d0Map.put(kk, d0Cell.get(kk));
                                            Map<String, Map<Integer, Matrix>> r = MMAPPH1FCFS.MMAPPH1FCFS(
                                                    Mmap_shorten.mmap_shorten(aggrArrivalAtNode),
                                                    pieMap, d0Map, 1, null, null, null, false, false, null, null);
                                            Qret = new HashMap<Integer, Matrix>();
                                            if (r.get("ncMoms") != null) Qret.putAll(r.get("ncMoms"));
                                        }
                                    } else if (isMapDc) {
                                        MatrixCell arrivalMAP = Mmap_shorten.mmap_shorten(aggrArrivalAtNode);
                                        Matrix D0_arr = arrivalMAP.get(0);
                                        Matrix D1_arr = arrivalMAP.get(1);
                                        double detServiceTime = S.get(ist, 0);
                                        int numServers = (int) sn.nservers.get(ist);

                                        try {
                                            boolean isPoisson = D0_arr.getNumRows() == 1 && D1_arr.getNumRows() == 1
                                                    && D1_arr.get(0, 0) > 0.0
                                                    && Math.abs(D0_arr.get(0, 0) + D1_arr.get(0, 0)) < 1e-12;
                                            if (isPoisson) {
                                                double lambdaPoisson = D1_arr.get(0, 0);
                                                jline.api.qsys.MDcCrommelinResult crom =
                                                        Qsys_mdc_crommelin.qsys_mdc_crommelin(lambdaPoisson, detServiceTime, numServers);
                                                Qret.put(0, Matrix.singleton(crom.getMeanQueueLength()));
                                                mapdcStations[ist] = true;
                                            } else {
                                                jline.api.qsys.QsysMapDcResult mapdcResult =
                                                        Qsys_mapdc.qsys_mapdc(D0_arr, D1_arr, detServiceTime, numServers);
                                                Qret.put(0, Matrix.singleton(mapdcResult.getMeanQueueLength()));
                                            }
                                        } catch (Exception e) {
                                            Map<Integer, Matrix> pieMap = new HashMap<Integer, Matrix>();
                                            Map<Integer, Matrix> d0Map = new HashMap<Integer, Matrix>();
                                            MatrixCell pieCell = pie.get(ist);
                                            MatrixCell d0Cell = D0.get(ist);
                                            for (int kk = 0; kk < pieCell.size(); kk++) pieMap.put(kk, pieCell.get(kk));
                                            for (int kk = 0; kk < d0Cell.size(); kk++) d0Map.put(kk, d0Cell.get(kk));
                                            Map<String, Map<Integer, Matrix>> r = MMAPPH1FCFS.MMAPPH1FCFS(
                                                    Mmap_shorten.mmap_shorten(aggrArrivalAtNode),
                                                    pieMap, d0Map, 1, null, null, null, false, false, null, null);
                                            Qret = new HashMap<Integer, Matrix>();
                                            if (r.get("ncMoms") != null) Qret.putAll(r.get("ncMoms"));
                                        }
                                    } else if (!closed) {
                                        // see _kb/06-solver-catalog.md for rationale
                                        MatrixCell svcMap = PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(0));
                                        // see _kb/06-solver-catalog.md for rationale
                                        boolean useRapRap1 = false;
                                        if (isMEorRAPService(sn, ist) && svcMap != null && svcMap.size() >= 2) {
                                            if (K == 1 && sn.nservers.get(ist) == 1.0) {
                                                useRapRap1 = true;
                                            } else {
                                                // see _kb/06-solver-catalog.md for rationale
                                                if (meWarned.add(ist)) {
                                                    line_warning_always(mfilename(new Object() {}),
                                                            "Station %s has a matrix-exponential or rational service process, "
                                                                    + "which the RAP/RAP/1 analysis supports only with a single class "
                                                                    + "at a single server (here %d classes, %s servers). Falling back "
                                                                    + "to the phase-type approximation MMAPPH1FCFS, which is not exact "
                                                                    + "for this service process.",
                                                            sn.stations.get(ist).getName(), K, fmtServers(sn.nservers.get(ist)));
                                                }
                                            }
                                        }
                                        boolean useMapMap1 = false;
                                        if (!useRapRap1 && K == 1 && sn.nservers.get(ist) == 1.0
                                                && svcMap != null && svcMap.size() >= 2) {
                                            double acf1 = Map_acf.map_acf(svcMap.get(0), svcMap.get(1), 1).get(0);
                                            useMapMap1 = Math.abs(acf1) > GlobalConstants.CoarseTol;
                                        }
                                        if (useRapRap1) {
                                            MatrixCell arvMap = Mmap_shorten.mmap_shorten(aggrArrivalAtNode);
                                            MatrixCell arvRAP = new MatrixCell(2);
                                            arvRAP.set(0, arvMap.get(0));
                                            arvRAP.set(1, arvMap.get(1));
                                            MatrixCell svcRAP = new MatrixCell(2);
                                            svcRAP.set(0, svcMap.get(0));
                                            svcRAP.set(1, svcMap.get(1));
                                            Map<String, Map<Integer, Matrix>> r =
                                                    solveRapRap1(arvRAP, svcRAP, ist, 1);
                                            Qret = new HashMap<Integer, Matrix>();
                                            Qret.putAll(r.get("ncMoms"));
                                        } else if (useMapMap1) {
                                            MatrixCell arvMap = Mmap_shorten.mmap_shorten(aggrArrivalAtNode);
                                            MAPMAP1Result mm = Q_CT_MAP_MAP_1.qCtMapMap1(
                                                    arvMap.get(0), arvMap.get(1), svcMap.get(0), svcMap.get(1),
                                                    new MAPMAP1Options("SylvesCR", 100000, 0));
                                            Matrix ql = mm.getQueueLength();
                                            double en = 0.0;
                                            for (int n = 0; n < ql.getNumElements(); n++) en += n * ql.get(n);
                                            Qret = new HashMap<Integer, Matrix>();
                                            Qret.put(0, Matrix.singleton(en));
                                        } else {
                                            Map<Integer, Matrix> pieMap = new HashMap<Integer, Matrix>();
                                            Map<Integer, Matrix> d0Map = new HashMap<Integer, Matrix>();
                                            MatrixCell pieCell = pie.get(ist);
                                            MatrixCell d0Cell = D0.get(ist);
                                            for (int kk = 0; kk < pieCell.size(); kk++) pieMap.put(kk, pieCell.get(kk));
                                            for (int kk = 0; kk < d0Cell.size(); kk++) d0Map.put(kk, d0Cell.get(kk));
                                            Map<String, Map<Integer, Matrix>> r = MMAPPH1FCFS.MMAPPH1FCFS(
                                                    Mmap_shorten.mmap_shorten(aggrArrivalAtNode),
                                                    pieMap, d0Map, 1, null, null, null, false, false, null, null);
                                            Qret = new HashMap<Integer, Matrix>();
                                            if (r.get("ncMoms") != null) Qret.putAll(r.get("ncMoms"));
                                        }
                                    } else {
                                        Matrix finite_N = N.copy();
                                        finite_N.removeInfinite();
                                        double maxLevel = finite_N.elementMax() + 1;
                                        MatrixCell D = Mmap_shorten.mmap_shorten(aggrArrivalAtNode);
                                        Map<Integer, Matrix> pdistr = new HashMap<Integer, Matrix>();
                                        if (Map_lambda.map_lambda(D.get(0), D.get(1)) < GlobalConstants.FineTol) {
                                            for (int k = 0; k < K; k++) {
                                                Matrix pdistrK = new Matrix(1, 2, 2);
                                                pdistrK.set(0, 1 - GlobalConstants.FineTol);
                                                pdistrK.set(1, GlobalConstants.FineTol);
                                                pdistr.put(k, pdistrK);
                                                Qret.put(k, Matrix.singleton(GlobalConstants.FineTol / sn.rates.get(ist)));
                                            }
                                        } else {
                                            Station station = sn.stations.get(ist);
                                            if (station instanceof Queue && ((Queue) station).isDelayOffEnabled()) {
                                                double alpharate = 1.0;
                                                double betarate = 1.0;
                                                double betascv = 1.0;
                                                for (int k = 0; k < K; k++) {
                                                    JobClass jobClass = sn.jobclasses.get(k);
                                                    Object setupDist = ((Queue) station).getSetupTime(jobClass);
                                                    Object delayOffDist = ((Queue) station).getDelayOffTime(jobClass);
                                                    if (setupDist != null && delayOffDist != null) {
                                                        try {
                                                            java.lang.reflect.Method getMean = setupDist.getClass().getMethod("getMean");
                                                            java.lang.reflect.Method getSCV = setupDist.getClass().getMethod("getSCV");
                                                            alpharate = 1.0 / ((Number) getMean.invoke(setupDist)).doubleValue();
                                                            betarate = 1.0 / ((Number) getMean.invoke(delayOffDist)).doubleValue();
                                                            betascv = ((Number) getSCV.invoke(delayOffDist)).doubleValue();
                                                        } catch (Exception ex) {
                                                            // ignore
                                                        }
                                                        break;
                                                    }
                                                }

                                                // see _kb/06-solver-catalog.md for rationale
                                                double setupMean = 1.0 / alpharate;
                                                for (int c = 0; c < C; c++) {
                                                    Matrix inchain = sn.inchain.get(c);
                                                    // see _kb/06-solver-catalog.md for rationale
                                                    double Vtot = 0.0;
                                                    for (int i = 0; i < inchain.length(); i++) {
                                                        Vtot += V.get(ist, (int) inchain.get(i));
                                                    }
                                                    double ZT = 0.0;
                                                    for (int ist2 = 0; ist2 < M; ist2++) {
                                                        if (Double.isInfinite(sn.nservers.get(ist2))) {
                                                            double d = Lchain.get(ist2, c);
                                                            if (!Double.isNaN(d)) {
                                                                ZT += d;
                                                            }
                                                        }
                                                    }
                                                    ZT = ZT / FastMath.max(Vtot, GlobalConstants.FineTol);
                                                    double nu = 1.0 / FastMath.max(ZT, GlobalConstants.FineTol);
                                                    double pcold;
                                                    if (betascv == 1.0) {
                                                        pcold = betarate / (betarate + nu);
                                                    } else {
                                                        // general delay-off: LST of the fitted APH
                                                        // at nu, pcold = alpha*(nu I - T)^-1*(-T 1)
                                                        APH betaAPH = APH.fitMeanAndSCV(1.0 / betarate, betascv);
                                                        Matrix Tb = (Matrix) betaAPH.getParam(3).getValue();
                                                        Matrix alphaB = (Matrix) betaAPH.getParam(2).getValue();
                                                        int nbb = Tb.getNumRows();
                                                        Matrix A = Matrix.zeros(nbb, nbb);
                                                        Matrix tvec = new Matrix(nbb, 1, nbb);
                                                        for (int i = 0; i < nbb; i++) {
                                                            double rowSum = 0.0;
                                                            for (int j = 0; j < nbb; j++) {
                                                                A.set(i, j, (i == j ? nu : 0.0) - Tb.get(i, j));
                                                                rowSum -= Tb.get(i, j);
                                                            }
                                                            tvec.set(i, 0, rowSum);
                                                        }
                                                        Matrix x = A.inv().mult(tvec);
                                                        pcold = 0.0;
                                                        for (int i = 0; i < nbb; i++) {
                                                            pcold += alphaB.get(i) * x.get(i);
                                                        }
                                                    }
                                                    for (int i = 0; i < inchain.length(); i++) {
                                                        int k = (int) inchain.get(i);
                                                        double lamK = rates.get(ist).get(c).get(k);
                                                        if (Double.isFinite(S.get(ist, k)) && lamK > 0) {
                                                            // see _kb/06-solver-catalog.md for rationale
                                                            Qret.put(k, Matrix.singleton(
                                                                    lamK * (pcold * setupMean
                                                                            + S.get(ist, k) / sn.nservers.get(ist))));
                                                        } else {
                                                            // see _kb/06-solver-catalog.md for rationale
                                                            Qret.put(k, Matrix.singleton(0.0));
                                                        }
                                                    }
                                                }
                                            } else {
                                                Map<Integer, Matrix> pieMap = new HashMap<Integer, Matrix>();
                                                Map<Integer, Matrix> d0Map = new HashMap<Integer, Matrix>();
                                                MatrixCell pieCell = pie.get(ist);
                                                MatrixCell d0Cell = D0.get(ist);
                                                for (int kk = 0; kk < pieCell.size(); kk++) pieMap.put(kk, pieCell.get(kk));
                                                for (int kk = 0; kk < d0Cell.size(); kk++) d0Map.put(kk, d0Cell.get(kk));
                                                Map<String, Map<Integer, Matrix>> r = MMAPPH1FCFS.MMAPPH1FCFS(
                                                        D, pieMap, d0Map, null, (int) maxLevel, null, null, false, false, null, null);
                                                pdistr = new HashMap<Integer, Matrix>();
                                                if (r.get("ncDistr") != null) pdistr.putAll(r.get("ncDistr"));
                                                for (int k = 0; k < K; k++) {
                                                    Matrix pdistrK = pdistr.get(k);
                                                    if (pdistrK == null) {
                                                        Qret.put(k, Matrix.singleton(GlobalConstants.FineTol / sn.rates.get(ist)));
                                                        continue;
                                                    }
                                                    pdistr.put(k, Matrix.extractRows(pdistrK.transpose(), 0, (int) N.get(k) + 1, null));
                                                    pdistr.get(k).absEq();
                                                    double sumP = 0.0;
                                                    for (int i = 0; i < pdistr.get(k).length() - 1; i++) {
                                                        sumP = sumP + pdistr.get(k).get(i);
                                                    }
                                                    pdistr.get(k).set(pdistr.get(k).length() - 1, Math.abs(1 - sumP));
                                                    pdistr.get(k).scaleEq(1 / pdistr.get(k).elementSum());
                                                    Matrix a = new Matrix(1, (int) N.get(k) + 1, (int) N.get(k) + 1);
                                                    for (int i = 0; i < a.length(); i++) {
                                                        a.set(i, (double) i);
                                                    }
                                                    Matrix b = new Matrix(1, (int) N.get(k) + 1, (int) N.get(k) + 1);
                                                    for (int i = 0; i < a.length(); i++) {
                                                        b.set(i, pdistr.get(k).get(i));
                                                    }
                                                    Qret.put(k, Matrix.singleton(Math.max(0.0, Math.min(N.get(k), a.mult(b.transpose()).get(0)))));
                                                }
                                            }
                                        }
                                    }
                                } else {
                                    for (int k = 0; k < K; k++) {
                                        double njobValue = sn.njobs.get(k);
                                        Qret.put(k, Matrix.singleton(njobValue));
                                    }
                                }
                            }
                            for (int i = 0; i < Qret.size(); i++) {
                                QN.set(ist, i, Qret.get(i).get(0));
                            }
                            boolean isFunctionStation = sn.stations.get(ist) instanceof Queue
                                    && ((Queue) sn.stations.get(ist)).isDelayOffEnabled();
                            for (int k = 0; k < K; k++) {
                                int c = 0;
                                for (int i = 0; i < sn.chains.getNumRows(); i++) {
                                    if (sn.chains.get(i, k) != 0.0) { c = i; break; }
                                }
                                TN.set(ist, k, rates.get(ist).get(c).get(k));
                                UN.set(ist, k, TN.get(ist, k) * S.get(ist, k) / sn.nservers.get(ist));
                                if (isFunctionStation && !Double.isFinite(UN.get(ist, k))) {
                                    // see _kb/06-solver-catalog.md for rationale
                                    UN.set(ist, k, 0.0);
                                }
                                QN.set(ist, k, Qret.get(k).get(0));
                                if (!mapdcStations[ist]) {
                                    // see _kb/06-solver-catalog.md for rationale
                                    if (Double.isFinite(S.get(ist, k))) {
                                        QN.set(ist, k,
                                                QN.get(ist, k) + TN.get(ist, k) * S.get(ist, k)
                                                        * (sn.nservers.get(ist) - 1) / sn.nservers.get(ist));
                                    }
                                }
                                if (V.get(ist, k) <= GlobalConstants.Zero) {
                                    RN.set(ist, k, 0.0);
                                } else {
                                    RN.set(ist, k, QN.get(ist, k) / TN.get(ist, k));
                                }
                            }
                        }
                    }
                } else {
                    // other node types handled by default traffic flow
                }
            }
            dif_matrix = TN.add(-1.0, TN_1);
            dif = dif_matrix.elementMaxAbs();
        }
        int totiter = it + 2;

        CN = RN.sumCols();
        QN.absEq();
        for (int pass = 1; pass <= 2; pass++) {
            for (int c = 0; c < C; c++) {
                Matrix inchain = sn.inchain.get(c);
                double Nc = 0.0;
                for (int j = 0; j < inchain.length(); j++) {
                    Nc = Nc + sn.njobs.get((int) inchain.get(j));
                }
                if (Double.isFinite(Nc)) {
                    double QNc = 0.0;
                    boolean hasNaN = false;
                    for (int j = 0; j < inchain.length(); j++) {
                        int colIdx = (int) inchain.get(j);
                        for (int rowIdx = 0; rowIdx < QN.getNumRows(); rowIdx++) {
                            double val_ij = QN.get(rowIdx, colIdx);
                            if (Double.isNaN(val_ij)) {
                                hasNaN = true;
                            } else if (Double.isFinite(val_ij)) {
                                QNc += val_ij;
                            }
                        }
                    }
                    if (hasNaN) QNc = Double.NaN;
                    if (Double.isFinite(QNc) && QNc > GlobalConstants.FineTol) {
                        for (int m = 0; m < inchain.length(); m++) {
                            for (int n = 0; n < QN.getNumRows(); n++) {
                                double oldVal = QN.get(n, (int) inchain.get(m));
                                if (Double.isFinite(oldVal)) {
                                    QN.set(n, (int) inchain.get(m), oldVal * (Nc / QNc));
                                }
                            }
                        }
                    } else if (Double.isNaN(QNc)) {
                        for (int m = 0; m < inchain.length(); m++) {
                            for (int n = 0; n < QN.getNumRows(); n++) {
                                double oldVal = QN.get(n, (int) inchain.get(m));
                                if (Double.isFinite(oldVal)) {
                                    QN.set(n, (int) inchain.get(m), Double.NaN);
                                }
                            }
                        }
                    }
                }

                for (int ind = 0; ind < I; ind++) {
                    for (int k = 0; k < inchain.length(); k++) {
                        if (sn.isstation.get(ind) == 1.0) {
                            int ist = (int) sn.nodeToStation.get(ind);
                            if (mapdcStations[ist]) continue;
                            if (V.get(ist, (int) inchain.get(k)) > GlobalConstants.Zero) {
                                if (Utils.isInf(sn.nservers.get(ist))) {
                                    RN.set(ist, (int) inchain.get(k), S.get(ist, (int) inchain.get(k)));
                                } else {
                                    double s_val = S.get(ist, (int) inchain.get(k));
                                    double tn_val = TN.get(ist, (int) inchain.get(k));
                                    double qn_tn;
                                    if (tn_val > GlobalConstants.Zero) {
                                        qn_tn = QN.get(ist, (int) inchain.get(k)) / tn_val;
                                    } else {
                                        qn_tn = Double.NaN;
                                    }
                                    double rn_val;
                                    if (Double.isNaN(qn_tn)) rn_val = s_val;
                                    else if (Double.isNaN(s_val)) rn_val = qn_tn;
                                    else rn_val = Math.max(s_val, qn_tn);
                                    RN.set(ist, (int) inchain.get(k), rn_val);
                                }
                            } else {
                                RN.set(ist, (int) inchain.get(k), 0);
                            }
                            QN.set(ist, (int) inchain.get(k),
                                    RN.get(ist, (int) inchain.get(k)) * TN.get(ist, (int) inchain.get(k)));
                        }
                    }
                }
                if (Nc == 0.0) {
                    for (int k = 0; k < QN.getNumRows(); k++) QN.set(k, c, 0);
                    for (int k = 0; k < UN.getNumRows(); k++) UN.set(k, c, 0);
                    for (int k = 0; k < RN.getNumRows(); k++) RN.set(k, c, 0);
                    for (int k = 0; k < TN.getNumRows(); k++) TN.set(k, c, 0);
                    CN.set(0, c, 0);
                    XN.set(0, c, 0);
                }
            }
        }

        // see _kb/06-solver-catalog.md for rationale
        boolean anyslc = false;
        if (sn.isslc != null) {
            for (int k = 0; k < K; k++) {
                if (sn.isslc.get(k) == 1.0) { anyslc = true; break; }
            }
        }
        if (anyslc) {
            // see _kb/06-solver-catalog.md for rationale
            for (int ist = 0; ist < M; ist++) {
                if (slcjobs[ist] > 0) {
                    for (int k = 0; k < K; k++) {
                        if (sn.isslc.get(k) != 1.0) {
                            UN.set(ist, k, Strue.get(ist, k) * TN.get(ist, k));
                        }
                    }
                }
            }
            for (int k = 0; k < K; k++) {
                if (sn.isslc.get(k) != 1.0) continue;
                for (int i = 0; i < M; i++) {
                    QN.set(i, k, 0.0);
                    UN.set(i, k, 0.0);
                    RN.set(i, k, 0.0);
                    TN.set(i, k, 0.0);
                }
            }
            for (int ist = 0; ist < M; ist++) {
                List<Integer> slck = new ArrayList<Integer>();
                for (int k = 0; k < K; k++) {
                    if (sn.isslc.get(k) == 1.0 && (int) sn.refstat.get(k, 0) == ist) {
                        slck.add(k);
                    }
                }
                if (slck.isEmpty()) continue;
                if (Utils.isInf(sn.nservers.get(ist))) {
                    // Delay station: no contention, every customer is always in
                    // service, so the class completes at its full aggregate rate.
                    for (int idx = 0; idx < slck.size(); idx++) {
                        int k = slck.get(idx);
                        QN.set(ist, k, sn.njobs.get(k));
                        TN.set(ist, k, sn.njobs.get(k) * sn.rates.get(ist, k));
                        RN.set(ist, k, Strue.get(ist, k));
                        UN.set(ist, k, Strue.get(ist, k) * TN.get(ist, k));
                    }
                } else {
                    double nsrv = sn.nservers.get(ist);
                    double Uother = 0.0;
                    for (int k = 0; k < K; k++) {
                        if (sn.isslc.get(k) != 1.0) Uother += UN.get(ist, k);
                    }
                    double Uleft = FastMath.max(0.0, 1 - Uother);
                    // Several self-looping classes at one station share the free
                    // capacity in proportion to the service rate they offer.
                    double wsum = 0.0;
                    double[] w = new double[slck.size()];
                    for (int idx = 0; idx < slck.size(); idx++) {
                        int k = slck.get(idx);
                        w[idx] = sn.njobs.get(k) * sn.rates.get(ist, k);
                        wsum += w[idx];
                    }
                    if (wsum <= 0) continue;
                    for (int idx = 0; idx < slck.size(); idx++) {
                        int k = slck.get(idx);
                        double ucap = FastMath.min(sn.njobs.get(k), nsrv) / nsrv;
                        double uk = FastMath.min(Uleft * w[idx] / wsum, ucap);
                        UN.set(ist, k, uk);
                        TN.set(ist, k, sn.rates.get(ist, k) * uk * nsrv);
                        QN.set(ist, k, sn.njobs.get(k));
                        if (TN.get(ist, k) > 0) {
                            RN.set(ist, k, QN.get(ist, k) / TN.get(ist, k));
                        } else {
                            RN.set(ist, k, 0.0);
                        }
                    }
                }
            }
            for (int k = 0; k < K; k++) {
                CN.set(0, k, RN.sumCols(k));
            }
        }

        MAMResult result = new MAMResult();
        result.QN = QN;
        result.UN = UN;
        result.RN = RN;
        result.TN = TN;
        result.CN = CN;
        result.XN = XN;
        result.WN = WN;
        result.AN = AN;
        result.iter = totiter;
        return result;
    }
}
