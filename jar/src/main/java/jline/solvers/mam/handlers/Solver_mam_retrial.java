/**
 * @file MAM Retrial Queue Solver
 *
 * @since LINE 3.0
 */
package jline.solvers.mam.handlers;

import jline.VerboseLevel;
import jline.api.qsys.Qsys_bmapphnn_retrial;
import jline.api.qsys.Qsys_is_retrial;
import jline.api.qsys.RetrialInfo;
import jline.io.InputOutput;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.RetrialPolicy;
import jline.lang.constant.ProcessType;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.solvers.mam.MAMResult;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Solver_mam_retrial {
    private Solver_mam_retrial() {}

    /**
     * Solves BMAP/PH/N/N bufferless retrial queues using the MAM retrial solver.
     */
    public static MAMResult solver_mam_retrial(NetworkStruct sn, SolverOptions options) {
        RetrialInfo retInfo = Qsys_is_retrial.qsys_is_retrial(sn);
        if (!retInfo.isRetrial()) {
            InputOutput.line_error("solver_mam_retrial",
                    "No valid retrial configuration detected: " + retInfo.getErrorMsg());
            throw new RuntimeException("No valid retrial configuration detected: " + retInfo.getErrorMsg());
        }

        int M = sn.nstations;
        int K = sn.nclasses;

        Matrix QN = new Matrix(M, K, M * K);
        Matrix UN = new Matrix(M, K, M * K);
        Matrix RN = new Matrix(M, K, M * K);
        Matrix TN = new Matrix(M, K, M * K);
        Matrix CN = new Matrix(1, K, K);
        Matrix XN = new Matrix(1, K, K);

        int stationIdx = retInfo.getStationIdx();
        int sourceIdx = retInfo.getSourceIdx();
        int classIdx = retInfo.getClassIdx();
        int N = retInfo.getN();

        // see _kb/06-solver-catalog.md for rationale
        assertMarkovian(sn, sourceIdx, classIdx, "arrival");
        assertMarkovian(sn, stationIdx, classIdx, "service");
        warnIfApproximated(sn, stationIdx, classIdx);

        Object sourceStation = sn.stations.get(sourceIdx);
        Object jobClass = sn.jobclasses.get(classIdx);
        java.util.Map<?, MatrixCell> srcProcMap = sn.proc.get(sourceStation);
        MatrixCell arrivalProc = srcProcMap == null ? null : srcProcMap.get(jobClass);
        if (arrivalProc == null) {
            throw new RuntimeException("No arrival process at source station");
        }

        int numMats = arrivalProc.size();
        if (numMats < 2) {
            throw new RuntimeException("Invalid arrival process: need at least D0 and D1 matrices");
        }
        Matrix[] D = new Matrix[numMats];
        for (int k = 0; k < numMats; k++) {
            D[k] = arrivalProc.get(k);
            if (D[k] == null) {
                throw new RuntimeException("Null matrix at index " + k + " in arrival process");
            }
        }

        Object queueStation = sn.stations.get(stationIdx);
        java.util.Map<?, MatrixCell> qProcMap = sn.proc.get(queueStation);
        MatrixCell serviceProc = qProcMap == null ? null : qProcMap.get(jobClass);
        if (serviceProc == null) {
            throw new RuntimeException("No service process at queue station");
        }

        Matrix S = serviceProc.get(0);
        Matrix D1_service = serviceProc.get(1);
        if (S == null) {
            throw new RuntimeException("No D0 matrix in service process");
        }
        if (D1_service == null) {
            throw new RuntimeException("No D1 matrix in service process");
        }

        int nPhases = S.getNumRows();

        Matrix S0 = new Matrix(nPhases, 1);
        for (int i = 0; i < nPhases; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < nPhases; j++) {
                rowSum += S.get(i, j);
            }
            S0.set(i, 0, -rowSum);
        }

        Matrix beta = new Matrix(1, nPhases);
        int foundIdx = -1;
        for (int i = 0; i < nPhases; i++) {
            if (S0.get(i, 0) > 1e-10) {
                foundIdx = i;
                break;
            }
        }

        if (foundIdx >= 0) {
            for (int j = 0; j < nPhases; j++) {
                beta.set(0, j, D1_service.get(foundIdx, j) / S0.get(foundIdx, 0));
            }
        } else {
            for (int j = 0; j < nPhases; j++) {
                beta.set(0, j, 1.0 / nPhases);
            }
        }

        double betaSum = 0.0;
        for (int j = 0; j < nPhases; j++) {
            betaSum += beta.get(0, j);
        }
        if (Math.abs(betaSum - 1.0) > 1e-6 && betaSum > 0) {
            for (int j = 0; j < nPhases; j++) {
                beta.set(0, j, beta.get(0, j) / betaSum);
            }
        }

        double alpha = 0.1;
        if (sn.retrialProc != null) {
            java.util.Map<?, MatrixCell> retProcMap = sn.retrialProc.get(queueStation);
            MatrixCell retrialDistCell = retProcMap == null ? null : retProcMap.get(jobClass);
            if (retrialDistCell != null && retrialDistCell.size() >= 2) {
                Matrix retrialD0 = retrialDistCell.get(0);
                if (retrialD0 != null && retrialD0.getNumRows() >= 1 && retrialD0.getNumCols() >= 1) {
                    alpha = -retrialD0.get(0, 0);
                }
            }
        }

        double gamma = retInfo.getGamma();
        if (sn.orbitImpatience != null) {
            java.util.Map<?, MatrixCell> orbitMap = sn.orbitImpatience.get(queueStation);
            MatrixCell orbitDistCell = orbitMap == null ? null : orbitMap.get(jobClass);
            if (orbitDistCell != null && orbitDistCell.size() >= 2) {
                Matrix orbitD0 = orbitDistCell.get(0);
                if (orbitD0 != null && orbitD0.getNumRows() >= 1 && orbitD0.getNumCols() >= 1) {
                    // For Exp(gamma), D0 = -gamma
                    gamma = -orbitD0.get(0, 0);
                }
            }
        }
        double p = retInfo.getP();
        int R = retInfo.getR();

        // see _kb/06-solver-catalog.md for rationale
        int maxLevel = -1;
        double tailTol = Qsys_bmapphnn_retrial.DEFAULT_TAIL_TOLERANCE;
        if (options.config != null) {
            if (options.config.orbit_maxlevel > 0) {
                maxLevel = options.config.orbit_maxlevel;
            }
            if (options.config.orbit_tailtol > 0) {
                tailTol = options.config.orbit_tailtol;
            }
        }
        // see _kb/06-solver-catalog.md for rationale
        if (sn.orbitMaxJobs != null && stationIdx < sn.stations.size()) {
            java.util.Map<jline.lang.JobClass, Integer> omap =
                    sn.orbitMaxJobs.get(sn.stations.get(stationIdx));
            Integer cap = (omap != null) ? omap.get(sn.jobclasses.get(classIdx)) : null;
            if (cap != null && cap.intValue() >= 0) {
                maxLevel = cap.intValue();
            }
        }

        int retrialPolicy = RetrialPolicy.LINEAR;
        if (sn.retrialPolicy != null && stationIdx < sn.stations.size()) {
            java.util.Map<jline.lang.JobClass, Integer> pmap =
                    sn.retrialPolicy.get(sn.stations.get(stationIdx));
            Integer pol = (pmap != null) ? pmap.get(sn.jobclasses.get(classIdx)) : null;
            if (pol != null && pol.intValue() > 0) {
                retrialPolicy = pol.intValue();
            }
        }
        double tol = (options.tol > 0) ? options.tol : 1e-10;
        boolean verbose = options.verbose != VerboseLevel.SILENT;

        InputOutput.line_debug(options.verbose, "Calling qsys_bmapphnn_retrial: N=" + N
                + ", alpha=" + alpha + ", gamma=" + gamma + ", p=" + p + ", R=" + R);

        jline.api.qsys.QsysRetrialResult perf = Qsys_bmapphnn_retrial.qsys_bmapphnn_retrial(
                D, beta, S, N, alpha, gamma, p, R, maxLevel, tol, verbose,
                tailTol, Qsys_bmapphnn_retrial.DEFAULT_MAX_DIM, Qsys_bmapphnn_retrial.DEFAULT_MAX_BLOCK_SIZE,
                retrialPolicy);

        QN.set(stationIdx, classIdx, perf.L_orbit + perf.N_server);
        UN.set(stationIdx, classIdx, perf.Utilization);
        TN.set(stationIdx, classIdx, perf.Throughput);

        if (perf.Throughput > 0) {
            RN.set(stationIdx, classIdx, QN.get(stationIdx, classIdx) / perf.Throughput);
        } else {
            RN.set(stationIdx, classIdx, Double.POSITIVE_INFINITY);
        }

        XN.set(0, classIdx, perf.Throughput);
        CN.set(0, classIdx, RN.get(stationIdx, classIdx));

        MAMResult result = new MAMResult();
        result.retrialInternals = perf;
        result.QN = QN;
        result.UN = UN;
        result.RN = RN;
        result.TN = TN;
        result.CN = CN;
        result.XN = XN;
        result.iter = perf.truncLevel;
        result.method = "retrial";

        return result;
    }

    /**
     * Rejects a station-class process that is not a valid Markovian (D0,D1,...)
     * representation. Non-Markovian distributions (Det, traces, NHPP rate
     * schedules) and disabled classes reach here as empty cells or as NaN-filled
     * matrices; they have no place in a matrix-analytic generator.
     */
    private static void assertMarkovian(NetworkStruct sn, int ist, int r, String role) {
        String stationName = sn.nodenames.get((int) sn.stationToNode.get(ist));
        String className = sn.classnames.get(r);
        Station station = sn.stations.get(ist);
        JobClass jobClass = sn.jobclasses.get(r);
        java.util.Map<JobClass, MatrixCell> procMap = sn.proc == null ? null : sn.proc.get(station);
        MatrixCell proc = procMap == null ? null : procMap.get(jobClass);

        if (proc == null || proc.size() < 2) {
            InputOutput.line_error("solver_mam_retrial", String.format(
                    "The %s process of class '%s' at station '%s' is not a Markovian "
                    + "(D0,D1) representation. The matrix-analytic retrial solver requires phase-type or "
                    + "Markovian-arrival distributions; use SolverCTMC, SolverSSA or SolverLDES instead.",
                    role, className, stationName));
            return;
        }

        int n0 = proc.get(0) == null ? -1 : proc.get(0).getNumRows();
        for (int e = 0; e < proc.size(); e++) {
            Matrix De = proc.get(e);
            if (De == null || De.getNumRows() != De.getNumCols() || De.getNumRows() != n0) {
                InputOutput.line_error("solver_mam_retrial", String.format(
                        "The %s process of class '%s' at station '%s' has "
                        + "inconsistent (D0,D1) block dimensions.", role, className, stationName));
            }
            if (De.hasNaN() || De.hasInfinite()) {
                InputOutput.line_error("solver_mam_retrial", String.format(
                        "The %s process of class '%s' at station '%s' contains "
                        + "NaN or Inf entries: the class is disabled at this station or the distribution has no "
                        + "Markovian representation.", role, className, stationName));
            }
        }
    }

    /**
     * Non-Markovian service distributions are replaced by an Erlang MAP of
     * matching mean. That is a documented approximation, not the requested
     * distribution, so say so.
     */
    private static void warnIfApproximated(NetworkStruct sn, int ist, int r) {
        if (sn.procid == null) {
            return;
        }
        Station station = sn.stations.get(ist);
        JobClass jobClass = sn.jobclasses.get(r);
        java.util.Map<JobClass, ProcessType> procidMap = sn.procid.get(station);
        ProcessType ptype = procidMap == null ? null : procidMap.get(jobClass);
        if (ptype == null) {
            return;
        }
        switch (ptype) {
            case DET:
            case REPLAYER:
            case UNIFORM:
            case GAMMA:
            case PARETO:
            case WEIBULL:
            case LOGNORMAL:
                String stationName = sn.nodenames.get((int) sn.stationToNode.get(ist));
                InputOutput.line_warning("solver_mam_retrial", String.format(
                        "Service distribution %s at station '%s' is not phase-type. "
                        + "The matrix-analytic retrial solver uses an Erlang approximation matching its mean and "
                        + "squared coefficient of variation (%d phases); results are approximate.",
                        ProcessType.toText(ptype), stationName, (int) sn.phases.get(ist, r)));
                break;
            default:
                break;
        }
    }
}
