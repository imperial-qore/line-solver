package jline.solvers.mva.handlers;

import java.util.Arrays;
import java.util.List;

import jline.api.mapqn.Mapqn_amva;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * SolverMVA method 'amva.mapqn': the horizontal-cut mean value analysis (Mapqn_amva) of a
 * closed multiclass model with one exponential delay station and one FCFS single-server
 * queue whose class-r service is a MAP. {@link #mapqnReason} is the one structural
 * predicate behind the method: listValidMethods drops the name when it is nonempty,
 * supportsModelMethod reports it, and {@link #solver_mapqn} raises it, so the offered,
 * reported and run answers cannot drift apart. Mirrors matlab/src/solvers/MVA/
 * mva_mapqn_reason.m and solver_mva_mapqn_analyzer.m.
 */
public final class Solver_mapqn {
    private Solver_mapqn() {}

    private static final List<ProcessType> MARKOVIAN = Arrays.asList(
            ProcessType.EXP, ProcessType.ERLANG, ProcessType.HYPEREXP, ProcessType.PH, ProcessType.APH,
            ProcessType.COXIAN, ProcessType.COX2, ProcessType.MAP, ProcessType.MMPP2);

    private static boolean isDelay(NetworkStruct sn, int i) {
        Station st = sn.stations.get(i);
        return sn.sched.get(st) == SchedStrategy.INF || Double.isInfinite(sn.nservers.get(i));
    }

    /** The reason 'amva.mapqn' cannot solve sn, or "" when it can. */
    public static String mapqnReason(NetworkStruct sn) {
        if (sn == null) return "Method 'amva.mapqn' needs a network model.";
        final int M = sn.nstations;
        final int R = sn.nclasses;
        if (M != 2)
            return "Method 'amva.mapqn' requires exactly two stations: one delay (infinite server) and one FCFS queue.";
        int id = -1, nDelay = 0;
        for (int i = 0; i < M; i++) if (isDelay(sn, i)) { id = i; nDelay++; }
        if (nDelay != 1)
            return "Method 'amva.mapqn' requires exactly one delay (infinite-server) station and one queue.";
        final int iq = 1 - id;
        final Station sq = sn.stations.get(iq);
        if (sn.sched.get(sq) != SchedStrategy.FCFS)
            return "Method 'amva.mapqn' requires FCFS scheduling at the queue; station " + (iq + 1) + " is "
                    + sn.sched.get(sq) + ".";
        if (sn.nservers.get(iq) != 1.0)
            return "Method 'amva.mapqn' supports a single-server queue only.";
        for (int r = 0; r < R; r++)
            if (Double.isInfinite(sn.njobs.get(r)))
                return "Method 'amva.mapqn' supports closed models only.";
        if (sn.nclosedjobs <= 0)
            return "Method 'amva.mapqn' supports closed models only.";
        final Station sd = sn.stations.get(id);
        final int fd = (int) sn.stationToStateful.get(id);
        final int fq = (int) sn.stationToStateful.get(iq);
        final double tol = 1e-12;
        for (int r = 0; r < R; r++) {
            if (sn.njobs.get(r) <= 0) continue;
            final JobClass jc = sn.jobclasses.get(r);
            final ProcessType pd = sn.procid.get(sd) == null ? null : sn.procid.get(sd).get(jc);
            if (pd != ProcessType.EXP)
                return "Method 'amva.mapqn' requires exponential think times; class " + (r + 1) + " has a " + pd
                        + " think time.";
            final ProcessType pq = sn.procid.get(sq) == null ? null : sn.procid.get(sq).get(jc);
            final MatrixCell cell = sn.proc.get(sq) == null ? null : sn.proc.get(sq).get(jc);
            if (pq == null || !MARKOVIAN.contains(pq) || cell == null || cell.size() < 2)
                return "Method 'amva.mapqn' requires a Markovian (MAP-representable) service process at the queue; class "
                        + (r + 1) + " is not.";
            if (Math.abs(sn.rt.get(fd * R + r, fq * R + r) - 1.0) > tol
                    || Math.abs(sn.rt.get(fq * R + r, fd * R + r) - 1.0) > tol)
                return "Method 'amva.mapqn' requires every class to cycle delay -> queue -> delay without class switching; class "
                        + (r + 1) + " does not.";
        }
        return "";
    }

    public static MVAResult solver_mapqn(NetworkStruct sn, SolverOptions options) {
        final long t0 = System.nanoTime();
        final String reason = mapqnReason(sn);
        if (!reason.isEmpty()) throw new RuntimeException(reason);
        final int M = sn.nstations;
        final int R = sn.nclasses;
        int id = -1;
        for (int i = 0; i < M; i++) if (isDelay(sn, i)) id = i;
        final int iq = 1 - id;
        final Station sq = sn.stations.get(iq);
        final int[] N = new int[R];
        final double[] mu = new double[R];
        final Matrix[] D0s = new Matrix[R];
        final Matrix[] D1s = new Matrix[R];
        for (int r = 0; r < R; r++) {
            N[r] = (int) Math.round(sn.njobs.get(r));
            if (N[r] > 0) {
                mu[r] = sn.rates.get(id, r);
                MatrixCell cell = sn.proc.get(sq).get(sn.jobclasses.get(r));
                D0s[r] = cell.get(0);
                D1s[r] = cell.get(1);
            } else {
                mu[r] = 1.0;                       // absent class: inert single phase
                D0s[r] = new Matrix(1, 1); D0s[r].set(0, 0, -1.0);
                D1s[r] = new Matrix(1, 1); D1s[r].set(0, 0, 1.0);
            }
        }
        final Mapqn_amva.Result res = Mapqn_amva.solve(mu, D0s, D1s, N);
        Matrix QN = new Matrix(M, R), UN = new Matrix(M, R), RN = new Matrix(M, R), TN = new Matrix(M, R);
        Matrix CN = new Matrix(1, R), XN = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            if (N[r] <= 0 || res.X[r] <= 0) continue;
            final double X = res.X[r];
            XN.set(0, r, X);
            QN.set(iq, r, res.Qq[r]); UN.set(iq, r, res.U[r]); TN.set(iq, r, X); RN.set(iq, r, res.Qq[r] / X);
            QN.set(id, r, X / mu[r]); UN.set(id, r, X / mu[r]); TN.set(id, r, X); RN.set(id, r, 1.0 / mu[r]);
            CN.set(0, r, N[r] / X);
        }
        MVAResult result = new MVAResult();
        result.QN = QN; result.UN = UN; result.RN = RN; result.TN = TN; result.CN = CN; result.XN = XN;
        int iter = 1;
        for (int r = 0; r < R; r++) iter *= (N[r] + 1);
        result.iter = iter;
        result.logNormConstAggr = Double.NaN;
        result.runtime = (System.nanoTime() - t0) / 1e9;
        result.method = "amva.mapqn";
        return result;
    }
}
