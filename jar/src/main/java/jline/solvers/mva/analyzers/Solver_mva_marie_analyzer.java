package jline.solvers.mva.analyzers;

import jline.api.pfqn.Pfqn_marie;
import jline.api.sn.SnDeaggregateChainResults;
import jline.api.sn.SnGetDemandsChain;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.util.matrix.Matrix;

/**
 * Marie's iterative aggregation-decomposition (Marie 1979/1980) for closed
 * networks with FCFS non-exponential (Coxian) service, wired as SolverMVA
 * method 'marie'. Infinite-server stations fold into per-chain think time;
 * Pfqn_marie is applied to the queueing stations, where only FCFS is service-
 * sensitive (its SCV is used) while insensitive product-form disciplines
 * (PS, LCFSPR) are forced to exponential.
 *
 * <p>Single chain: aggregate is the exact load-dependent product-form solve, so
 * exponential service is exact. Multiple chains use the multiclass QD-AMVA with
 * class-dependent scaling ({@code Pfqn_marie.pfqn_marie_multi}); single-server
 * queueing stations only.
 *
 * Ported from matlab/src/solvers/MVA/solver_mva_marie_analyzer.m.
 */
public final class Solver_mva_marie_analyzer {

    private Solver_mva_marie_analyzer() {}

    public static MVAResult solver_mva_marie_analyzer(NetworkStruct sn, SolverOptions options) {
        long startTime = System.nanoTime();
        int M = sn.nstations;
        int C = sn.nchains;

        Ret.snGetDemands dem = SnGetDemandsChain.snGetDemandsChain(sn);
        Matrix Lchain = dem.Dchain;
        Matrix STchain = dem.STchain;
        Matrix Vchain = dem.Vchain;
        Matrix alpha = dem.alpha;
        Matrix Nchain = dem.Nchain;
        Matrix SCVchain = dem.SCVchain;

        // Gating: closed models only.
        boolean hasSource = false;
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Source) { hasSource = true; break; }
        }
        for (int c = 0; c < C; c++) {
            if (Double.isInfinite(Nchain.get(c))) hasSource = true;
        }
        if (hasSource) {
            throw new RuntimeException("The 'marie' method supports closed models only; this model has "
                    + "open classes. Use another SolverMVA method (e.g. 'default').");
        }

        // Delay/infinite-server vs queueing stations; scheduling support check.
        boolean[] isDelay = new boolean[M];
        for (int i = 0; i < M; i++) {
            SchedStrategy s = sn.sched.get(sn.stations.get(i));
            boolean delay = Double.isInfinite(sn.nservers.get(i)) || s == SchedStrategy.INF;
            isDelay[i] = delay;
            boolean ok = delay || s == SchedStrategy.FCFS || s == SchedStrategy.PS || s == SchedStrategy.LCFSPR;
            if (!ok) {
                throw new RuntimeException("The 'marie' method supports FCFS, PS, LCFSPR and Delay stations "
                        + "only; station " + (i + 1) + " has an unsupported scheduling strategy.");
            }
        }

        int Mq = 0;
        for (int i = 0; i < M; i++) if (!isDelay[i]) Mq++;
        int[] queueRows = new int[Mq];
        int qi = 0;
        for (int i = 0; i < M; i++) if (!isDelay[i]) queueRows[qi++] = i;

        // Per-chain think time from delay stations.
        Matrix Z = new Matrix(1, C);
        for (int c = 0; c < C; c++) {
            double z = 0.0;
            for (int i = 0; i < M; i++) if (isDelay[i]) z += Lchain.get(i, c);
            Z.set(0, c, z);
        }

        // Queueing-station demands and effective SCV (only FCFS is sensitive).
        Matrix L = new Matrix(Mq, C);
        Matrix SCV = new Matrix(Mq, C);
        Matrix nservers = new Matrix(Mq, 1);
        for (int j = 0; j < Mq; j++) {
            int i = queueRows[j];
            boolean fcfs = sn.sched.get(sn.stations.get(i)) == SchedStrategy.FCFS;
            for (int c = 0; c < C; c++) {
                L.set(j, c, Lchain.get(i, c));
                double scv = SCVchain.get(i, c);
                if (!fcfs) scv = 1.0;
                if (!Double.isFinite(scv) || scv <= 0) scv = 1.0;
                SCV.set(j, c, scv);
            }
            double ns = sn.nservers.get(i);
            if (!Double.isFinite(ns)) ns = 1.0;
            nservers.set(j, 0, ns);
        }

        Matrix Xchain = new Matrix(1, C);
        Matrix Qchain = new Matrix(M, C);
        Matrix Uchain = new Matrix(M, C);
        Matrix Rchain = new Matrix(M, C);
        Matrix Tchain = new Matrix(M, C);
        int lastiter;

        if (Mq == 0) {
            // Nothing to isolate: with every station an infinite server the
            // aggregation-decomposition degenerates to the exact delay solution
            // X_c = N_c / Z_c, and Pfqn_marie would be handed a zero-row demand
            // matrix.
            lastiter = 1;
            for (int c = 0; c < C; c++) {
                double z = Z.get(0, c);
                double x = (z > 0) ? Nchain.get(c) / z : 0.0;
                Xchain.set(0, c, x);
                for (int i = 0; i < M; i++) Tchain.set(i, c, x * Vchain.get(i, c));
            }
        } else if (C == 1) {
            double N = Nchain.get(0);
            Pfqn_marie.Result mr = Pfqn_marie.pfqn_marie(L, N, Z.get(0, 0), SCV, 1e-8, 1000, nservers);
            lastiter = mr.iter;
            Xchain.set(0, 0, mr.X);
            for (int j = 0; j < Mq; j++) {
                int i = queueRows[j];
                Qchain.set(i, 0, mr.Q.get(j));
                Uchain.set(i, 0, mr.U.get(j));
            }
            for (int i = 0; i < M; i++) Tchain.set(i, 0, mr.X * Vchain.get(i, 0));
        } else {
            // Multiclass: QD-AMVA + cd-scaling. Multiserver queueing stations are
            // not supported by the multiclass isolation sub-model.
            for (int j = 0; j < Mq; j++) {
                if (nservers.get(j, 0) > 1) {
                    throw new UnsupportedOperationException("Multiclass 'marie' supports single-server "
                            + "queueing stations only; station has " + (int) nservers.get(j, 0) + " servers.");
                }
            }
            Matrix Nrow = new Matrix(1, C); for (int c = 0; c < C; c++) Nrow.set(0, c, Nchain.get(c));
            Pfqn_marie.MultiResult mr = Pfqn_marie.pfqn_marie_multi(L, Nrow, Z, SCV, 1e-8, 1000);
            lastiter = mr.iter;
            for (int c = 0; c < C; c++) Xchain.set(0, c, mr.X.get(0, c));
            for (int j = 0; j < Mq; j++) {
                int i = queueRows[j];
                for (int c = 0; c < C; c++) {
                    Qchain.set(i, c, mr.Q.get(j, c));
                    Uchain.set(i, c, mr.U.get(j, c));
                }
            }
            for (int i = 0; i < M; i++)
                for (int c = 0; c < C; c++) Tchain.set(i, c, mr.X.get(0, c) * Vchain.get(i, c));
        }

        for (int i = 0; i < M; i++) {
            for (int c = 0; c < C; c++) {
                double t = Tchain.get(i, c);
                if (t > 0) Rchain.set(i, c, Qchain.get(i, c) / t);
            }
        }
        // Delay stations: Q = U = T*S, R = S.
        for (int i = 0; i < M; i++) {
            if (isDelay[i]) {
                for (int c = 0; c < C; c++) {
                    double ts = Tchain.get(i, c) * STchain.get(i, c);
                    Qchain.set(i, c, ts);
                    Uchain.set(i, c, ts);
                    Rchain.set(i, c, STchain.get(i, c));
                }
            }
        }

        Ret.snDeaggregateChainResults dre = SnDeaggregateChainResults.snDeaggregateChainResults(
                sn, Lchain, null, STchain, Vchain, alpha, Qchain, Uchain, Rchain, Tchain, null, Xchain);

        MVAResult res = new MVAResult();
        res.QN = dre.Q;
        res.UN = dre.U;
        res.RN = dre.R;
        res.TN = dre.T;
        res.CN = dre.C;
        res.XN = dre.X;
        res.AN = new Matrix(0, 0);
        res.WN = new Matrix(0, 0);
        res.logNormConstAggr = Double.NaN;
        res.runtime = (System.nanoTime() - startTime) / 1000000000.0;
        res.iter = lastiter;
        res.method = "marie";
        return res;
    }
}
