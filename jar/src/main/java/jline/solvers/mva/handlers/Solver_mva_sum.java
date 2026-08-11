package jline.solvers.mva.handlers;

import java.util.ArrayList;
import java.util.List;

import jline.api.sn.SnDeaggregateChainResults;
import jline.api.sn.SnGetDemandsChain;
import jline.api.sum.Sum_closed;
import jline.api.sum.Sum_closing;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.util.matrix.Matrix;

/**
 * Summation method (SUM/ESUM) analyzer handler. Closed models are solved
 * with Sum_closed; open and mixed models with Sum_closing (closing method,
 * Kclosed=5000). FCFS and SIRO stations use their service SCVs (ESUM
 * corrections for scv!=1); PS, LCFS-PR and infinite-server stations are
 * insensitive and are passed scv=1. See jline.api.sum and Bolch et al.,
 * Secs. 9.2, 10.1.4.4 and 10.1.5. Mirrors MATLAB solver_mva_sum.m.
 */
public final class Solver_mva_sum {
    private Solver_mva_sum() {}

    public static MVAResult solver_mva_sum(NetworkStruct sn, SolverOptions options) {
        Ret.snGetDemands ret = SnGetDemandsChain.snGetDemandsChain(sn);
        Matrix Lchain = ret.Dchain;
        Matrix STchain = ret.STchain;
        Matrix Vchain = ret.Vchain;
        Matrix alpha = ret.alpha;
        Matrix Nchain = ret.Nchain;
        Matrix SCVchain = ret.SCVchain;
        Matrix refstatchain = ret.refstatchain;
        int M = sn.nstations;
        int K = sn.nchains;

        // station rows passed to the summation method (all but the source)
        List<Integer> rows = new ArrayList<Integer>();
        List<Double> miList = new ArrayList<Double>();
        for (int ist = 0; ist < M; ist++) {
            SchedStrategy s = sn.sched.get(sn.stations.get(ist));
            if (s == SchedStrategy.EXT) {
                // external world handled by lambda/sum_closing
            } else if (s == SchedStrategy.INF) {
                rows.add(ist);
                miList.add(Double.POSITIVE_INFINITY);
            } else if (s == SchedStrategy.PS || s == SchedStrategy.LCFSPR
                    || s == SchedStrategy.FCFS || s == SchedStrategy.SIRO) {
                rows.add(ist);
                miList.add(sn.nservers.get(ist));
            } else {
                throw new RuntimeException("The summation method does not support "
                        + SchedStrategy.toText(s) + " scheduling.");
            }
        }

        int Mr = rows.size();
        Matrix L = new Matrix(Mr, K);
        Matrix scv = new Matrix(Mr, K);
        Matrix mi = new Matrix(Mr, 1);
        for (int j = 0; j < Mr; j++) {
            int ist = rows.get(j);
            SchedStrategy s = sn.sched.get(sn.stations.get(ist));
            boolean scvSensitive = (s == SchedStrategy.FCFS || s == SchedStrategy.SIRO);
            mi.set(j, 0, miList.get(j));
            for (int r = 0; r < K; r++) {
                L.set(j, r, STchain.get(ist, r) * Vchain.get(ist, r));
                double c2 = 1;
                if (scvSensitive && Double.isFinite(SCVchain.get(ist, r)) && SCVchain.get(ist, r) > 0) {
                    c2 = SCVchain.get(ist, r);
                }
                scv.set(j, r, c2);
            }
        }

        Matrix Z = new Matrix(1, K);
        boolean hasOpen = false;
        for (int r = 0; r < K; r++) {
            if (Double.isInfinite(Nchain.get(r))) {
                hasOpen = true;
            }
        }

        Matrix Xchain;
        Matrix Qrows;
        Matrix Urows;
        int iter;
        if (!hasOpen) {
            Sum_closed.Result sres = Sum_closed.sum_closed(L, Nchain, Z, mi, scv,
                    options.iter_tol, (int) options.iter_max);
            Xchain = sres.XN;
            Qrows = sres.QN;
            Urows = sres.UN;
            iter = sres.it;
        } else {
            Matrix lambda = new Matrix(1, K);
            Matrix scva = new Matrix(1, K);
            for (int r = 0; r < K; r++) {
                scva.set(0, r, 1);
                if (Double.isInfinite(Nchain.get(r))) {
                    int refstat = (int) refstatchain.get(r);
                    lambda.set(0, r, 1 / STchain.get(refstat, r));
                    if (Double.isFinite(SCVchain.get(refstat, r)) && SCVchain.get(refstat, r) > 0) {
                        scva.set(0, r, SCVchain.get(refstat, r)); // interarrival SCV at the source
                    }
                }
            }
            Sum_closing.Result cres = Sum_closing.sum_closing(lambda, scva, L, mi, scv,
                    Nchain, Z, 5000, options.iter_tol, (int) options.iter_max);
            Xchain = cres.XN;
            Qrows = cres.QN;
            Urows = cres.UN;
            iter = cres.it;
        }

        Matrix Qchain = new Matrix(M, K);
        Matrix Uchain = new Matrix(M, K);
        Matrix Tchain = new Matrix(M, K);
        for (int j = 0; j < Mr; j++) {
            int ist = rows.get(j);
            for (int r = 0; r < K; r++) {
                Qchain.set(ist, r, Qrows.get(j, r));
                Uchain.set(ist, r, Urows.get(j, r));
            }
        }
        Matrix Rchain = new Matrix(M, K);
        for (int k = 0; k < M; k++) {
            for (int r = 0; r < K; r++) {
                double t = Xchain.get(r) * Vchain.get(k, r);
                Tchain.set(k, r, t);
                double w = (t > 0) ? Qchain.get(k, r) / t : 0;
                Rchain.set(k, r, Double.isFinite(w) ? w : 0);
            }
        }
        for (int r = 0; r < K; r++) {
            if (!Double.isFinite(Xchain.get(r)) || Nchain.get(r) == 0) {
                Xchain.set(0, r, 0);
                for (int k = 0; k < M; k++) {
                    Qchain.set(k, r, 0);
                    Uchain.set(k, r, 0);
                    Rchain.set(k, r, 0);
                    Tchain.set(k, r, 0);
                }
            }
        }

        Ret.snDeaggregateChainResults ret2 = SnDeaggregateChainResults.snDeaggregateChainResults(
                sn, Lchain, null, STchain, Vchain, alpha, null, null, Rchain, Tchain, null, Xchain);

        MVAResult res = new MVAResult();
        res.QN = ret2.Q;
        res.UN = ret2.U;
        res.RN = ret2.R;
        res.TN = ret2.T;
        res.CN = ret2.C;
        res.XN = ret2.X;
        res.logNormConstAggr = Double.NaN;
        res.iter = iter;
        return res;
    }
}
