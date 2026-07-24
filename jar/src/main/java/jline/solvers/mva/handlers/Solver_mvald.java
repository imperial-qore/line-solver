/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.handlers;

import jline.api.pfqn.ld.Pfqn_mvaldmx;
import jline.api.sn.SnDeaggregateChainResults;
import jline.api.sn.SnGetDemandsChain;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.util.Utils;
import jline.util.matrix.Matrix;

public final class Solver_mvald {
    private Solver_mvald() {}

    /**
     * Handler for the solver_mvald function.
     */
    public static MVAResult solver_mvald(NetworkStruct sn, SolverOptions options) {
        Ret.snGetDemands chainReturn = SnGetDemandsChain.snGetDemandsChain(sn);
        Matrix Lchain = chainReturn.Dchain;
        Matrix STchain = chainReturn.STchain;
        Matrix Vchain = chainReturn.Vchain;
        Matrix alpha = chainReturn.alpha;
        Matrix Nchain = chainReturn.Nchain;
        Matrix ST = new Matrix(sn.rates.getNumRows(), sn.rates.getNumCols());
        for (int i = 0; i < ST.getNumRows(); i++) {
            for (int j = 0; j < ST.getNumCols(); j++) {
                double serviceTime = 1.0 / sn.rates.get(i, j);
                if (!Double.isNaN(serviceTime)) {
                    ST.set(i, j, serviceTime);
                }
            }
        }
        int M = STchain.getNumRows();
        int C = sn.nchains;
        Matrix S = sn.nservers;
        Matrix N = sn.njobs;
        double NchainSum = 0.0;
        for (int i = 0; i < Nchain.getNumCols(); i++) {
            if (Double.isFinite(Nchain.get(i))) {
                NchainSum += Nchain.get(i);
            }
        }
        Matrix mu_chain = Matrix.ones(M, (int) NchainSum);
        for (int i = 0; i < M; i++) {
            if (Utils.isInf(S.get(i))) {
                int j = 0;
                while (j < NchainSum) {
                    mu_chain.set(i, j, j + 1);
                    j++;
                }
            } else if (sn.lldscaling != null && !sn.lldscaling.isEmpty()) {
                int j = 0;
                while (j < NchainSum) {
                    mu_chain.set(i, j, sn.lldscaling.get(i, j));
                    j++;
                }
            }
        }
        Matrix lambda = new Matrix(1, C);
        for (int c = 0; c < sn.nchains; c++) {
            for (int r = 0; r < N.getNumCols(); r++) {
                if (Utils.isInf(N.get(r)) && sn.chains.get(c, r) == 1.0) {
                    Nchain.set(0, c, Double.POSITIVE_INFINITY);
                    lambda.set(0, c, lambda.get(0, c) + 1.0 / ST.get((int) sn.refstat.get(r), r));
                }
            }
        }
        Ret.pfqnMVALDMX ret = Pfqn_mvaldmx.pfqn_mvaldmx(lambda, Lchain, Nchain,
                new Matrix(Nchain.getNumRows(), Nchain.getNumCols()), mu_chain, S);
        Matrix Xchain = ret.X;
        Matrix Qchain = ret.Q;
        Matrix Uchain = ret.U;
        Matrix Tchain = Xchain.repmat(M, 1);
        for (int i = 0; i < Tchain.getNumRows(); i++) {
            for (int j = 0; j < Tchain.getNumCols(); j++) {
                Tchain.set(i, j, Tchain.get(i, j) * Vchain.get(i, j));
            }
        }
        Matrix Rchain = Xchain.repmat(M, 1);
        for (int i = 0; i < Rchain.getNumRows(); i++) {
            for (int j = 0; j < Rchain.getNumCols(); j++) {
                Rchain.set(i, j, Qchain.get(i, j) / Rchain.get(i, j));
            }
        }
        double lG = Double.NaN;

        // This is likely wrong as it uses Little's law for the utilization computation
        Ret.snDeaggregateChainResults deAggregateReturn =
                SnDeaggregateChainResults.snDeaggregateChainResults(sn, Lchain, null, STchain, Vchain, alpha, null,
                        Uchain, Rchain, Tchain, null, Xchain);
        int iter = 1;
        // see _kb/06-solver-catalog.md for rationale
        if (sn.lldscaling != null && !sn.lldscaling.isEmpty()) {
            Matrix U = deAggregateReturn.U;
            Matrix T = deAggregateReturn.T;
            for (int ist = 0; ist < Math.min(U.getNumRows(), sn.lldscaling.getNumRows()); ist++) {
                if (Double.isFinite(S.get(ist))) {
                    double ceff = S.get(ist);
                    for (int j = 0; j < sn.lldscaling.getNumCols(); j++) {
                        ceff = Math.max(ceff, sn.lldscaling.get(ist, j));
                    }
                    for (int r = 0; r < U.getNumCols(); r++) {
                        double stir = ST.get(ist, r);
                        if (Double.isFinite(stir) && stir > 0) {
                            U.set(ist, r, T.get(ist, r) * stir / ceff);
                        }
                    }
                }
            }
        }
        MVAResult res = new MVAResult();
        res.QN = deAggregateReturn.Q;
        res.UN = deAggregateReturn.U;
        res.RN = deAggregateReturn.R;
        res.TN = deAggregateReturn.T;
        res.CN = deAggregateReturn.C;
        res.XN = deAggregateReturn.X;
        res.logNormConstAggr = lG;
        res.iter = iter;
        res.method = options.method;
        return res;
    }
}
