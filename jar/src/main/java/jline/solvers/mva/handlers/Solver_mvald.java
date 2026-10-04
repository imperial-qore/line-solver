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
        Matrix refstatchain = chainReturn.refstatchain;
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
        double NchainSum = 0.0;
        for (int i = 0; i < Nchain.getNumCols(); i++) {
            if (Double.isFinite(Nchain.get(i))) {
                NchainSum += Nchain.get(i);
            }
        }
        int Nt = (int) NchainSum;

        // Chain arrival rates. The reference station of an open chain is its Source
        // and STchain holds one over the SUM of the class arrival rates there, so
        // reading the rate at chain level also covers a chain whose classes arrive at
        // several rates, or one carrying a class only ever reached by a switch (whose
        // own service time at the source is zero, making the per-class 1/ST infinite).
        Matrix lambda = new Matrix(1, C);
        boolean[] openChain = new boolean[C];
        int nOpenChains = 0;
        for (int c = 0; c < C; c++) {
            if (!Utils.isInf(Nchain.get(c))) {
                continue;
            }
            openChain[c] = true;
            nOpenChains++;
            int rst = (int) refstatchain.get(c, 0);
            if (STchain.get(rst, c) > 0) {
                lambda.set(0, c, 1.0 / STchain.get(rst, c));
            }
        }

        Matrix Xchain;
        Matrix Qchain;
        Matrix Uchain;
        if (nOpenChains == 0) {
            // PURELY CLOSED. Every station enters the recursion, an infinite server as
            // the load-dependent rate mu(n)=n, which is exact because n cannot then
            // exceed the closed population.
            Matrix mu_chain = Matrix.ones(M, Nt);
            for (int i = 0; i < M; i++) {
                if (Utils.isInf(S.get(i))) {
                    for (int j = 0; j < Nt; j++) {
                        mu_chain.set(i, j, j + 1);
                    }
                } else if (sn.lldscaling != null && !sn.lldscaling.isEmpty()) {
                    for (int j = 0; j < Nt; j++) {
                        mu_chain.set(i, j, sn.lldscaling.get(i, j));
                    }
                }
            }
            Ret.pfqnMVALDMX ret = Pfqn_mvaldmx.pfqn_mvaldmx(lambda, Lchain, Nchain,
                    new Matrix(Nchain.getNumRows(), Nchain.getNumCols()), mu_chain, S);
            Xchain = ret.X;
            Qchain = ret.Q;
            Uchain = ret.U;
        } else {
            // MIXED OR PURELY OPEN. Three kinds of row are not the same thing to
            // pfqn_mvaldmx and have to be separated before it is called. This is the
            // partition Solver_ncld makes for the same recursion.
            //  - THE SOURCE IS NOT A STATION. Its chain demand is the interarrival time
            //    1/lambda, so it carries offered load Lo=1 exactly, Pfqn_ldmx_ec then
            //    forms 1/(1-Lo/mu)=Inf and the 0*NaN in the residence-time sum turned
            //    every chain into NaN, which the deaggregation reported as zeros.
            //  - A DELAY IS AN INFINITE SERVER FOR THE OPEN CHAINS TOO. mu(n)=n cut at
            //    the closed population declares it saturated at Nt jobs. It enters as
            //    chain think time instead and its queue length is X*L, which is exact.
            //  - A QUEUEING STATION KEEPS ITS WHOLE RATE ROW. Pfqn_ldmx_ec reads the
            //    limited-load-dependence level b off the row itself, so a row cut at
            //    the closed population is read as a slower station, and with no closed
            //    class at all it collapses to mu(1), a single fixed-rate server.
            boolean[] sourceStation = new boolean[M];
            for (int c = 0; c < C; c++) {
                if (openChain[c]) {
                    sourceStation[(int) refstatchain.get(c, 0)] = true;
                }
            }
            boolean[] delayStation = new boolean[M];
            int nq = 0;
            for (int i = 0; i < M; i++) {
                if (sourceStation[i]) {
                    continue;
                }
                if (Utils.isInf(S.get(i))) {
                    delayStation[i] = true;
                } else {
                    nq++;
                }
            }
            int[] queueStations = new int[nq];
            int qi = 0;
            for (int i = 0; i < M; i++) {
                if (!sourceStation[i] && !delayStation[i]) {
                    queueStations[qi++] = i;
                }
            }
            Matrix Zchain = new Matrix(Nchain.getNumRows(), Nchain.getNumCols());
            for (int c = 0; c < C; c++) {
                double z = 0.0;
                for (int i = 0; i < M; i++) {
                    if (delayStation[i]) {
                        z += Lchain.get(i, c);
                    }
                }
                Zchain.set(0, c, z);
            }
            int lldWidth = (sn.lldscaling == null || sn.lldscaling.isEmpty()) ? 0 : sn.lldscaling.getNumCols();
            int ncol = Math.max(1, Nt);
            for (int k = 0; k < nq; k++) {
                // first column of the trailing constant run, the level b of Pfqn_ldmx_ec
                int b = lldWidth;
                while (b > 1 && sn.lldscaling.get(queueStations[k], b - 2)
                        == sn.lldscaling.get(queueStations[k], b - 1)) {
                    b--;
                }
                ncol = Math.max(ncol, b);
            }
            Matrix mu_chain = Matrix.ones(nq, ncol);
            for (int k = 0; k < nq; k++) {
                if (lldWidth == 0) {
                    continue;
                }
                for (int j = 0; j < ncol; j++) {
                    // saturated tail past the end of the caller's own row
                    mu_chain.set(k, j, sn.lldscaling.get(queueStations[k], Math.min(j, lldWidth - 1)));
                }
            }
            Matrix Dq = new Matrix(nq, C);
            Matrix Sq = new Matrix(nq, 1);
            for (int k = 0; k < nq; k++) {
                for (int c = 0; c < C; c++) {
                    Dq.set(k, c, Lchain.get(queueStations[k], c));
                }
                Sq.set(k, 0, S.get(queueStations[k]));
            }
            Ret.pfqnMVALDMX ret = Pfqn_mvaldmx.pfqn_mvaldmx(lambda, Dq, Nchain, Zchain, mu_chain, Sq);
            Xchain = ret.X;
            Qchain = new Matrix(M, C);
            Uchain = new Matrix(M, C);
            for (int k = 0; k < nq; k++) {
                for (int c = 0; c < C; c++) {
                    Qchain.set(queueStations[k], c, ret.Q.get(k, c));
                    Uchain.set(queueStations[k], c, ret.U.get(k, c));
                }
            }
            for (int i = 0; i < M; i++) {
                if (delayStation[i]) {
                    for (int c = 0; c < C; c++) {
                        // infinite server: X*L for a closed chain, lambda*L for an open one
                        Qchain.set(i, c, Lchain.get(i, c) * Xchain.get(c));
                    }
                }
            }
        }
        Matrix Tchain = Xchain.repmat(M, 1);
        for (int i = 0; i < Tchain.getNumRows(); i++) {
            for (int j = 0; j < Tchain.getNumCols(); j++) {
                Tchain.set(i, j, Tchain.get(i, j) * Vchain.get(i, j));
            }
        }
        // Rchain is the PER-VISIT response time Qchain./Tchain, not the residence time
        // Qchain./Xchain: snDeaggregateChainResults multiplies the visit ratio back in,
        // so dividing by Xchain here counts it twice and the class queue lengths stop
        // summing to the population (MATLAB solver_mvald.m carries the same note)
        Matrix Rchain = new Matrix(M, Tchain.getNumCols());
        for (int i = 0; i < Rchain.getNumRows(); i++) {
            for (int j = 0; j < Rchain.getNumCols(); j++) {
                double t = Tchain.get(i, j);
                double r = Qchain.get(i, j) / t;
                Rchain.set(i, j, Double.isFinite(r) ? r : 0.0);
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
