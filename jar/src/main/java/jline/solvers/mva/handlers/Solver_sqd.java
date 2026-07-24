package jline.solvers.mva.handlers;

import jline.api.npfqn.Npfqn_sqd;
import jline.api.sn.SnDeaggregateChainResults;
import jline.api.sn.SnGetDemandsChain;
import jline.io.InputOutput;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.util.matrix.Matrix;

/**
 * Handler for the Blocking-After-Service (BAS) approximate MVA method.
 *
 * <p>Wraps the {@link Npfqn_sqd#npfqn_sqd(NetworkStruct, int)} chain-aggregated
 * approximation: solves the single closed chain, then disaggregates the chain-level
 * throughput, queue length and utilization back to per-class results via
 * {@link SnDeaggregateChainResults}. Supports single-chain closed networks (one or
 * more classes connected by class switching); multi-chain models warn and return empty, since
 * the BAS approximation models a single circulating population.</p>
 */
public final class Solver_sqd {
    private Solver_sqd() {}

    private static Matrix nanMatrix(int rows, int cols) {
        Matrix m = new Matrix(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                m.set(i, j, Double.NaN);
            }
        }
        return m;
    }

    public static MVAResult solver_sqd(NetworkStruct sn, SolverOptions options) {
        int C = sn.nchains;
        if (C != 1) {
            InputOutput.line_warning("solver_sqd",
                    "SQD (Smith Queue Decomposition) supports single-chain closed networks only; "
                    + "this model is multichain (nchains=" + C + ") — returning empty results.");
            int M0 = sn.nstations;
            int K0 = sn.nclasses;
            MVAResult empty = new MVAResult();
            empty.QN = nanMatrix(M0, K0);
            empty.UN = nanMatrix(M0, K0);
            empty.RN = nanMatrix(M0, K0);
            empty.TN = nanMatrix(M0, K0);
            empty.CN = nanMatrix(1, K0);
            empty.XN = nanMatrix(1, K0);
            empty.logNormConstAggr = Double.NaN;
            empty.iter = 0;
            return empty;
        }

        Ret.snGetDemands dem = SnGetDemandsChain.snGetDemandsChain(sn);
        Matrix Lchain  = dem.Dchain;
        Matrix STchain = dem.STchain;
        Matrix Vchain  = dem.Vchain;
        Matrix alpha   = dem.alpha;
        Matrix refstatchain = dem.refstatchain;

        int M = sn.nstations;
        int N = sn.nclosedjobs;

        Ret.pfqnMVA bas = Npfqn_sqd.npfqn_sqd(sn, N);

        // Assemble chain-level (M x 1) matrices for the single closed chain.
        int refstat = (int) refstatchain.get(0, 0);
        Matrix Xchain = new Matrix(1, 1);
        Matrix Tchain = new Matrix(M, 1);
        Matrix Qchain = new Matrix(M, 1);
        Matrix Uchain = new Matrix(M, 1);
        Matrix Rchain = new Matrix(M, 1);
        // Vchain is normalized to 1 at the reference station, so the per-station
        // throughput there equals the chain reference throughput.
        Xchain.set(0, 0, bas.X.get(refstat, 0));
        for (int i = 0; i < M; i++) {
            Tchain.set(i, 0, bas.X.get(i, 0));
            Qchain.set(i, 0, bas.Q.get(i, 0));
            Uchain.set(i, 0, bas.U.get(i, 0));
            Rchain.set(i, 0, bas.R.get(i, 0));
        }

        Ret.snDeaggregateChainResults dr = SnDeaggregateChainResults.snDeaggregateChainResults(
                sn, Lchain, null, STchain, Vchain, alpha, Qchain, Uchain, Rchain, Tchain, null, Xchain);

        MVAResult res = new MVAResult();
        res.QN = dr.Q;
        res.UN = dr.U;
        res.RN = dr.R;
        res.TN = dr.T;
        res.CN = dr.C;
        res.XN = dr.X;
        res.logNormConstAggr = Double.NaN;
        res.iter = 0;
        return res;
    }
}
