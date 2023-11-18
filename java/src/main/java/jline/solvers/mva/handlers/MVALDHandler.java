package jline.solvers.mva.handlers;

import jline.api.PFQN;
import jline.api.SN;
import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVAResult;
import jline.util.Matrix;

/**
 * Handler for the solver_mvald function.
 */
public class MVALDHandler implements MVASolverHandler {

    @Override
    public SolverMVAResult solve(NetworkStruct sn, SolverOptions options) {
        SN.snGetDemandsChainReturn chainReturn = SN.snGetDemandsChain(sn);
        Matrix Lchain = chainReturn.Lchain;
        Matrix STchain = chainReturn.STchain;
        Matrix Vchain = chainReturn.Vchain;
        Matrix alpha = chainReturn.alpha;
        Matrix Nchain = chainReturn.Nchain;
        Matrix ST = new Matrix(sn.rates.getNumRows(), sn.rates.getNumCols());
        for(int i = 0; i < ST.getNumRows(); i++){
            for(int j = 0; j < ST.getNumCols(); j++){
                double val = 1 / sn.rates.get(i, j);
                if(!Double.isNaN(val)){
                   ST.set(i, j, val);
                }
            }
        }
        int M = STchain.getNumRows();
        int C = sn.nchains;
        Matrix S = sn.nservers;
        Matrix N = sn.njobs;
        double NchainSum = 0;
        for(int i = 0; i < Nchain.getNumCols(); i++){
            if(Double.isFinite(Nchain.get(i))){
                NchainSum += Nchain.get(i);
            }
        }
        Matrix mu_chain = Matrix.ones(M, (int) NchainSum);
        for(int i = 0; i < M; i++){
            if(Double.isInfinite(S.get(i))){
                for(int j = 0; j < NchainSum; j++){
                    mu_chain.set(i, j, j + 1);
                }
            } else if(!sn.lldscaling.isEmpty()){
                for(int j = 0; j < NchainSum; j++){
                    mu_chain.set(i, j, sn.lldscaling.get(i, j));
                }
            } /* else {
            mu_chain(i,1:sum(Nchain)) = ones(1,sum(Nchain(isfinite(Nchain))));
            } is redundant, mu_chain is already set to Matrix.ones
            */
        }
        Matrix lambda = new Matrix(1, C);
        for(int c = 0; c < sn.nchains; c++){
            for(int r = 0; r < N.getNumCols(); r++){
                if(Double.isInfinite(N.get(r)) && sn.chains.get(c, r) == 1){
                    Nchain.set(0, c, Double.POSITIVE_INFINITY);
                    lambda.set(0, c, lambda.get(0, c) + 1 / ST.get((int) sn.refstat.get(r), r));
                }
            }
        }
        PFQN.pfqnMVALDMXReturn ret = PFQN.pfqn_mvaldmx(lambda, Lchain, Nchain,
                new Matrix(Nchain.getNumRows(), Nchain.getNumCols()), mu_chain, S);
        Matrix Xchain = ret.XN;
        Matrix Qchain = ret.QN;
        Matrix Uchain = ret.UN;
        Matrix Tchain = Xchain.repmat(M, 1);
        for(int i = 0; i < Tchain.getNumRows(); i++){
            for(int j = 0; j < Tchain.getNumCols(); j++){
                Tchain.set(i, j, Tchain.get(i, j) * Vchain.get(i, j));
            }
        }
        Matrix Rchain = Xchain.repmat(M, 1);
        for(int i = 0; i < Rchain.getNumRows(); i++){
            for(int j = 0; j < Rchain.getNumCols(); j++){
                Rchain.set(i, j, Qchain.get(i, j) / Rchain.get(i, j));
            }
        }
        double lG = Double.NaN;

        // This is likely wrong as it uses Little's law for the utilization computation
        SN.snDeaggregateChainResultsReturn deAggregateReturn = SN.snDeaggregateChainResults(sn, Lchain,
                null, STchain, Vchain, alpha,  null, Uchain, Rchain,
                Tchain, null, Xchain);
        int iter = 1;
        SolverMVAResult res = new SolverMVAResult();
        res.QN = deAggregateReturn.Q;
        res.UN = deAggregateReturn.U;
        res.RN = deAggregateReturn.R;
        res.TN = deAggregateReturn.T;
        res.CN = deAggregateReturn.C;
        res.XN = deAggregateReturn.X;
        res.logNormConstAggr = lG;
        res.iter = iter;
        return res;
    }
}
