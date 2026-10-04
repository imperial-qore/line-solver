package jline.solvers.nc.handlers;

import org.apache.commons.math3.util.FastMath;

import jline.api.pfqn.ld.Pfqn_ncld;
import jline.api.sn.SnGetDemandsChain;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.state.State;
import jline.lang.state.ToMarginal;
import jline.solvers.SolverOptions;
import jline.solvers.nc.SolverNC;
import jline.util.Utils;
import jline.util.matrix.Matrix;

public final class Solver_nc_joint {
    private Solver_nc_joint() {}

    public static SolverNC.SolverNCJointReturn solver_nc_joint(NetworkStruct sn, SolverOptions options) {
        Matrix S = sn.nservers;
        Matrix rates = sn.rates;
        Matrix ST = rates.copy();
        for (int i = 0; i < ST.getNumRows(); i++) {
            for (int j = 0; j < ST.getNumCols(); j++) {
                ST.set(i, j, 1.0 / ST.get(i, j));
            }
        }
        ST.removeNaN();
        Ret.snGetDemands ret = SnGetDemandsChain.snGetDemandsChain(sn);
        Matrix STchain = ret.STchain;
        Matrix Vchain = ret.Vchain;
        Matrix alpha = ret.alpha;
        Matrix Nchain = ret.Nchain;
        long startTimeMillis = System.nanoTime();
        int M = STchain.getNumRows();
        int K = ST.getNumCols();
        Matrix Lchain = new Matrix(0, K);

        Matrix mu_chain = new Matrix(0, (int) Nchain.elementSum());
        for (int i = 0; i < M; i++) {
            Matrix tmp = new Matrix(1, (int) Nchain.elementSum());
            if (Utils.isInf(S.get(i))) {
                for (int j = 0; j < tmp.length(); j++) {
                    tmp.set(j, (double) (j + 1));
                }
            } else {
                for (int j = 0; j < tmp.length(); j++) {
                    tmp.set(j, FastMath.min((double) (j + 1), S.get(i)));
                }
            }
            mu_chain = Matrix.concatRows(mu_chain, tmp, null);

            Matrix tmp1 = Matrix.extractRows(STchain, i, i + 1, null)
                    .elementMult(Matrix.extractRows(Vchain, i, i + 1, null), null);
            Lchain = Matrix.concatRows(Lchain, tmp1, null);
        }

        Matrix Z_tmp = Nchain.copy();
        Z_tmp.fill(0.0);
        double lG = Pfqn_ncld.pfqn_ncld(Lchain, Nchain, Z_tmp, mu_chain, options).lG;
        double lPr = 0.0;

        for (int i = 0; i < M; i++) {
            int isf = (int) sn.stationToStateful.get(i);
            State.StateMarginalStatistics ret1 = ToMarginal.toMarginal(
                    sn, i, sn.state.get(sn.stateful.get(isf)), null, null, null, null, null);
            Matrix nivec = ret1.nir;
            Matrix nivec_chain = nivec.mult(sn.chains.transpose());
            Matrix Lchain_i = Matrix.extractRows(Lchain, i, i + 1, null);
            Matrix ST_i = Matrix.extractRows(ST, i, i + 1, null);
            Matrix alpha_i = Matrix.extractRows(alpha, i, i + 1, null);
            Matrix STchain_i = Matrix.extractRows(STchain, i, i + 1, null);
            Matrix mu_chain_i = Matrix.extractRows(mu_chain, i, i + 1, null);
            Matrix Zvec_chain = nivec_chain.copy();
            Zvec_chain.fill(0.0);
            Matrix Zvec = nivec.copy();
            Zvec.fill(0.0);

            double lF_i = Pfqn_ncld.pfqn_ncld(Lchain_i, nivec_chain, Zvec_chain, mu_chain_i, options).lG;
            double lg0_i = Pfqn_ncld.pfqn_ncld(ST_i.elementMult(alpha_i, null), nivec, Zvec, mu_chain_i, options).lG;
            double lG0_i = Pfqn_ncld.pfqn_ncld(STchain_i, nivec_chain, Zvec_chain, mu_chain_i, options).lG;
            lPr = lPr + lF_i + (lg0_i - lG0_i);
        }
        double Pr = FastMath.exp(lPr - lG);
        long endTimeMillis = System.nanoTime();
        double runtime = (endTimeMillis - startTimeMillis) / 1000000000.0;
        double G = FastMath.exp(lG);
        return new SolverNC.SolverNCJointReturn(Pr, G, lG, runtime);
    }
}
