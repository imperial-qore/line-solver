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

public final class Solver_nc_jointaggr_ld {
    private Solver_nc_jointaggr_ld() {}

    /**
     * Computes load-dependent aggregated joint probabilities for the NC solver.
     */
    public static SolverNC.SolverNCJointReturn solver_nc_jointaggr_ld(NetworkStruct sn, SolverOptions options) {
        long startTimeMillis = System.nanoTime();

        // Initialization
        Matrix S = sn.nservers;
        Matrix rates = sn.rates;

        // Determine service times
        Matrix ST = rates.copy();
        for (int i = 0; i < ST.getNumRows(); i++) {
            for (int j = 0; j < ST.getNumCols(); j++) {
                ST.set(i, j, 1.0 / ST.get(i, j));
            }
        }
        ST.removeNaN();

        Ret.snGetDemands ret = SnGetDemandsChain.snGetDemandsChain(sn);
        Matrix Lchain = ret.Dchain;
        Matrix STchain = ret.STchain;
        Matrix alpha = ret.alpha;
        Matrix Nchain = ret.Nchain;

        int M = STchain.getNumRows();
        int K = STchain.getNumCols();

        // Build load-dependent scaling matrix mu_chain
        Matrix mu_chain = new Matrix(M, (int) Nchain.elementSum());
        for (int ist = 0; ist < M; ist++) {
            if (Utils.isInf(S.get(ist))) {
                for (int j = 0; j < (int) Nchain.elementSum(); j++) {
                    mu_chain.set(ist, j, (double) (j + 1));
                }
            } else {
                for (int j = 0; j < (int) Nchain.elementSum(); j++) {
                    mu_chain.set(ist, j, FastMath.min((double) (j + 1), S.get(ist)));
                }
            }
        }

        // Compute global normalization constant
        Matrix Zchain = Nchain.copy();
        Zchain.fill(0.0);
        double lG = Pfqn_ncld.pfqn_ncld(Lchain, Nchain, Zchain, mu_chain, options).lG;

        // Compute joint probability in log space
        double lPr = 0.0;
        for (int ist = 0; ist < M; ist++) {
            int isf = (int) sn.stationToStateful.get(ist);
            State.StateMarginalStatistics marginalResult = ToMarginal.toMarginal(sn, ist, sn.state.get(sn.stateful.get(isf)), null, null, null, null, null);
            Matrix nivec = marginalResult.nir;

            // Convert to chain populations
            Matrix nivec_chain = nivec.mult(sn.chains.transpose());

            // Extract parameters for station ist
            Matrix Lchain_ist = Matrix.extractRows(Lchain, ist, ist + 1, null);
            Matrix mu_chain_ist = Matrix.extractRows(mu_chain, ist, ist + 1, null);
            Matrix Znivec_chain = nivec_chain.copy();
            Znivec_chain.fill(0.0);

            // Compute lF_i: normalization constant for Lchain at this station
            double lF_i = Pfqn_ncld.pfqn_ncld(Lchain_ist, nivec_chain, Znivec_chain, mu_chain_ist, options).lG;

            // Build alpha-scaled service times for station ist
            Matrix ST_alpha_ist = new Matrix(1, sn.nclasses);
            for (int k = 0; k < sn.nclasses; k++) {
                ST_alpha_ist.set(0, k, ST.get(ist, k) * alpha.get(ist, k));
            }
            Matrix Znivec = nivec.copy();
            Znivec.fill(0.0);

            // Compute lg0_i: normalization constant for alpha-scaled service times
            double lg0_i = Pfqn_ncld.pfqn_ncld(ST_alpha_ist, nivec, Znivec, mu_chain_ist, options).lG;

            // Extract STchain for station ist
            Matrix STchain_ist = Matrix.extractRows(STchain, ist, ist + 1, null);

            // Compute lG0_i: normalization constant for STchain
            double lG0_i = Pfqn_ncld.pfqn_ncld(STchain_ist, nivec_chain, Znivec_chain, mu_chain_ist, options).lG;

            // Accumulate log probability: lF_i + (lg0_i - lG0_i)
            lPr += lF_i + (lg0_i - lG0_i);
        }

        // Final probability computation
        double Pr = FastMath.exp(lPr - lG);
        double G = FastMath.exp(lG);

        long endTimeMillis = System.nanoTime();
        double runtime = (endTimeMillis - startTimeMillis) / 1000000000.0;

        return new SolverNC.SolverNCJointReturn(Pr, G, lG, runtime);
    }
}
