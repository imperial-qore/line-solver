package jline.solvers.nc.handlers;

import org.apache.commons.math3.util.FastMath;

import jline.api.mam.Map_lambda;
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
import jline.util.matrix.MatrixCell;

public final class Solver_nc_margaggr {
    private Solver_nc_margaggr() {}

    /**
     * Computes aggregated marginal probabilities for the NC solver.
     */
    public static SolverNC.SolverNCMargReturn solver_nc_margaggr(NetworkStruct sn, SolverOptions options, Double lG) {
        long startTimeMillis = System.nanoTime();

        int M = sn.nstations;
        int K = sn.nclasses;
        java.util.Map<jline.lang.nodes.StatefulNode, Matrix> state = sn.state;
        Matrix S = sn.nservers;
        Matrix NK = sn.njobs.transpose();
        int C = sn.nchains;
        java.util.Map<jline.lang.nodes.Station, java.util.Map<jline.lang.JobClass, MatrixCell>> PH = sn.proc;

        Matrix ST = new Matrix(M, K);
        for (int k = 0; k < K; k++) {
            for (int ist = 0; ist < M; ist++) {
                MatrixCell processCell = PH.get(sn.stations.get(ist)).get(sn.jobclasses.get(k));
                double lambda = Map_lambda.map_lambda(processCell);
                ST.set(ist, k, 1.0 / lambda);
            }
        }
        ST.removeNaN();

        Ret.snGetDemands ret = SnGetDemandsChain.snGetDemandsChain(sn);
        Matrix Lchain = ret.Dchain;
        Matrix STchain = ret.STchain;
        Matrix Nchain = ret.Nchain;

        Matrix V = new Matrix(sn.nstations, sn.nclasses);
        for (int c = 0; c < sn.nchains; c++) {
            Matrix inchain = sn.inchain.get(Integer.valueOf(c));
            for (int ist = 0; ist < sn.nstations; ist++) {
                for (int col = 0; col < inchain.getNumCols(); col++) {
                    int k = (int) inchain.get(0, col);
                    Matrix vc = sn.visits.get(Integer.valueOf(c));
                    V.set(ist, k, vc != null ? vc.get(ist, k) : 0.0);
                }
            }
        }

        int M_chain = STchain.getNumRows();

        Matrix mu = new Matrix(M_chain, (int) Nchain.elementSum());
        for (int ist = 0; ist < M_chain; ist++) {
            if (Utils.isInf(S.get(ist))) {
                for (int j = 0; j < (int) Nchain.elementSum(); j++) {
                    mu.set(ist, j, (double) (j + 1));
                }
            } else {
                for (int j = 0; j < (int) Nchain.elementSum(); j++) {
                    mu.set(ist, j, FastMath.min((double) (j + 1), S.get(ist)));
                }
            }
        }

        double lG_val;
        if (lG != null) {
            lG_val = lG.doubleValue();
        } else {
            Matrix Zchain_tmp = Nchain.copy();
            Zchain_tmp.fill(0.0);
            lG_val = (Double) Pfqn_ncld.pfqn_ncld(Lchain, Nchain, Zchain_tmp, mu, options).lG;
        }
        double G = FastMath.exp(lG_val);

        Matrix lPr = new Matrix(sn.nstations, 1);
        for (int ist = 0; ist < sn.nstations; ist++) {
            int ind = (int) sn.stationToNode.get(ist);
            int isf = (int) sn.stationToStateful.get(ist);
            State.StateMarginalStatistics marginalResult =
                    ToMarginal.toMarginal(sn, ind, state.get(sn.stateful.get(isf)), null, null, null, null, null);
            Matrix nivec = marginalResult.nir;

            if (nivec.elementMin() < 0) {
                lPr.set(ist, Double.NaN);
            } else {
                Matrix Lchain_minus_i = new Matrix(0, Lchain.getNumCols());
                Matrix mu_minus_i = new Matrix(0, mu.getNumCols());
                for (int j = 0; j < sn.nstations; j++) {
                    if (j != ist) {
                        Matrix Lchain_row = Matrix.extractRows(Lchain, j, j + 1, null);
                        Matrix mu_row = Matrix.extractRows(mu, j, j + 1, null);
                        Lchain_minus_i = Matrix.concatRows(Lchain_minus_i, Lchain_row, null);
                        mu_minus_i = Matrix.concatRows(mu_minus_i, mu_row, null);
                    }
                }

                Matrix nivec_chain = nivec.mult(sn.chains.transpose());
                Matrix Nchain_minus_i = Nchain.copy();
                for (int j = 0; j < Nchain_minus_i.length(); j++) {
                    Nchain_minus_i.set(j, Nchain_minus_i.get(j) - nivec_chain.get(j));
                }
                Matrix Zchain_minus_i = Nchain.copy();
                Zchain_minus_i.fill(0.0);

                double lG_minus_i = (Double) Pfqn_ncld.pfqn_ncld(Lchain_minus_i, Nchain_minus_i, Zchain_minus_i, mu_minus_i, options).lG;

                Matrix ST_ist = new Matrix(1, K);
                Matrix V_ist = new Matrix(1, K);
                for (int k = 0; k < K; k++) {
                    ST_ist.set(0, k, ST.get(ist, k));
                    V_ist.set(0, k, V.get(ist, k));
                }
                Matrix ST_V_ist = ST_ist.elementMult(V_ist);
                Matrix mu_ist = Matrix.extractRows(mu, ist, ist + 1, null);
                Matrix Znivec = nivec.copy();
                Znivec.fill(0.0);

                double lF_i = (Double) Pfqn_ncld.pfqn_ncld(ST_V_ist, nivec, Znivec, mu_ist, options).lG;

                lPr.set(ist, lF_i + lG_minus_i - lG_val);
            }
        }

        Matrix Pr = lPr.copy();
        for (int i = 0; i < Pr.getNumRows(); i++) {
            for (int j = 0; j < Pr.getNumCols(); j++) {
                Pr.set(i, j, FastMath.exp(Pr.get(i, j)));
            }
        }
        Pr.removeNaN();

        long endTimeMillis = System.nanoTime();
        double runtime = (double) (endTimeMillis - startTimeMillis) / 1000000000.0;

        return new SolverNC.SolverNCMargReturn(Pr, G, lG != null ? lG.doubleValue() : Double.MIN_VALUE, runtime);
    }
}
