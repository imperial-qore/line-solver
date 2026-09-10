/**
 * @file NC solver: aggregated joint probability handler
 *
 * @since LINE 3.0
 */
package jline.solvers.nc.handlers;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.util.FastMath;

import jline.api.pfqn.ld.Pfqn_ncld;
import jline.api.sn.SnGetDemandsChain;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.state.ToMarginal;
import jline.solvers.SolverOptions;
import jline.solvers.nc.SolverNC;
import jline.util.Utils;
import jline.util.matrix.Matrix;

public final class Solver_nc_jointaggr {
    private Solver_nc_jointaggr() {}

    /**
     * Computes aggregated joint probabilities for the NC solver.
     */
    public static SolverNC.SolverNCJointReturn solver_nc_jointaggr(NetworkStruct sn, SolverOptions options) {
        long startTimeMillis = System.nanoTime();

        Matrix V = new Matrix(sn.nstateful, sn.nclasses);
        for (int i = 0; i < sn.visits.size(); i++) {
            V = V.add(1.0, sn.visits.get(i));
        }

        int M = sn.nstations;
        int C = sn.nchains;
        Matrix Nchain = new Matrix(1, C);

        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            double chainPop = 0.0;
            for (int col = 0; col < inchain.getNumCols(); col++) {
                int k = (int) inchain.get(0, col);
                chainPop += sn.njobs.get(k);
            }
            Nchain.set(c, chainPop);
        }

        Matrix nservers = sn.nservers;
        Matrix mu = new Matrix(M, (int) Nchain.elementSum());
        for (int ist = 0; ist < M; ist++) {
            if (Utils.isInf(nservers.get(ist))) {
                for (int j = 0; j < (int) Nchain.elementSum(); j++) {
                    mu.set(ist, j, (double) (j + 1));
                }
            } else {
                for (int j = 0; j < (int) Nchain.elementSum(); j++) {
                    mu.set(ist, j, FastMath.min((double) (j + 1), nservers.get(ist)));
                }
            }
        }
        Object state = sn.state;

        double lG;
        Matrix ST;

        if ("exact".equals(options.method)) {
            Ret.snGetDemands ret = SnGetDemandsChain.snGetDemandsChain(sn);
            Matrix Lchain = ret.Dchain;
            Matrix Nchain_ret = ret.Nchain;
            Matrix Zchain = Nchain_ret.copy();
            Zchain.fill(0.0);
            lG = (Double) Pfqn_ncld.pfqn_ncld(Lchain, Nchain_ret, Zchain, mu, options).lG;

            ST = sn.rates.copy();
            for (int i = 0; i < ST.getNumRows(); i++) {
                for (int j = 0; j < ST.getNumCols(); j++) {
                    ST.set(i, j, 1.0 / ST.get(i, j));
                }
            }
            ST.removeNaN();
        } else {
            jline.solvers.nc.SolverNC.SolverNCReturn ncResult = Solver_nc.solver_nc(sn, options);
            lG = ncResult.lG;
            ST = ncResult.STeff;
        }

        double G = FastMath.exp(lG);
        double lPr = 0.0;

        for (int ist = 0; ist < M; ist++) {
            int isf = (int) sn.stationToStateful.get(ist);
            int ind = (int) sn.stationToNode.get(ist);
            Matrix stateMatrix = (Matrix) ((java.util.Map<?, ?>) state).get(sn.stateful.get(isf));
            jline.lang.state.State.StateMarginalStatistics marginalResult = ToMarginal.toMarginal(sn, ind, stateMatrix, null, null, null, null, null);
            Matrix nivec = marginalResult.nir;

            Matrix uniqueNivec = getUniqueRows(nivec);

            for (int row = 0; row < uniqueNivec.getNumRows(); row++) {
                Matrix currentNivec = Matrix.extractRows(uniqueNivec, row, row + 1, null);
                Matrix nivec_chain = currentNivec.mult(sn.chains.transpose());

                boolean hasPositivePop = false;
                for (int j = 0; j < nivec_chain.length(); j++) {
                    if (nivec_chain.get(j) > 0) {
                        hasPositivePop = true;
                        break;
                    }
                }

                if (hasPositivePop) {
                    Matrix ST_ist = new Matrix(1, sn.nclasses);
                    Matrix V_ist = new Matrix(1, sn.nclasses);
                    for (int k = 0; k < sn.nclasses; k++) {
                        ST_ist.set(0, k, ST.get(ist, k));
                        V_ist.set(0, k, V.get(ist, k));
                    }
                    Matrix ST_V_ist = ST_ist.elementMult(V_ist);
                    Matrix mu_ist = Matrix.extractRows(mu, ist, ist + 1, null);
                    Matrix Znivec = currentNivec.copy();
                    Znivec.fill(0.0);

                    double lF_i = (Double) Pfqn_ncld.pfqn_ncld(ST_V_ist, currentNivec, Znivec, mu_ist, options).lG;
                    lPr += lF_i;
                }
            }
        }

        lPr -= lG;
        double Pr = FastMath.exp(lPr);

        long endTimeMillis = System.nanoTime();
        double runtime = (endTimeMillis - startTimeMillis) / 1000000000.0;

        return new SolverNC.SolverNCJointReturn(Pr, G, lG, runtime);
    }

    private static Matrix getUniqueRows(Matrix matrix) {
        List<Matrix> uniqueRowsList = new ArrayList<Matrix>();
        Set<String> seenRows = new HashSet<String>();

        for (int i = 0; i < matrix.getNumRows(); i++) {
            Matrix row = Matrix.extractRows(matrix, i, i + 1, null);
            StringBuilder sb = new StringBuilder();
            for (int c = 0; c < row.getNumCols(); c++) {
                if (c > 0) sb.append(',');
                sb.append(row.get(0, c));
            }
            String rowString = sb.toString();

            if (!seenRows.contains(rowString)) {
                seenRows.add(rowString);
                uniqueRowsList.add(row);
            }
        }

        if (uniqueRowsList.isEmpty()) {
            return new Matrix(0, matrix.getNumCols());
        }

        Matrix result = uniqueRowsList.get(0);
        for (int i = 1; i < uniqueRowsList.size(); i++) {
            result = Matrix.concatRows(result, uniqueRowsList.get(i), null);
        }

        return result;
    }
}
