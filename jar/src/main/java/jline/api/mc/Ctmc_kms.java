/**
 * @file Koury-McAllister-Stewart aggregation-disaggregation method for CTMCs
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import java.util.ArrayList;
import java.util.List;

import jline.util.Pair;
import jline.util.Triple;

import jline.util.matrix.Matrix;

public final class Ctmc_kms {
    private Ctmc_kms() {}

    /**
     * Koury-McAllister-Stewart aggregation-disaggregation method for CTMCs
     */
    public static Triple<Matrix, Double, Double> ctmc_kms(Matrix Q, List<List<Integer>> MS, int numSteps) {
        int nMacroStates = MS.size();
        int nStates = Q.getNumRows();

        // Start from Courtois decomposition solution (q derived as 1.05*max|Qperm|)
        Triple<Matrix, Double, Double> courtoisResult = Ctmc_courtois.ctmc_courtois(Q, MS);
        Matrix pcourt = courtoisResult.getFirst();
        double eps = courtoisResult.getSecond().doubleValue();
        double epsMAX = courtoisResult.getThird().doubleValue();

        // see _kb/03-api-layer.md for rationale
        double qVal = 1.05 * Q.elementMaxAbs();
        Pair<Matrix, Double> randResult = Ctmc_randomization.ctmc_randomization(Q, qVal);
        Matrix P = randResult.getFirst();

        // Rearrange P according to macrostates
        List<Integer> v = new ArrayList<Integer>();
        for (List<Integer> macroState : MS) {
            v.addAll(macroState);
        }
        Matrix Pperm = Matrix.zeros(nStates, nStates);
        for (int i = 0; i < nStates; i++) {
            for (int j = 0; j < nStates; j++) {
                Pperm.set(i, j, P.get(v.get(i), v.get(j)));
            }
        }

        // see _kb/03-api-layer.md for rationale
        Matrix pn = Matrix.zeros(1, nStates);
        for (int i = 0; i < nStates; i++) {
            pn.set(0, i, pcourt.get(v.get(i), 0));
        }

        // Main loop
        for (int n = 0; n < numSteps; n++) {
            Matrix pn_1 = pn.transpose();

            // Aggregation step
            Matrix pcondn_1 = Matrix.zeros(1, nStates);
            int procRows = 0;
            for (int I = 0; I < nMacroStates; I++) {
                int blockSize = MS.get(I).size();
                double sum = 0.0;
                for (int j = 0; j < blockSize; j++) {
                    sum += pn_1.get(procRows + j);
                }
                if (sum > 0) {
                    for (int j = 0; j < blockSize; j++) {
                        pcondn_1.set(0, procRows + j, pn_1.get(procRows + j) / sum);
                    }
                }
                procRows += blockSize;
            }

            // Compute aggregated transition matrix G
            Matrix G = Matrix.zeros(nMacroStates, nMacroStates);
            int procCols = 0;
            for (int I = 0; I < nMacroStates; I++) {
                procRows = 0;
                for (int J = 0; J < nMacroStates; J++) {
                    double sum = 0.0;
                    for (int j = 0; j < MS.get(J).size(); j++) {
                        for (int i = 0; i < MS.get(I).size(); i++) {
                            sum += pcondn_1.get(0, procRows + j) * Pperm.get(procRows + j, procCols + i);
                        }
                    }
                    G.set(I, J, sum);
                    procRows += MS.get(J).size();
                }
                procCols += MS.get(I).size();
            }

            Matrix w = Dtmc_solve.dtmc_solve(G.transpose());

            // Disaggregation step
            Matrix z = pcondn_1.copy();
            Matrix zn = Matrix.zeros(1, nStates);
            Matrix L = Matrix.zeros(nStates, nStates);
            Matrix D = Matrix.zeros(nStates, nStates);
            Matrix U = Matrix.zeros(nStates, nStates);

            procRows = 0;
            for (int I = 0; I < nMacroStates; I++) {
                int blockSizeI = MS.get(I).size();
                for (int j = 0; j < blockSizeI; j++) {
                    zn.set(0, procRows + j, w.get(I) * z.get(0, procRows + j));
                }

                procCols = 0;
                for (int J = 0; J < nMacroStates; J++) {
                    int blockSizeJ = MS.get(J).size();
                    if (I > J) {
                        for (int i = 0; i < blockSizeI; i++) {
                            for (int j = 0; j < blockSizeJ; j++) {
                                L.set(procRows + i, procCols + j, Pperm.get(procRows + i, procCols + j));
                            }
                        }
                    } else if (I == J) {
                        for (int i = 0; i < blockSizeI; i++) {
                            for (int j = 0; j < blockSizeJ; j++) {
                                if (i == j) {
                                    D.set(procRows + i, procCols + j, 1.0 - Pperm.get(procRows + i, procCols + j));
                                } else {
                                    D.set(procRows + i, procCols + j, -Pperm.get(procRows + i, procCols + j));
                                }
                            }
                        }
                    } else {
                        // I < J
                        for (int i = 0; i < blockSizeI; i++) {
                            for (int j = 0; j < blockSizeJ; j++) {
                                U.set(procRows + i, procCols + j, Pperm.get(procRows + i, procCols + j));
                            }
                        }
                    }
                    procCols += blockSizeJ;
                }
                procRows += blockSizeI;
            }

            // see _kb/03-api-layer.md for rationale
            Matrix M = D.sub(U);
            Matrix rhs = zn.mult(L).transpose();
            Matrix x = Matrix.zeros(nStates, 1);
            // Above the dispatch threshold the direct factorization is what limits the
            // model rather than the arithmetic; it remains the fallback.
            boolean solved = false;
            if (nStates > Ctmc_solve.GMRES_MIN_STATES) {
                Ctmc_gmres.GmresResult g =
                        Ctmc_gmres.ctmc_gmres(M.transpose(), rhs, 0.0, 0, 0, null);
                if (g.flag == 0) {
                    x = g.x;
                    solved = true;
                } else {
                    // Short-recurrence retry before the cubic factorization, as in Ctmc_solve.
                    Ctmc_bicgstab.BicgstabResult bs =
                            Ctmc_bicgstab.ctmc_bicgstab(M.transpose(), rhs, 0.0, 0, null);
                    if (bs.flag == 0) {
                        x = bs.x;
                        solved = true;
                    }
                }
            }
            if (!solved) {
                Matrix.solve(M.transpose(), rhs, x);
            }
            pn = x.transpose();

            // Normalize for numerical stability of the fixed point
            double pnSum = 0.0;
            for (int i = 0; i < nStates; i++) {
                pnSum += pn.get(0, i);
            }
            if (pnSum != 0.0) {
                for (int i = 0; i < nStates; i++) {
                    pn.set(0, i, pn.get(0, i) / pnSum);
                }
            }
        }

        // Reorder back to original ordering
        Matrix pout = Matrix.zeros(nStates, 1);
        for (int i = 0; i < nStates; i++) {
            pout.set(v.get(i), 0, pn.get(0, i));
        }

        return new Triple<Matrix, Double, Double>(pout, eps, epsMAX);
    }
}
