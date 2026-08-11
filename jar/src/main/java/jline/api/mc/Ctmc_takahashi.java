package jline.api.mc;

import jline.util.Pair;
import jline.util.Triple;

import jline.util.matrix.Matrix;

import java.util.List;

public final class Ctmc_takahashi {
    private Ctmc_takahashi() {}

    /**
     * Takahashi's aggregation-disaggregation method for CTMCs
     */
    public static Triple<Matrix, Double, Double> ctmc_takahashi(Matrix Q, List<List<Integer>> MS, int numSteps) {
        int nMacroStates = MS.size();
        int nStates = Q.getNumRows();

        Triple<Matrix, Double, Double> courtoisResult = Ctmc_courtois.ctmc_courtois(Q, MS);
        Matrix pcourt = courtoisResult.getFirst();
        double eps = courtoisResult.getSecond().doubleValue();
        double epsMAX = courtoisResult.getThird().doubleValue();
        // see _kb/03-api-layer.md for rationale
        Pair<Matrix, Double> randResult = Ctmc_randomization.ctmc_randomization(Q, 1.05 * Q.elementMaxAbs());
        Matrix P = randResult.getFirst();

        Matrix pn = pcourt.copy();

        for (int n = 0; n < numSteps; n++) {
            Matrix pn_1 = pn.copy();

            // see _kb/03-api-layer.md for rationale
            Matrix G = Matrix.zeros(nMacroStates, nMacroStates);
            for (int I = 0; I < nMacroStates; I++) {
                List<Integer> idxI = MS.get(I);
                double S = 0.0;
                for (int i = 0; i < idxI.size(); i++) {
                    S += pn_1.get(idxI.get(i), 0);
                }

                for (int J = 0; J < nMacroStates; J++) {
                    List<Integer> idxJ = MS.get(J);
                    if (I != J) {
                        for (int i = 0; i < idxI.size(); i++) {
                            for (int j = 0; j < idxJ.size(); j++) {
                                if (S > 1e-14) {
                                    G.set(I, J, G.get(I, J) + P.get(idxI.get(i), idxJ.get(j)) * pn_1.get(idxI.get(i), 0) / S);
                                }
                            }
                        }
                    }
                }
            }

            for (int i = 0; i < nMacroStates; i++) {
                double rowSum = 0.0;
                for (int j = 0; j < nMacroStates; j++) {
                    if (i != j) {
                        rowSum += G.get(i, j);
                    }
                }
                G.set(i, i, 1.0 - rowSum);
            }

            Matrix gamma = Dtmc_solve.dtmc_solve(G);

            Matrix GI = Matrix.zeros(nMacroStates, nStates);
            for (int I = 0; I < nMacroStates; I++) {
                List<Integer> idxI = MS.get(I);
                double S = 0.0;
                for (int i = 0; i < idxI.size(); i++) {
                    S += pn_1.get(idxI.get(i), 0);
                }

                for (int j = 0; j < nStates; j++) {
                    double sum = 0.0;
                    for (int i = 0; i < idxI.size(); i++) {
                        sum += P.get(idxI.get(i), j) * pn_1.get(idxI.get(i), 0);
                    }
                    if (S > 0) {
                        GI.set(I, j, sum / S);
                    }
                }
            }

            for (int I = 0; I < nMacroStates; I++) {
                List<Integer> idxI = MS.get(I);
                int blockSize = idxI.size();
                Matrix A = Matrix.eye(blockSize);
                Matrix b = Matrix.zeros(blockSize, 1);

                for (int i = 0; i < blockSize; i++) {
                    for (int j = 0; j < blockSize; j++) {
                        A.set(i, j, A.get(i, j) - P.get(idxI.get(j), idxI.get(i)));
                    }

                    for (int K = 0; K < nMacroStates; K++) {
                        if (K != I) {
                            b.set(i, 0, b.get(i, 0) + gamma.get(K) * GI.get(K, idxI.get(i)));
                        }
                    }
                }

                Matrix x = Matrix.zeros(blockSize, 1);
                // Above the dispatch threshold the direct factorization is what limits
                // the model rather than the arithmetic; it remains the fallback.
                boolean solved = false;
                if (blockSize > Ctmc_solve.GMRES_MIN_STATES) {
                    Ctmc_gmres.GmresResult g = Ctmc_gmres.ctmc_gmres(A, b, 0.0, 0, 0, null);
                    if (g.flag == 0) {
                        x = g.x;
                        solved = true;
                    }
                }
                if (!solved) {
                    Matrix.solve(A, b, x);
                }
                for (int i = 0; i < blockSize; i++) {
                    pn.set(idxI.get(i), 0, x.get(i, 0));
                }
            }

            double sum = 0.0;
            for (int i = 0; i < nStates; i++) {
                sum += pn.get(i, 0);
            }
            for (int i = 0; i < nStates; i++) {
                pn.set(i, 0, pn.get(i, 0) / sum);
            }
        }

        return new Triple<Matrix, Double, Double>(pn, Double.valueOf(eps), Double.valueOf(epsMAX));
    }
}
