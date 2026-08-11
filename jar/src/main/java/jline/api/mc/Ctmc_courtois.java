/**
 * @file Courtois decomposition for nearly completely decomposable CTMCs
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import java.util.ArrayList;
import java.util.List;

import jline.util.Pair;
import jline.util.Triple;

import jline.util.matrix.Matrix;

public final class Ctmc_courtois {
    private Ctmc_courtois() {}

    /**
     * Courtois decomposition for nearly completely decomposable CTMCs
     *
     * @param Q  infinitesimal generator matrix
     * @param MS macro-states: MS.get(i) is the set of rows of Q in macro-state i
     * @param q  randomization coefficient, or null to derive it as 1.05*max|Qperm|
     *           (the MATLAB ctmc_courtois.m default). Must be a valid uniformization
     *           rate (q >= max|q_ii|) or P=I+Qperm/q is not stochastic.
     */
    public static Triple<Matrix, Double, Double> ctmc_courtois(Matrix Q, List<List<Integer>> MS, Double q) {
        int nMacroStates = MS.size();
        int n = Q.getNumRows();

        List<Integer> v = new ArrayList<Integer>();
        for (List<Integer> macroState : MS) {
            v.addAll(macroState);
        }

        Matrix Qperm = Matrix.zeros(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Qperm.set(i, j, Q.get(v.get(i), v.get(j)));
            }
        }

        Matrix Qdec = Qperm.copy();
        int procRows = 0;
        for (int i = 0; i < nMacroStates; i++) {
            int blockSize = MS.get(i).size();
            if (procRows > 0) {
                for (int row = procRows; row < procRows + blockSize; row++) {
                    for (int col = 0; col < procRows; col++) {
                        Qdec.set(row, col, 0.0);
                    }
                }
            }
            for (int row = procRows; row < procRows + blockSize; row++) {
                for (int col = procRows + blockSize; col < n; col++) {
                    Qdec.set(row, col, 0.0);
                }
            }
            procRows += blockSize;
        }

        Matrix QdecFixed = Ctmc_makeinfgen.ctmc_makeinfgen(Qdec);

        // A null q selects the derived default. The former "q == 1.0" sentinel could
        // not express a genuine caller-supplied q of 1.0.
        double qVal;
        if (q == null) {
            qVal = 1.05 * Qperm.elementMaxAbs();
        } else {
            qVal = q.doubleValue();
        }

        Pair<Matrix, Double> randResult = Ctmc_randomization.ctmc_randomization(Qperm, qVal);
        Matrix P = randResult.getFirst();

        Matrix A = P.copy();
        procRows = 0;
        for (int i = 0; i < nMacroStates; i++) {
            int blockSize = MS.get(i).size();
            if (procRows > 0) {
                for (int row = procRows; row < procRows + blockSize; row++) {
                    for (int col = 0; col < procRows; col++) {
                        A.set(row, col, 0.0);
                    }
                }
            }
            for (int row = procRows; row < procRows + blockSize; row++) {
                for (int col = procRows + blockSize; col < n; col++) {
                    A.set(row, col, 0.0);
                }
            }
            procRows += blockSize;
        }

        Matrix B = P.sub(A);
        double eps = B.sumRows().elementMax();

        double epsMAX = computeEpsMax(A, MS);

        // Compute microprobabilities
        Matrix pmicro = Matrix.zeros(n, 1);
        procRows = 0;
        for (int i = 0; i < nMacroStates; i++) {
            int blockSize = MS.get(i).size();
            Matrix Qmicrostate = Matrix.zeros(blockSize, blockSize);
            for (int row = 0; row < blockSize; row++) {
                for (int col = 0; col < blockSize; col++) {
                    Qmicrostate.set(row, col, Qdec.get(procRows + row, procRows + col));
                }
            }
            Matrix microProb = Ctmc_solve.ctmc_solve(Qmicrostate);
            // ctmc_solve returns a 1xn row vector; read across columns
            for (int j = 0; j < blockSize; j++) {
                pmicro.set(procRows + j, 0, microProb.get(0, j));
            }
            procRows += blockSize;
        }

        Matrix G = Matrix.zeros(nMacroStates, nMacroStates);
        procRows = 0;
        for (int i = 0; i < nMacroStates; i++) {
            int blockSizeI = MS.get(i).size();
            int procCols = 0;
            for (int j = 0; j < nMacroStates; j++) {
                int blockSizeJ = MS.get(j).size();
                if (i != j) {
                    for (int iState = 0; iState < blockSizeI; iState++) {
                        double sum = 0.0;
                        for (int jState = 0; jState < blockSizeJ; jState++) {
                            sum += P.get(procRows + iState, procCols + jState);
                        }
                        G.set(i, j, G.get(i, j) + pmicro.get(procRows + iState, 0) * sum);
                    }
                }
                procCols += blockSizeJ;
            }
            procRows += blockSizeI;
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

        // dtmc_solve returns a 1xnMacroStates row vector; read across columns
        Matrix pMacro = Dtmc_solve.dtmc_solve(G);

        Matrix p = Matrix.zeros(n, 1);
        procRows = 0;
        for (int i = 0; i < nMacroStates; i++) {
            int blockSize = MS.get(i).size();
            for (int j = 0; j < blockSize; j++) {
                p.set(procRows + j, 0, pMacro.get(0, i) * pmicro.get(procRows + j, 0));
            }
            procRows += blockSize;
        }

        Matrix pout = Matrix.zeros(n, 1);
        for (int i = 0; i < n; i++) {
            pout.set(v.get(i), 0, p.get(i, 0));
        }

        return new Triple<Matrix, Double, Double>(pout, eps, epsMAX);
    }

    public static Triple<Matrix, Double, Double> ctmc_courtois(Matrix Q, List<List<Integer>> MS) {
        return ctmc_courtois(Q, MS, (Double) null);
    }

    /**
     * Maximum acceptable NCD index epsMAX. Port of the {@code nargout>3} branch of MATLAB
     * {@code ctmc_courtois.m}: each diagonal block of the block-diagonalized randomized
     * matrix A is made stochastic (diagonal absorbs the within-block row deficit), and
     * epsMAX is derived from the largest per-block second-largest eigenvalue modulus.
     */
    static double computeEpsMax(Matrix A, List<List<Integer>> MS) {
        int nMacroStates = MS.size();
        double maxSecondEig = 0.0;
        int procRows = 0;

        for (int i = 0; i < nMacroStates; i++) {
            int blockSize = MS.get(i).size();

            // Make each row of the diagonal block stochastic by placing the
            // normalization condition on the diagonal position.
            Matrix block = Matrix.zeros(blockSize, blockSize);
            for (int row = 0; row < blockSize; row++) {
                double offDiagSum = 0.0;
                for (int col = 0; col < blockSize; col++) {
                    if (col != row) {
                        double val = A.get(procRows + row, procRows + col);
                        block.set(row, col, val);
                        offDiagSum += val;
                    }
                }
                block.set(row, row, 1.0 - offDiagSum);
            }

            // Second-largest eigenvalue modulus of the block (0 if scalar block).
            if (blockSize > 1) {
                List<org.apache.commons.math3.complex.Complex> eigs = block.eig();
                double[] mods = new double[eigs.size()];
                for (int e = 0; e < eigs.size(); e++) {
                    mods[e] = eigs.get(e).abs();
                }
                java.util.Arrays.sort(mods); // ascending
                double secondEig = mods[mods.length - 2];
                if (secondEig > maxSecondEig) {
                    maxSecondEig = secondEig;
                }
            }

            procRows += blockSize;
        }
        return (1.0 - maxSecondEig) / 2.0;
    }
}
