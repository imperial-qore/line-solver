/**
 * @file M/G/1-type Explicit G computation
 *
 * Determines G directly when rank(A0)=1 for M/G/1-type Markov chains.
 * This is an optimization that avoids iterative computation when possible.
 *
 * Based on the SMC Solver implementation by Benny Van Houdt.
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc;

import jline.util.matrix.Matrix;

public final class MG1_EG {
    private MG1_EG() {}

    /**
     * Determines G directly if rank(A0)=1 for M/G/1-type Markov chains.
     *
     * @param A The block matrix [A0 A1 A2 ... A_max] with m rows and m*(max+1) columns
     * @param verbose When true, prints residual error
     * @return G matrix if explicit solution exists, null otherwise
     */
    public static Matrix mg1_eg(Matrix A, boolean verbose) {
        int m = A.getNumRows();
        int dega = A.getNumCols() / m - 1;

        Matrix sumA = A.extractCols(dega * m, (dega + 1) * m);
        Matrix beta = sumA.sumRows();

        for (int i = dega - 1; i >= 1; i--) {
            sumA = sumA.add(A.extractCols(i * m, (i + 1) * m));
            beta = beta.add(sumA.sumRows());
        }
        sumA = sumA.add(A.extractCols(0, m));

        Matrix theta = Stat.stat(sumA);
        double drift = theta.mult(beta).get(0, 0);

        Matrix G = null;

        if (drift < 1) {
            Matrix A0 = A.extractCols(0, m);
            if (matrixRank(A0) == 1) {
                Matrix rowSums = A0.sumRows();
                int temp = -1;
                for (int i = 0; i < m; i++) {
                    if (rowSums.get(i, 0) > 0) {
                        temp = i;
                        break;
                    }
                }
                if (temp >= 0) {
                    double rowSum = A0.extractRows(temp, temp + 1).elementSum();
                    Matrix betaVec = A0.extractRows(temp, temp + 1).scale(1.0 / rowSum);
                    G = Matrix.ones(m, 1).mult(betaVec);
                }
            }
        } else if (drift > 1) {
            Matrix A0 = A.extractCols(0, m);
            if (matrixRank(A0) == 1) {
                Matrix Atransformed = new Matrix(m, A.getNumCols());
                for (int i = 0; i <= dega; i++) {
                    Matrix Ai = A.extractCols(i * m, (i + 1) * m);
                    Matrix thetaInvDiag = Matrix.diag(theta.scale(-1.0).exp().getRow(0).toArray1D());
                    Matrix thetaDiag = Matrix.diag(theta.getRow(0).toArray1D());
                    Matrix transformedAi = thetaInvDiag.mult(Ai.transpose()).mult(thetaDiag);
                    for (int r = 0; r < m; r++) {
                        for (int c = 0; c < m; c++) {
                            Atransformed.set(r, i * m + c, transformedAi.get(r, c));
                        }
                    }
                }

                Double etahat = gim1_caudal_scalar(Atransformed);

                if (etahat != null) {
                    Matrix temp = Atransformed.extractCols(dega * m, (dega + 1) * m);
                    for (int i = dega - 1; i >= 1; i--) {
                        temp = temp.scale(etahat.doubleValue()).add(Atransformed.extractCols(i * m, (i + 1) * m));
                    }

                    Matrix A0t = Atransformed.extractCols(0, m);
                    Matrix ImT = Matrix.eye(m).sub(temp);
                    Matrix invImT = ImT.inv();
                    Matrix product = A0t.mult(invImT).transpose();

                    Matrix thetaInvDiag2 = Matrix.diag(theta.scale(-1.0).exp().getRow(0).toArray1D());
                    Matrix thetaDiag2 = Matrix.diag(theta.getRow(0).toArray1D());
                    G = thetaInvDiag2.mult(product).mult(thetaDiag2);
                }
            }
        }

        return G;
    }

    public static Matrix mg1_eg(Matrix A) {
        return mg1_eg(A, false);
    }

    /**
     * Compute matrix rank (simplified estimation)
     */
    private static int matrixRank(Matrix A) {
        jline.io.Ret.SVD svd = A.svd();
        Matrix singularValues = svd.s;
        double tol = 1e-10 * singularValues.get(0, 0);
        int rank = 0;
        int min = Math.min(A.getNumRows(), A.getNumCols());
        for (int i = 0; i < min; i++) {
            if (singularValues.get(i, 0) > tol) {
                rank++;
            }
        }
        return rank;
    }

    /**
     * Wrapper for GIM1_Caudal function. Returns the scalar eta value, or null on failure.
     */
    private static Double gim1_caudal_scalar(Matrix A) {
        try {
            GIM1CaudalResult result = GIM1_Caudal.gim1_caudal(A, false, false);
            double eta = result.getEta();
            if (Double.isNaN(eta) || eta <= 0.0 || eta >= 1.0) {
                return null;
            }
            return Double.valueOf(eta);
        } catch (Exception e) {
            return null;
        }
    }
}
