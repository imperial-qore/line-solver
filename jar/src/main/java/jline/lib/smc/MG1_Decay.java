/**
 * @file M/G/1-type Decay Rate
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;

import jline.util.matrix.Matrix;

public final class MG1_Decay {
    private MG1_Decay() {}

    public static MG1DecayResult mg1_decay(Matrix A) {
        return mg1_decay(A, true);
    }

    /**
     * Computes the decay rate of a recurrent M/G/1 type Markov chain.
     */
    public static MG1DecayResult mg1_decay(Matrix A, boolean computeEigenvector) {
        int m = A.getNumRows();
        int dega = A.getNumCols() / m - 1;

        // Integer expansion phase
        double eta = 1.0;
        double newEta = 0.0;
        while (newEta - eta < 0) {
            eta += 1.0;
            Matrix temp = A.extractCols(dega * m, (dega + 1) * m).copy();
            for (int i = dega - 1; i >= 0; i--) {
                temp = temp.scale(eta).add(A.extractCols(i * m, (i + 1) * m));
            }
            newEta = computeSpectralRadiusMG1(temp);
        }

        double etaMin = eta - 1.0;
        double etaMax = eta;
        eta = etaMin + 0.5;

        while (etaMax - etaMin > 1e-15) {
            Matrix temp = A.extractCols(dega * m, (dega + 1) * m).copy();
            for (int i = dega - 1; i >= 0; i--) {
                temp = temp.scale(eta).add(A.extractCols(i * m, (i + 1) * m));
            }
            newEta = computeSpectralRadiusMG1(temp);
            if (newEta < eta) {
                etaMin = eta;
            } else {
                etaMax = eta;
            }
            eta = (etaMin + etaMax) / 2.0;
        }

        Matrix uT = null;
        if (computeEigenvector) {
            Matrix temp = A.extractCols(dega * m, (dega + 1) * m).copy();
            for (int i = dega - 1; i >= 0; i--) {
                temp = temp.scale(eta).add(A.extractCols(i * m, (i + 1) * m));
            }
            Matrix v = computeDominantEigenvectorMG1(temp.transpose());
            uT = v.transpose();
        }

        return new MG1DecayResult(eta, uT);
    }

    private static double computeSpectralRadiusMG1(Matrix A) {
        try {
            int m = A.getNumRows();
            double[][] array = new double[m][m];
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    array[i][j] = A.get(i, j);
                }
            }
            Array2DRowRealMatrix realMatrix = new Array2DRowRealMatrix(array);
            EigenDecomposition decomposition = new EigenDecomposition(realMatrix);

            double maxReal = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < m; i++) {
                double real = decomposition.getRealEigenvalues()[i];
                if (real > maxReal) {
                    maxReal = real;
                }
            }
            return maxReal;
        } catch (Exception e) {
            int m = A.getNumRows();
            Matrix v = Matrix.ones(m, 1).scale(1.0 / m);
            for (int iter = 0; iter < 100; iter++) {
                Matrix Av = A.mult(v);
                double norm = Av.infinityNorm();
                if (norm > 0) {
                    v = Av.scale(1.0 / norm);
                }
            }
            Matrix Av = A.mult(v);
            return v.transpose().mult(Av).get(0, 0) / v.transpose().mult(v).get(0, 0);
        }
    }

    private static Matrix computeDominantEigenvectorMG1(Matrix A) {
        try {
            int m = A.getNumRows();
            double[][] array = new double[m][m];
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    array[i][j] = A.get(i, j);
                }
            }
            Array2DRowRealMatrix realMatrix = new Array2DRowRealMatrix(array);
            EigenDecomposition decomposition = new EigenDecomposition(realMatrix);

            int maxIdx = 0;
            double maxReal = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < m; i++) {
                double real = decomposition.getRealEigenvalues()[i];
                if (real > maxReal) {
                    maxReal = real;
                    maxIdx = i;
                }
            }

            org.apache.commons.math3.linear.RealVector eigenvector = decomposition.getEigenvector(maxIdx);
            Matrix result = new Matrix(m, 1);
            for (int i = 0; i < m; i++) {
                result.set(i, 0, eigenvector.getEntry(i));
            }
            return result;
        } catch (Exception e) {
            int m = A.getNumRows();
            Matrix v = Matrix.ones(m, 1).scale(1.0 / m);
            for (int iter = 0; iter < 100; iter++) {
                Matrix Av = A.mult(v);
                double norm = Av.infinityNorm();
                if (norm > 0) {
                    v = Av.scale(1.0 / norm);
                }
            }
            return v;
        }
    }
}
