/**
 * @file GI/M/1-type Caudal Characteristic
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;

import jline.util.matrix.Matrix;

public final class GIM1_Caudal {
    private GIM1_Caudal() {}

    public static GIM1CaudalResult gim1_caudal(Matrix A) {
        return gim1_caudal(A, false, false);
    }

    public static GIM1CaudalResult gim1_caudal(Matrix A, boolean dual) {
        return gim1_caudal(A, dual, false);
    }

    /**
     * Computes the dominant eigenvalue of the R matrix for GI/M/1-type chains.
     */
    public static GIM1CaudalResult gim1_caudal(Matrix A, boolean dual, boolean computeEigenvector) {
        int m = A.getNumRows();
        int dega = A.getNumCols() / m - 1;

        double etaMin = 0.0;
        double etaMax = 1.0;
        double eta = 0.5;

        while (etaMax - etaMin > 1e-15) {
            Matrix temp = A.extractCols(dega * m, (dega + 1) * m).copy();
            for (int i = dega - 1; i >= 0; i--) {
                temp = temp.scale(eta).add(A.extractCols(i * m, (i + 1) * m));
            }

            double newEta = computeSpectralRadius(temp);

            if (newEta > eta) {
                etaMin = eta;
            } else {
                etaMax = eta;
            }
            eta = (etaMin + etaMax) / 2;
        }

        Matrix v = null;
        if (computeEigenvector) {
            Matrix temp = A.extractCols(dega * m, (dega + 1) * m).copy();
            for (int i = dega - 1; i >= 0; i--) {
                temp = temp.scale(eta).add(A.extractCols(i * m, (i + 1) * m));
            }
            v = computeDominantEigenvector(temp);
        }

        return new GIM1CaudalResult(eta, v);
    }

    private static double computeSpectralRadius(Matrix A) {
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

            double maxEigenvalue = 0.0;
            for (int i = 0; i < m; i++) {
                double real = decomposition.getRealEigenvalues()[i];
                if (real > maxEigenvalue) {
                    maxEigenvalue = real;
                }
            }
            return maxEigenvalue;
        } catch (Exception e) {
            return powerIterationSpectralRadius(A, 100);
        }
    }

    private static double powerIterationSpectralRadius(Matrix A, int maxIter) {
        int m = A.getNumRows();
        Matrix v = Matrix.ones(m, 1).scale(1.0 / m);

        for (int iter = 0; iter < maxIter; iter++) {
            Matrix Av = A.mult(v);
            double norm = Av.infinityNorm();
            if (norm > 0) {
                v = Av.scale(1.0 / norm);
            }
        }

        Matrix Av = A.mult(v);
        return v.transpose().mult(Av).get(0, 0) / v.transpose().mult(v).get(0, 0);
    }

    private static Matrix computeDominantEigenvector(Matrix A) {
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
