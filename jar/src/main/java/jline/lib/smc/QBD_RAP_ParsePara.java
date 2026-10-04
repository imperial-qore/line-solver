package jline.lib.smc;

import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.interfaces.decomposition.EigenDecomposition_F64;

import jline.util.matrix.Matrix;

public final class QBD_RAP_ParsePara {
    private QBD_RAP_ParsePara() {}

    /**
     * QBD_RAP_ParsePara checks the validity of the input matrices A0, A1 and A2
     * for a QBD with RAP components. The evaluated conditions are necessary but not
     * sufficient.
     */
    public static void QBD_RAP_ParsePara(Matrix A0, Matrix A1, Matrix A2) {
        // Check dimensions
        if (A0.getNumRows() != A0.getNumCols()) {
            throw new IllegalArgumentException("A0 is not a square matrix");
        }
        if (A1.getNumRows() != A1.getNumCols()) {
            throw new IllegalArgumentException("A1 is not a square matrix");
        }
        if (A2.getNumRows() != A2.getNumCols()) {
            throw new IllegalArgumentException("A2 is not a square matrix");
        }
        if (A0.getNumRows() != A1.getNumRows()) {
            throw new IllegalArgumentException("The matrices A0 and A1 do not have the same dimension");
        }
        if (A0.getNumRows() != A2.getNumRows()) {
            throw new IllegalArgumentException("The matrices A0 and A2 do not have the same dimension");
        }

        // Check zero row sum
        Matrix sum_matrix = A0.add(1.0, A1).add(1.0, A2);
        Matrix row_sums = sum_matrix.sumRows();
        double max_row_sum = row_sums.elementMax();
        double min_row_sum = row_sums.elementMin();

        if (max_row_sum > 1e-14 || min_row_sum < -1e-14) {
            throw new IllegalArgumentException("The matrix A0+A1+A2 must have zero row sum");
        }

        // Check dominant eigenvalue of A1: it must have negative real part,
        // otherwise the local block is not a proper RAP sub-generator and the
        // cyclic reduction below does not converge. Mirrors QBD_RAP_ParsePara.m.
        double maxRealA1 = maxRealEigenvalue(A1);
        if (maxRealA1 > -1e-14) {
            throw new IllegalArgumentException(
                    "The dominant eigenvalue of the matrix A1 must have negative real part");
        }

        // A0+A1+A2 is the generator of the phase process, so its dominant
        // eigenvalue must be zero. Mirrors QBD_RAP_ParsePara.m.
        double maxRealSum = maxRealEigenvalue(sum_matrix);
        if (maxRealSum > 1e-14 || maxRealSum < -1e-14) {
            throw new IllegalArgumentException(
                    "The dominant eigenvalue of the matrix A0+A1+A2 must have zero real part");
        }
    }

    /**
     * Largest real part over the eigenvalues of a square matrix.
     */
    private static double maxRealEigenvalue(Matrix A) {
        int n = A.getNumRows();
        DMatrixRMaj dm = new DMatrixRMaj(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                dm.set(i, j, A.get(i, j));
            }
        }
        EigenDecomposition_F64<DMatrixRMaj> evd = DecompositionFactory_DDRM.eig(n, false);
        if (!evd.decompose(dm)) {
            throw new IllegalArgumentException(
                    "Eigendecomposition failed while validating the QBD blocks with RAP components");
        }
        double maxReal = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < evd.getNumberOfEigenvalues(); i++) {
            double re = evd.getEigenvalue(i).real;
            if (re > maxReal) {
                maxReal = re;
            }
        }
        return maxReal;
    }
}
