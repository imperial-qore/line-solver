/**
 * @file Compute T-matrix using Sylvester equation approach
 *
 * @since LINE 3.0
 */
package jline.lib.fjcodes;

import java.lang.reflect.Constructor;
import java.lang.reflect.Method;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import jline.util.matrix.Matrix;

public final class ComputeT_Sylvester {
    private ComputeT_Sylvester() {}

    /**
     * Compute T-matrix via Sylvester equation method
     */
    public static Matrix computeT_Sylvester(Matrix D0, Matrix D1, Matrix S, Matrix A_jump) {
        int d0 = D0.getNumRows();
        int ms = S.getNumRows() / d0;
        int m = ms * d0;

        Matrix Ims = Matrix.eye(ms);
        Matrix ID0 = Ims.kron(D0);

        Matrix A_jump_Arr = A_jump.kron(Matrix.eye(d0));
        Matrix DS = Ims.kron(D1).mult(A_jump_Arr);

        Array2DRowRealMatrix ID0_commons = new Array2DRowRealMatrix(ID0.toArray2D());

        EigenDecomposition eigen = new EigenDecomposition(ID0_commons);
        RealMatrix U_commons = eigen.getV();
        RealMatrix D_commons = eigen.getD();

        Matrix U = new Matrix(U_commons.getData());
        Matrix Tr = new Matrix(D_commons.getData());

        Matrix Tnew = S.copy();
        Matrix Told = Matrix.zeros(m, m);
        double tolerance = 1e-12;
        int maxIter = 200;
        int iter = 0;

        while (iter < maxIter) {
            double maxDiff = 0.0;
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    maxDiff = Math.max(maxDiff, Math.abs(Told.get(i, j) - Tnew.get(i, j)));
                }
            }

            if (maxDiff <= tolerance) {
                break;
            }

            Told = Tnew.copy();

            Matrix L = qSylvest(U, Tr, Tnew);

            Tnew = S.add(1.0, L.mult(DS));

            iter++;
        }

        Matrix L_final = qSylvest(U, Tr, Tnew);
        Matrix termTL = Tnew.mult(L_final);
        Matrix termLID0 = L_final.mult(ID0);
        Matrix termI = Matrix.eye(m);

        Matrix residual = termTL.add(1.0, termLID0).add(1.0, termI);
        double residualNorm = normInf(residual);
        if (residualNorm > 1e-6) {
            // suppressed warning per source
        }

        return Tnew;
    }

    /**
     * Solve Sylvester equation X*kron(A,I) + BX = -I using Hessenberg decomposition
     */
    private static Matrix qSylvest(Matrix U, Matrix T, Matrix B) {
        int n = U.getNumCols();

        Array2DRowRealMatrix B_commons = new Array2DRowRealMatrix(B.toArray2D());

        try {
            Class<?> hessClass = Class.forName("org.apache.commons.math3.linear.HessenbergTransformer");
            Constructor<?> constructor = hessClass.getDeclaredConstructor(RealMatrix.class);
            constructor.setAccessible(true);
            Object hessTransformer = constructor.newInstance(B_commons);

            Method getH = hessClass.getDeclaredMethod("getH");
            getH.setAccessible(true);
            RealMatrix NBAR_commons = (RealMatrix) getH.invoke(hessTransformer);

            Method getP = hessClass.getDeclaredMethod("getP");
            getP.setAccessible(true);
            RealMatrix V_commons = (RealMatrix) getP.invoke(hessTransformer);

            Matrix V = new Matrix(V_commons.getData());
            Matrix NBAR = new Matrix(NBAR_commons.getData());

            Matrix F = V.transpose().mult(U).scale(-1.0);

            Matrix Y = Matrix.zeros(n, n);

            for (int k = 0; k < n; k++) {
                Matrix temp = Matrix.zeros(n, 1);

                if (k == 0) {
                    for (int i = 0; i < n; i++) {
                        temp.set(i, 0, F.get(i, k));
                    }
                } else {
                    for (int i = 0; i < n; i++) {
                        double sum = F.get(i, k);
                        for (int j = 0; j < k; j++) {
                            sum -= Y.get(i, j) * T.get(j, k);
                        }
                        temp.set(i, 0, sum);
                    }
                }

                Matrix A_sys = NBAR.add(T.get(k, k), Matrix.eye(n));

                Array2DRowRealMatrix A_sys_commons = new Array2DRowRealMatrix(A_sys.toArray2D());
                Array2DRowRealMatrix temp_commons = new Array2DRowRealMatrix(temp.toArray2D());

                DecompositionSolver solver = new LUDecomposition(A_sys_commons).getSolver();
                RealMatrix y_k_commons = solver.solve(temp_commons);

                for (int i = 0; i < n; i++) {
                    Y.set(i, k, y_k_commons.getEntry(i, 0));
                }
            }

            Matrix X = V.mult(Y).mult(U.transpose());

            return X;
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Compute infinity norm (max row sum of absolute values)
     */
    private static double normInf(Matrix A) {
        double maxSum = 0.0;
        for (int i = 0; i < A.getNumRows(); i++) {
            double rowSum = 0.0;
            for (int j = 0; j < A.getNumCols(); j++) {
                rowSum += Math.abs(A.get(i, j));
            }
            maxSum = Math.max(maxSum, rowSum);
        }
        return maxSum;
    }
}
