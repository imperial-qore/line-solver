/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 *
 * Reference:
 * P. Buchholz, M. Telek, "On minimal representation of rational arrival
 * processes." Madrid Conference on Queueing theory (MCQT), June 2010.
 */
package jline.lib.butools.reptrans;

import java.util.ArrayList;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.SingularValueDecomposition;

import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class MStaircase {
    private MStaircase() {}

    public static Pair<Matrix, Integer> mStaircase(java.util.List<Matrix> Y, Matrix Z) {
        return mStaircase(Y, Z, 1e-12);
    }

    /**
     * Computes a smaller representation using the staircase algorithm.
     */
    public static Pair<Matrix, Integer> mStaircase(java.util.List<Matrix> Y, Matrix Z, double precision) {
        int MCount = Y.size();
        int m = Y.get(0).getNumRows();

        ArrayList<Matrix> X = new ArrayList<Matrix>(MCount);
        for (Matrix y : Y) {
            X.add(y.copy());
        }
        Matrix Zwork = Z.copy();

        Matrix U = Matrix.eye(m);
        int ranksum = 0;
        boolean crit = true;

        while (crit) {
            int r = matrixRank(Zwork, precision);
            ranksum += r;

            Matrix Ui = fullSvdU(Zwork);
            int mCurrent = Ui.getNumRows();

            Matrix Transf = Matrix.eye(m);
            int offset = m - mCurrent;
            for (int i = 0; i < mCurrent; i++) {
                for (int j = 0; j < mCurrent; j++) {
                    Transf.set(offset + i, offset + j, Ui.get(j, i));
                }
            }

            Matrix Uold = U.copy();
            Matrix newU = new Matrix(m, m);
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    double sum = 0.0;
                    for (int k = 0; k < m; k++) {
                        sum += Uold.get(i, k) * Transf.get(j, k);
                    }
                    newU.set(i, j, sum);
                }
            }
            U = newU;

            Matrix newZ = null;
            for (int i = 0; i < MCount; i++) {
                Matrix TEMP = Ui.transpose().mult(X.get(i)).mult(Ui);
                int currentSize = TEMP.getNumRows();

                int newSize = currentSize - r;
                Matrix newXi = new Matrix(newSize, newSize);
                for (int ri = 0; ri < newSize; ri++) {
                    for (int ci = 0; ci < newSize; ci++) {
                        newXi.set(ri, ci, TEMP.get(r + ri, r + ci));
                    }
                }
                X.set(i, newXi);

                Matrix block = new Matrix(newSize, r);
                for (int ri = 0; ri < newSize; ri++) {
                    for (int ci = 0; ci < r; ci++) {
                        block.set(ri, ci, TEMP.get(r + ri, ci));
                    }
                }

                if (i == 0) {
                    newZ = block;
                } else {
                    int oldCols = newZ.getNumCols();
                    Matrix combined = new Matrix(newSize, oldCols + r);
                    for (int ri = 0; ri < newSize; ri++) {
                        for (int ci = 0; ci < oldCols; ci++) {
                            combined.set(ri, ci, newZ.get(ri, ci));
                        }
                        for (int ci = 0; ci < r; ci++) {
                            combined.set(ri, oldCols + ci, block.get(ri, ci));
                        }
                    }
                    newZ = combined;
                }
            }
            Zwork = newZ;

            if (Zwork.norm() < precision || matrixRank(Zwork, precision) == m - ranksum) {
                crit = false;
            }
        }

        int n = ranksum;

        if (Zwork.norm() < precision) {
            double[] xFull = new double[m];
            for (int j = 0; j < m; j++) {
                double colSum = 0.0;
                for (int i = 0; i < m; i++) {
                    colSum += U.get(i, j);
                }
                xFull[j] = colSum;
            }

            boolean hasZero = false;
            boolean[] zeroloc = new boolean[n];
            int nonzeroIdx = -1;
            for (int l = 0; l < n; l++) {
                if (Math.abs(xFull[l]) < precision) {
                    hasZero = true;
                    zeroloc[l] = true;
                } else if (nonzeroIdx == -1) {
                    nonzeroIdx = l;
                }
            }

            Matrix R = Matrix.eye(n);
            if (hasZero && nonzeroIdx >= 0) {
                for (int l = 0; l < n; l++) {
                    if (zeroloc[l]) {
                        R.set(l, nonzeroIdx, 1.0);
                    }
                }
            }

            Matrix xVec = new Matrix(n, 1);
            for (int i = 0; i < n; i++) {
                xVec.set(i, 0, xFull[i]);
            }
            Matrix yVec = R.mult(xVec);

            double[] gammaVals = new double[n];
            for (int i = 0; i < n; i++) {
                gammaVals[i] = yVec.get(i, 0);
            }
            Matrix Gamma = Matrix.diag(gammaVals);

            Matrix TEMP1 = Matrix.eye(m);
            Matrix GammaInv = Gamma.inv();
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    TEMP1.set(i, j, GammaInv.get(i, j));
                }
            }

            Matrix TEMP2 = Matrix.eye(m);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    TEMP2.set(i, j, R.get(i, j));
                }
            }

            Matrix UT = U.transpose();
            Matrix B = TEMP1.mult(TEMP2).mult(UT).inv();

            return new Pair<Matrix, Integer>(B, Integer.valueOf(n));
        } else {
            return new Pair<Matrix, Integer>(Matrix.eye(m), Integer.valueOf(m));
        }
    }

    private static int matrixRank(Matrix A, double tol) {
        int rows = A.getNumRows();
        int cols = A.getNumCols();
        if (rows == 0 || cols == 0) return 0;

        org.apache.commons.math3.linear.RealMatrix realMatrix =
                MatrixUtils.createRealMatrix(A.toArray2D());
        SingularValueDecomposition svd = new SingularValueDecomposition(realMatrix);
        double[] singularValues = svd.getSingularValues();
        int rank = 0;
        for (double sv : singularValues) {
            if (sv > tol) {
                rank++;
            }
        }
        return rank;
    }

    private static Matrix fullSvdU(Matrix A) {
        int rows = A.getNumRows();
        int cols = A.getNumCols();
        int p = Math.min(rows, cols);

        org.apache.commons.math3.linear.RealMatrix realMatrix =
                MatrixUtils.createRealMatrix(A.toArray2D());
        SingularValueDecomposition svd = new SingularValueDecomposition(realMatrix);
        org.apache.commons.math3.linear.RealMatrix uCompact = svd.getU();

        if (p >= rows) {
            Matrix result = new Matrix(rows, rows);
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < rows; j++) {
                    result.set(i, j, uCompact.getEntry(i, j));
                }
            }
            return result;
        }

        Matrix result = new Matrix(rows, rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < p; j++) {
                result.set(i, j, uCompact.getEntry(i, j));
            }
        }

        int colIdx = p;
        for (int candidate = 0; candidate < rows; candidate++) {
            if (colIdx >= rows) break;

            double[] v = new double[rows];
            v[candidate] = 1.0;

            for (int j = 0; j < colIdx; j++) {
                double dot = 0.0;
                for (int i = 0; i < rows; i++) {
                    dot += v[i] * result.get(i, j);
                }
                for (int i = 0; i < rows; i++) {
                    v[i] -= dot * result.get(i, j);
                }
            }

            double norm = 0.0;
            for (int i = 0; i < rows; i++) {
                norm += v[i] * v[i];
            }
            norm = Math.sqrt(norm);

            if (norm > 1e-14) {
                for (int i = 0; i < rows; i++) {
                    result.set(i, colIdx, v[i] / norm);
                }
                colIdx++;
            }
        }

        return result;
    }
}
