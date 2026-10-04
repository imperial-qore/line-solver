/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference.util;

import java.util.ArrayList;
import java.util.List;

import jline.util.matrix.Matrix;

/**
 * Non-negative least squares solver implementing the Lawson-Hanson algorithm.
 *
 * Solves: min ||Ax - b||^2  subject to x &gt;= 0.
 */
public final class NnlsSolver {
    private NnlsSolver() {}

    /**
     * Solve the NNLS problem: min ||Ax - b||^2 subject to x &gt;= 0.
     */
    public static Matrix lsqnonneg(Matrix A, Matrix b) {
        int m = A.getNumRows();
        int n = A.getNumCols();
        Matrix x = new Matrix(n, 1);

        // Passive and active sets
        boolean[] passive = new boolean[n];
        int maxIter = 3 * n;

        // Compute initial gradient w = A'(b - Ax)
        Matrix w = computeGradient(A, b, x, m, n);

        int iter = 0;
        while (iter < maxIter) {
            // Find index of largest w among active (non-passive) variables
            int tMax = -1;
            double wMax = 0.0;
            for (int j = 0; j < n; j++) {
                if (!passive[j] && w.get(j, 0) > wMax) {
                    wMax = w.get(j, 0);
                    tMax = j;
                }
            }

            // If no positive gradient found, we're done
            if (tMax < 0 || wMax <= 1e-12) break;

            // Move variable tMax to passive set
            passive[tMax] = true;

            // Solve unconstrained LS on passive set
            int innerIter = 0;
            while (innerIter < maxIter) {
                List<Integer> passiveIdx = new ArrayList<Integer>();
                for (int idx = 0; idx < n; idx++) {
                    if (passive[idx]) passiveIdx.add(idx);
                }
                if (passiveIdx.isEmpty()) break;

                Matrix Ap = extractColumns(A, passiveIdx, m);
                Matrix sp = solveLeastSquares(Ap, b);

                // Check if all passive variables are positive
                boolean allPositive = true;
                for (int i = 0; i < passiveIdx.size(); i++) {
                    if (sp.get(i, 0) <= 0.0) {
                        allPositive = false;
                        break;
                    }
                }

                if (allPositive) {
                    for (int i = 0; i < passiveIdx.size(); i++) {
                        x.set(passiveIdx.get(i), 0, sp.get(i, 0));
                    }
                    break;
                }

                // Find alpha: min ratio x[j]/(x[j]-s[j]) for s[j] <= 0
                double alpha = Double.MAX_VALUE;
                for (int i = 0; i < passiveIdx.size(); i++) {
                    if (sp.get(i, 0) <= 0.0) {
                        double ratio = x.get(passiveIdx.get(i), 0) /
                                (x.get(passiveIdx.get(i), 0) - sp.get(i, 0));
                        if (ratio < alpha) {
                            alpha = ratio;
                        }
                    }
                }

                // Update x = x + alpha * (s - x)
                for (int i = 0; i < passiveIdx.size(); i++) {
                    int j = passiveIdx.get(i);
                    x.set(j, 0, x.get(j, 0) + alpha * (sp.get(i, 0) - x.get(j, 0)));
                }

                // Move variables with zero x back to active set
                for (int i = 0; i < passiveIdx.size(); i++) {
                    if (Math.abs(x.get(passiveIdx.get(i), 0)) < 1e-12) {
                        passive[passiveIdx.get(i)] = false;
                        x.set(passiveIdx.get(i), 0, 0.0);
                    }
                }

                innerIter++;
            }

            w = computeGradient(A, b, x, m, n);
            iter++;
        }

        return x;
    }

    /**
     * Solve NNLS with double array inputs: min ||Ax - b||^2 subject to x &gt;= 0.
     */
    public static double[] lsqnonneg(double[][] A, double[] b) {
        int m = A.length;
        int n = (m > 0) ? A[0].length : 0;
        Matrix matA = new Matrix(m, n);
        Matrix vecB = new Matrix(m, 1);
        for (int i = 0; i < m; i++) {
            vecB.set(i, 0, b[i]);
            for (int j = 0; j < n; j++) {
                matA.set(i, j, A[i][j]);
            }
        }
        Matrix result = lsqnonneg(matA, vecB);
        double[] out = new double[n];
        for (int i = 0; i < n; i++) out[i] = result.get(i, 0);
        return out;
    }

    private static Matrix computeGradient(Matrix A, Matrix b, Matrix x, int m, int n) {
        // w = A'(b - Ax)
        Matrix residual = new Matrix(m, 1);
        for (int i = 0; i < m; i++) {
            double ax = 0.0;
            for (int j = 0; j < n; j++) {
                ax += A.get(i, j) * x.get(j, 0);
            }
            residual.set(i, 0, b.get(i, 0) - ax);
        }
        Matrix w = new Matrix(n, 1);
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int i = 0; i < m; i++) {
                sum += A.get(i, j) * residual.get(i, 0);
            }
            w.set(j, 0, sum);
        }
        return w;
    }

    private static Matrix extractColumns(Matrix A, List<Integer> cols, int m) {
        Matrix result = new Matrix(m, cols.size());
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < cols.size(); j++) {
                result.set(i, j, A.get(i, cols.get(j)));
            }
        }
        return result;
    }

    private static Matrix solveLeastSquares(Matrix A, Matrix b) {
        // Solve A'Ax = A'b using Cholesky or direct method
        int m = A.getNumRows();
        int n = A.getNumCols();
        Matrix AtA = new Matrix(n, n);
        Matrix Atb = new Matrix(n, 1);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double sum = 0.0;
                for (int k = 0; k < m; k++) {
                    sum += A.get(k, i) * A.get(k, j);
                }
                AtA.set(i, j, sum);
            }
            double sum = 0.0;
            for (int k = 0; k < m; k++) {
                sum += A.get(k, i) * b.get(k, 0);
            }
            Atb.set(i, 0, sum);
        }

        // Solve using Gaussian elimination with partial pivoting
        return solveLinear(AtA, Atb);
    }

    private static Matrix solveLinear(Matrix A, Matrix b) {
        int n = A.getNumRows();
        double[][] aug = new double[n][n + 1];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n + 1; j++) {
                aug[i][j] = (j < n) ? A.get(i, j) : b.get(i, 0);
            }
        }

        // Forward elimination with partial pivoting
        for (int col = 0; col < n; col++) {
            int maxRow = col;
            double maxVal = Math.abs(aug[col][col]);
            for (int row = col + 1; row < n; row++) {
                if (Math.abs(aug[row][col]) > maxVal) {
                    maxVal = Math.abs(aug[row][col]);
                    maxRow = row;
                }
            }
            double[] temp = aug[col];
            aug[col] = aug[maxRow];
            aug[maxRow] = temp;

            double pivot = aug[col][col];
            if (Math.abs(pivot) < 1e-14) continue;

            for (int row = col + 1; row < n; row++) {
                double factor = aug[row][col] / pivot;
                for (int j = col; j < n + 1; j++) {
                    aug[row][j] -= factor * aug[col][j];
                }
            }
        }

        // Back substitution
        Matrix x = new Matrix(n, 1);
        for (int i = n - 1; i >= 0; i--) {
            double sum = aug[i][n];
            for (int j = i + 1; j < n; j++) {
                sum -= aug[i][j] * x.get(j, 0);
            }
            if (Math.abs(aug[i][i]) > 1e-14) {
                x.set(i, 0, sum / aug[i][i]);
            }
        }
        return x;
    }
}
