/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.perm;

import jline.util.matrix.Matrix;

/**
 * Implementation of Sum Product Algorithm (SPA) to approximate the Bethe permanent.
 */
public class BethePermanent extends PermSolver {

    private final double epsilon;
    private final int maxIteration;
    private final Matrix matrixSqrt;
    private static final double MIN_VALUE = 2.220446049250314e-16;

    public BethePermanent(Matrix matrix) {
        this(matrix, 0.001, 200000, false);
    }

    public BethePermanent(Matrix matrix, double epsilon, int maxIteration) {
        this(matrix, epsilon, maxIteration, false);
    }

    public BethePermanent(Matrix matrix, double epsilon, int maxIteration, boolean solve) {
        super(matrix);
        this.epsilon = epsilon;
        this.maxIteration = maxIteration;

        // Ensure all matrix elements are at least minValue to avoid numerical issues
        double[][] data = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                data[i][j] = Math.max(matrix.get(i, j), MIN_VALUE);
            }
        }
        Matrix sqrtMat = new Matrix(data);
        for (int i = 0; i < sqrtMat.getNumRows(); i++) {
            for (int j = 0; j < sqrtMat.getNumCols(); j++) {
                sqrtMat.set(i, j, Math.sqrt(sqrtMat.get(i, j)));
            }
        }
        this.matrixSqrt = sqrtMat;

        if (solve) {
            solve();
        }
    }

    @Override
    public void compute() {
        value = spa();
    }

    /**
     * Sum Product Algorithm implementation for Bethe permanent approximation.
     */
    private double spa() {
        Matrix rPast = Matrix.ones(n, n);
        Matrix lPast = Matrix.ones(n, n);

        Matrix[] rl = update(lPast);
        Matrix r = rl[0];
        Matrix l = rl[1];

        int iteration = 0;
        while (convergence(rPast, lPast, r, l) > epsilon && iteration < maxIteration) {
            iteration++;
            rPast = r.copy();
            lPast = l.copy();

            Matrix[] newRL = update(lPast);
            r = newRL[0];
            l = newRL[1];
        }

        return bethe(l, r);
    }

    /**
     * Update the right going message and the left going message.
     */
    private Matrix[] update(Matrix l) {
        Matrix r1 = Matrix.zeros(n, n);
        Matrix l1 = Matrix.zeros(n, n);

        for (int i = 0; i < n; i++) {
            double denomSum = 0.0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    denomSum += matrixSqrt.get(i, j) * l.get(i, j);
                }
            }
            for (int j = 0; j < n; j++) {
                r1.set(i, j, matrixSqrt.get(i, j) / denomSum);
            }
        }

        for (int j = 0; j < n; j++) {
            double denomSum = 0.0;
            for (int i = 0; i < n; i++) {
                if (i != j) {
                    denomSum += matrixSqrt.get(i, j) * r1.get(i, j);
                }
            }
            for (int i = 0; i < n; i++) {
                l1.set(i, j, matrixSqrt.get(i, j) / denomSum);
            }
        }

        return new Matrix[] { r1, l1 };
    }

    /**
     * Calculate the squared difference between past and present messages to measure convergence.
     */
    private double convergence(Matrix r0, Matrix l0, Matrix r1, Matrix l1) {
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double rDiff = r0.get(i, j) - r1.get(i, j);
                double lDiff = l0.get(i, j) - l1.get(i, j);
                sum += rDiff * rDiff + lDiff * lDiff;
            }
        }
        return sum;
    }

    /**
     * Compute the Bethe permanent from the left and right going messages.
     */
    private double bethe(Matrix l, Matrix r) {
        double[] term1 = new double[n];
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                sum += matrixSqrt.get(i, j) * l.get(i, j);
            }
            term1[i] = Math.max(sum, MIN_VALUE);
        }

        double[] term2 = new double[n];
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int i = 0; i < n; i++) {
                sum += matrixSqrt.get(i, j) * r.get(i, j);
            }
            term2[j] = Math.max(sum, MIN_VALUE);
        }

        double[][] term3 = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                term3[i][j] = Math.max(r.get(i, j) * l.get(i, j) + 1.0, MIN_VALUE);
            }
        }

        double lterm1 = 0.0;
        for (int i = 0; i < n; i++) {
            lterm1 -= Math.log(term1[i]);
        }

        double lterm2 = 0.0;
        for (int j = 0; j < n; j++) {
            lterm2 -= Math.log(term2[j]);
        }

        double lterm3 = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                lterm3 += Math.log(term3[i][j]);
            }
        }

        double result = Math.exp(-(lterm1 + lterm2 + lterm3));
        if (Double.isNaN(result) || Double.isInfinite(result)) {
            return 0.0;
        }
        return result;
    }
}
