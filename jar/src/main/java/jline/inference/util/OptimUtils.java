/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference.util;

import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.analysis.MultivariateFunction;

import jline.util.Pair;
import jline.util.matrix.Matrix;

/**
 * Optimization utilities for inference algorithms.
 */
public final class OptimUtils {
    private OptimUtils() {}

    /**
     * Bound-constrained minimization using BOBYQA (derivative-free).
     */
    public static Pair<double[], Double> fmincon(final java.util.function.Function<double[], Double> objFun,
                                                 double[] x0, double[] lb, double[] ub, int maxIter) {
        int n = x0.length;
        int interpPoints = 2 * n + 1;

        // Ensure x0 is within bounds
        final double[] x0Safe = new double[n];
        for (int i = 0; i < n; i++) {
            x0Safe[i] = Math.max(lb[i] + 1e-10, Math.min(ub[i] - 1e-10, x0[i]));
        }

        // Compute initial trust region radius
        double initialRadius = 0.0;
        for (int i = 0; i < n; i++) {
            double range = ub[i] - lb[i];
            if (Double.isFinite(range)) {
                initialRadius = Math.max(initialRadius, range * 0.1);
            } else {
                initialRadius = Math.max(initialRadius, Math.abs(x0Safe[i]) * 0.1 + 1.0);
            }
        }
        if (initialRadius < 1e-8) initialRadius = 1.0;

        BOBYQAOptimizer optimizer = new BOBYQAOptimizer(interpPoints, initialRadius, 1e-8);

        try {
            org.apache.commons.math3.optim.PointValuePair result = optimizer.optimize(
                    new MaxEval(maxIter * 10),
                    new ObjectiveFunction(new MultivariateFunction() {
                        @Override
                        public double value(double[] point) {
                            return objFun.apply(point);
                        }
                    }),
                    GoalType.MINIMIZE,
                    new SimpleBounds(lb, ub),
                    new InitialGuess(x0Safe)
            );
            return new Pair<double[], Double>(result.getPoint(), result.getValue());
        } catch (Exception e) {
            // If BOBYQA fails, try coordinate descent
            return coordinateDescent(objFun, x0Safe, lb, ub, maxIter);
        }
    }

    public static Pair<double[], Double> fmincon(java.util.function.Function<double[], Double> objFun,
                                                 double[] x0, double[] lb, double[] ub) {
        return fmincon(objFun, x0, lb, ub, 10000);
    }

    /**
     * Simple coordinate descent fallback optimizer.
     */
    private static Pair<double[], Double> coordinateDescent(java.util.function.Function<double[], Double> objFun,
                                                            double[] x0, double[] lb, double[] ub,
                                                            int maxIter) {
        int n = x0.length;
        double[] x = x0.clone();
        double fBest = objFun.apply(x);

        for (int iter = 0; iter < maxIter; iter++) {
            boolean improved = false;
            for (int i = 0; i < n; i++) {
                double delta = Math.max(1e-6, Math.abs(x[i]) * 0.01);
                int[] signs = new int[]{-1, 1};
                for (int sign : signs) {
                    double[] xTrial = x.clone();
                    xTrial[i] = Math.max(lb[i], Math.min(ub[i], x[i] + sign * delta));
                    double fTrial = objFun.apply(xTrial);
                    if (fTrial < fBest) {
                        x[i] = xTrial[i];
                        fBest = fTrial;
                        improved = true;
                    }
                }
            }
            if (!improved) break;
        }
        return new Pair<double[], Double>(x, fBest);
    }

    /**
     * Quadratic program solver with non-negative constraints.
     * Solves: min 0.5 x'Hx + f'x  subject to x &gt;= lb.
     */
    public static Pair<Matrix, Double> quadprog(Matrix H, Matrix f, Matrix lb) {
        int n = H.getNumRows();

        // Shift to handle non-zero lower bounds: y = x - lb
        Matrix fShifted = new Matrix(n, 1);
        for (int i = 0; i < n; i++) {
            double sum = f.get(i, 0);
            for (int j = 0; j < n; j++) {
                sum += H.get(i, j) * lb.get(j, 0);
            }
            fShifted.set(i, 0, sum);
        }

        // Cholesky decomposition H = L L'
        Matrix L = choleskyDecompose(H);

        if (L != null) {
            // Convert QP to NNLS: min ||L'y - b||^2 s.t. y >= 0, where b = -L^{-1}f
            Matrix LInvF = forwardSolve(L, fShifted);
            Matrix rhs = new Matrix(n, 1);
            for (int i = 0; i < n; i++) {
                rhs.set(i, 0, -LInvF.get(i, 0));
            }
            // Use L transpose (H = LL', objective = 0.5||L'y + L^{-1}f||^2)
            Matrix Lt = new Matrix(n, n);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    Lt.set(i, j, L.get(j, i));
                }
            }
            Matrix y = NnlsSolver.lsqnonneg(Lt, rhs);

            // x = y + lb
            Matrix x = new Matrix(n, 1);
            for (int i = 0; i < n; i++) {
                x.set(i, 0, y.get(i, 0) + lb.get(i, 0));
            }

            // Compute objective value
            double objVal = 0.0;
            for (int i = 0; i < n; i++) {
                objVal += f.get(i, 0) * x.get(i, 0);
                for (int j = 0; j < n; j++) {
                    objVal += 0.5 * x.get(i, 0) * H.get(i, j) * x.get(j, 0);
                }
            }
            return new Pair<Matrix, Double>(x, objVal);
        } else {
            // Fallback: use projected gradient descent
            return projectedGradientQP(H, f, lb, n);
        }
    }

    private static Pair<Matrix, Double> projectedGradientQP(Matrix H, Matrix f, Matrix lb, int n) {
        Matrix x = new Matrix(n, 1);
        for (int i = 0; i < n; i++) {
            x.set(i, 0, lb.get(i, 0) + 1.0);
        }

        int maxIter = 1000;
        double stepSize = 0.001;

        for (int iter = 0; iter < maxIter; iter++) {
            // gradient = Hx + f
            Matrix grad = new Matrix(n, 1);
            for (int i = 0; i < n; i++) {
                double sum = f.get(i, 0);
                for (int j = 0; j < n; j++) {
                    sum += H.get(i, j) * x.get(j, 0);
                }
                grad.set(i, 0, sum);
            }

            // Project: x = max(lb, x - step * grad)
            boolean changed = false;
            for (int i = 0; i < n; i++) {
                double newVal = Math.max(lb.get(i, 0), x.get(i, 0) - stepSize * grad.get(i, 0));
                if (Math.abs(newVal - x.get(i, 0)) > 1e-12) changed = true;
                x.set(i, 0, newVal);
            }
            if (!changed) break;
        }

        double objVal = 0.0;
        for (int i = 0; i < n; i++) {
            objVal += f.get(i, 0) * x.get(i, 0);
            for (int j = 0; j < n; j++) {
                objVal += 0.5 * x.get(i, 0) * H.get(i, j) * x.get(j, 0);
            }
        }
        return new Pair<Matrix, Double>(x, objVal);
    }

    private static Matrix choleskyDecompose(Matrix A) {
        int n = A.getNumRows();
        Matrix L = new Matrix(n, n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                double sum = 0.0;
                for (int k = 0; k < j; k++) {
                    sum += L.get(i, k) * L.get(j, k);
                }
                if (i == j) {
                    double diag = A.get(i, i) - sum;
                    if (diag <= 0) return null; // Not positive definite
                    L.set(i, j, Math.sqrt(diag));
                } else {
                    double lJJ = L.get(j, j);
                    if (Math.abs(lJJ) < 1e-14) return null;
                    L.set(i, j, (A.get(i, j) - sum) / lJJ);
                }
            }
        }
        return L;
    }

    private static Matrix forwardSolve(Matrix L, Matrix b) {
        int n = L.getNumRows();
        Matrix x = new Matrix(n, 1);
        for (int i = 0; i < n; i++) {
            double sum = b.get(i, 0);
            for (int j = 0; j < i; j++) {
                sum -= L.get(i, j) * x.get(j, 0);
            }
            x.set(i, 0, sum / L.get(i, i));
        }
        return x;
    }
}
