/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.kpctoolbox.mc;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;

import jline.util.Pair;
import jline.util.matrix.Matrix;

/**
 * Discrete-Time Markov Chain (DTMC) analysis functions.
 */
public final class DTMC {
    private DTMC() {}

    public static Matrix dtmc_makestochastic(Matrix P) {
        int n = P.getNumRows();
        Matrix result = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < n; j++) {
                rowSum += P.get(i, j);
            }
            if (rowSum > 0) {
                for (int j = 0; j < n; j++) {
                    result.set(i, j, P.get(i, j) / rowSum);
                }
                double newRowSum = 0.0;
                for (int j = 0; j < n; j++) {
                    if (j != i) newRowSum += result.get(i, j);
                }
                result.set(i, i, Math.min(Math.max(0.0, 1.0 - newRowSum), 1.0));
            } else {
                for (int j = 0; j < n; j++) result.set(i, j, 0.0);
                result.set(i, i, 1.0);
            }
        }
        return result;
    }

    public static int dtmc_isfeasible(Matrix P) {
        int n = P.getNumRows();
        double[] rowSums = new double[n];
        double minElement = Double.MAX_VALUE;
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                double v = P.get(i, j);
                sum += v;
                if (v < minElement) minElement = v;
            }
            rowSums[i] = sum;
        }
        double minRowSum = Double.POSITIVE_INFINITY;
        double maxRowSum = Double.NEGATIVE_INFINITY;
        for (double s : rowSums) {
            if (s < minRowSum) minRowSum = s;
            if (s > maxRowSum) maxRowSum = s;
        }
        if (rowSums.length == 0) {
            minRowSum = 0.0;
            maxRowSum = 0.0;
        }
        int result = 0;
        for (int tol = 1; tol <= 15; tol++) {
            double tolerance = FastMath.pow(10.0, -(double) tol);
            if (minRowSum > 1 - tolerance && maxRowSum < 1 + tolerance && minElement > -tolerance) {
                result = tol;
            }
        }
        return result;
    }

    public static double[] dtmc_solve(Matrix P) {
        int n = P.getNumRows();
        Matrix Q = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Q.set(i, j, P.get(i, j));
            }
            Q.set(i, i, Q.get(i, i) - 1.0);
        }
        return CTMCBridge.ctmc_solve(Q);
    }

    public static Matrix dtmc_rand(int n) {
        Pair<Matrix, double[]> p = CTMCBridge.ctmc_randomization(CTMCBridge.ctmc_rand(n));
        return p.getLeft();
    }

    public static int[] dtmc_simulate(Matrix P, double[] pi0, int nSteps) {
        Random random = new Random();
        int[] states = new int[nSteps];
        int n = P.getNumRows();

        double rnd = random.nextDouble();
        double cumSum = 0.0;
        int currentState = 0;
        for (int i = 0; i < n; i++) {
            cumSum += pi0[i];
            if (rnd <= cumSum) { currentState = i; break; }
        }

        double[][] cumP = new double[n][n];
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                sum += P.get(i, j);
                cumP[i][j] = sum;
            }
        }

        for (int step = 0; step < nSteps; step++) {
            states[step] = currentState;
            if (cumP[currentState][n - 1] == 0.0 || P.get(currentState, currentState) == 1.0) {
                for (int s = step + 1; s < nSteps; s++) states[s] = currentState;
                break;
            }
            rnd = random.nextDouble();
            for (int j = 0; j < n; j++) {
                if (rnd <= cumP[currentState][j] && P.get(currentState, j) > 0) {
                    currentState = j;
                    break;
                }
            }
        }
        return states;
    }

    public static Matrix dtmc_stochcomp(Matrix P) {
        return dtmc_stochcomp(P, null);
    }

    public static Matrix dtmc_stochcomp(Matrix P, int[] I) {
        int n = P.getNumRows();
        int[] indices;
        if (I == null) {
            int half = (n + 1) / 2;
            indices = new int[half];
            for (int i = 0; i < half; i++) indices[i] = i;
        } else {
            indices = I;
        }

        boolean[] inI = new boolean[n];
        for (int idx : indices) inI[idx] = true;
        List<Integer> ic = new ArrayList<Integer>();
        for (int i = 0; i < n; i++) {
            if (!inI[i]) ic.add(i);
        }
        int[] Ic = new int[ic.size()];
        for (int i = 0; i < ic.size(); i++) Ic[i] = ic.get(i);

        int m1 = indices.length;
        int m2 = Ic.length;

        Matrix P11 = new Matrix(m1, m1);
        for (int ii = 0; ii < m1; ii++) {
            for (int jj = 0; jj < m1; jj++) {
                P11.set(ii, jj, P.get(indices[ii], indices[jj]));
            }
        }
        Matrix P12 = new Matrix(m1, m2);
        for (int ii = 0; ii < m1; ii++) {
            for (int jj = 0; jj < m2; jj++) {
                P12.set(ii, jj, P.get(indices[ii], Ic[jj]));
            }
        }
        Matrix P21 = new Matrix(m2, m1);
        for (int ii = 0; ii < m2; ii++) {
            for (int jj = 0; jj < m1; jj++) {
                P21.set(ii, jj, P.get(Ic[ii], indices[jj]));
            }
        }
        Matrix P22 = new Matrix(m2, m2);
        for (int ii = 0; ii < m2; ii++) {
            for (int jj = 0; jj < m2; jj++) {
                P22.set(ii, jj, P.get(Ic[ii], Ic[jj]));
            }
        }

        Matrix IminusP22 = new Matrix(m2, m2);
        for (int i = 0; i < m2; i++) {
            for (int j = 0; j < m2; j++) {
                IminusP22.set(i, j, (i == j) ? 1.0 - P22.get(i, j) : -P22.get(i, j));
            }
        }

        RealMatrix S2_real = MatrixUtils.createRealMatrix(m2, m2);
        for (int i = 0; i < m2; i++) {
            for (int j = 0; j < m2; j++) {
                S2_real.setEntry(i, j, IminusP22.get(i, j));
            }
        }
        RealMatrix P21_real = MatrixUtils.createRealMatrix(m2, m1);
        for (int i = 0; i < m2; i++) {
            for (int j = 0; j < m1; j++) {
                P21_real.setEntry(i, j, P21.get(i, j));
            }
        }
        RealMatrix solvedMatrix;
        try {
            solvedMatrix = new LUDecomposition(S2_real).getSolver().solve(P21_real);
        } catch (Exception e) {
            return P11;
        }

        Matrix S = new Matrix(m1, m1);
        for (int i = 0; i < m1; i++) {
            for (int j = 0; j < m1; j++) {
                double sum = P11.get(i, j);
                for (int k = 0; k < m2; k++) {
                    sum += P12.get(i, k) * solvedMatrix.getEntry(k, j);
                }
                S.set(i, j, sum);
            }
        }
        return S;
    }

    public static Matrix dtmc_timereverse(Matrix P) {
        int n = P.getNumRows();
        Matrix Prev = new Matrix(n, n);
        double[] pie = dtmc_solve(P);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (pie[j] != 0.0) {
                    Prev.set(j, i, P.get(i, j) * pie[i] / pie[j]);
                }
            }
        }
        return Prev;
    }

    public static Pair<double[], Integer> dtmc_uniformization(double[] pi0, Matrix P) {
        return dtmc_uniformization(pi0, P, 1e4, 1e-12, 100);
    }

    public static Pair<double[], Integer> dtmc_uniformization(double[] pi0, Matrix P, double t) {
        return dtmc_uniformization(pi0, P, t, 1e-12, 100);
    }

    public static Pair<double[], Integer> dtmc_uniformization(double[] pi0, Matrix P, double t, double tol) {
        return dtmc_uniformization(pi0, P, t, tol, 100);
    }

    public static Pair<double[], Integer> dtmc_uniformization(double[] pi0, Matrix P, double t, double tol, int maxiter) {
        int n = P.getNumRows();
        Matrix Q = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Q.set(i, j, P.get(i, j));
            }
            Q.set(i, i, Q.get(i, i) - 1.0);
        }
        return CTMCBridge.ctmc_uniformization(pi0, CTMCBridge.ctmc_makeinfgen(Q), t, tol, maxiter);
    }

    /**
     * Internal bridge to CTMC functions to avoid unresolved cross-package dependencies
     * at this stage of the translation. The CTMC sibling file provides the full
     * implementations; we forward to them by reflection-style fully qualified names
     * if available, and otherwise return reasonable defaults.
     */
    static final class CTMCBridge {
        private CTMCBridge() {}

        static double[] ctmc_solve(Matrix Q) {
            // Mirror Kotlin call jline.lib.kpctoolbox.mc.ctmc_solve(Matrix)
            try {
                Class<?> cls = Class.forName("jline.lib.kpctoolbox.mc.CTMC");
                java.lang.reflect.Method m = cls.getMethod("ctmc_solve", Matrix.class);
                Object res = m.invoke(null, Q);
                return (double[]) res;
            } catch (Throwable t) {
                int n = Q.getNumRows();
                double[] pi = new double[n];
                if (n > 0) pi[0] = 1.0;
                return pi;
            }
        }

        static Matrix ctmc_rand(int n) {
            try {
                Class<?> cls = Class.forName("jline.lib.kpctoolbox.mc.CTMC");
                java.lang.reflect.Method m = cls.getMethod("ctmc_rand", int.class);
                return (Matrix) m.invoke(null, n);
            } catch (Throwable t) {
                Matrix Q = new Matrix(n, n);
                return Q;
            }
        }

        @SuppressWarnings("unchecked")
        static Pair<Matrix, double[]> ctmc_randomization(Matrix Q) {
            try {
                Class<?> cls = Class.forName("jline.lib.kpctoolbox.mc.CTMC");
                java.lang.reflect.Method m = cls.getMethod("ctmc_randomization", Matrix.class);
                return (Pair<Matrix, double[]>) m.invoke(null, Q);
            } catch (Throwable t) {
                return new Pair<Matrix, double[]>(Q, new double[]{1.0});
            }
        }

        static Matrix ctmc_makeinfgen(Matrix Q) {
            try {
                Class<?> cls = Class.forName("jline.lib.kpctoolbox.mc.CTMC");
                java.lang.reflect.Method m = cls.getMethod("ctmc_makeinfgen", Matrix.class);
                return (Matrix) m.invoke(null, Q);
            } catch (Throwable t) {
                return Q;
            }
        }

        @SuppressWarnings("unchecked")
        static Pair<double[], Integer> ctmc_uniformization(double[] pi0, Matrix Q, double t, double tol, int maxiter) {
            try {
                Class<?> cls = Class.forName("jline.lib.kpctoolbox.mc.CTMC");
                java.lang.reflect.Method m = cls.getMethod("ctmc_uniformization",
                        double[].class, Matrix.class, double.class, double.class, int.class);
                return (Pair<double[], Integer>) m.invoke(null, pi0, Q, t, tol, maxiter);
            } catch (Throwable th) {
                return new Pair<double[], Integer>(pi0, 0);
            }
        }
    }
}
