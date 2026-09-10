/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.fitting;

import jline.lib.butools.ph.PHRepresentation;
import jline.util.Pair;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.List;

public final class PHFromTrace {
    private PHFromTrace() {}

    public static class PHFitResult {
        public final Matrix alpha;
        public final Matrix A;
        public final double logli;
        public PHFitResult(Matrix alpha, Matrix A, double logli) {
            this.alpha = alpha;
            this.A = A;
            this.logli = logli;
        }
    }

    public static PHFitResult phFromTrace(double[] trace, int[] orders) {
        return phFromTrace(trace, orders, 200, 1e-7, null, null);
    }

    public static PHFitResult phFromTrace(double[] trace, int[] orders, int maxIter) {
        return phFromTrace(trace, orders, maxIter, 1e-7, null, null);
    }

    public static PHFitResult phFromTrace(double[] trace, int[] orders, int maxIter, double stopCond) {
        return phFromTrace(trace, orders, maxIter, stopCond, null, null);
    }

    public static PHFitResult phFromTrace(double[] trace, int[] orders, int maxIter, double stopCond,
                                          PHRepresentation initialGuess) {
        return phFromTrace(trace, orders, maxIter, stopCond, initialGuess, null);
    }

    public static PHFitResult phFromTrace(double[] trace, int[] orders, int maxIter, double stopCond,
                                          PHRepresentation initialGuess, String resultFormat) {
        Pair<PHRepresentation, Double> result = phFromTraceInternal(trace, orders, maxIter, stopCond, initialGuess);
        return new PHFitResult(result.getLeft().alpha, result.getLeft().A, result.getRight());
    }

    public static PHFitResult phFromTrace(double[] trace, int totalOrder) {
        return phFromTrace(trace, totalOrder, 200, 1e-7);
    }

    public static PHFitResult phFromTrace(double[] trace, int totalOrder, int maxIter) {
        return phFromTrace(trace, totalOrder, maxIter, 1e-7);
    }

    public static PHFitResult phFromTrace(double[] trace, int totalOrder, int maxIter, double stopCond) {
        Pair<PHRepresentation, Double> result = phFromTraceAllOrders(trace, totalOrder, maxIter, stopCond);
        return new PHFitResult(result.getLeft().alpha, result.getLeft().A, result.getRight());
    }

    private static List<int[]> allOrders(int branches, int sumOrders) {
        List<int[]> result = new ArrayList<int[]>();
        if (branches == 1) {
            result.add(new int[]{sumOrders});
            return result;
        }
        for (int i = 0; i < sumOrders - branches + 1; i++) {
            List<int[]> subPartitions = allOrders(branches - 1, sumOrders - i - 1);
            for (int[] sub : subPartitions) {
                int[] combined = new int[sub.length + 1];
                System.arraycopy(sub, 0, combined, 0, sub.length);
                combined[sub.length] = i + 1;
                java.util.Arrays.sort(combined);
                boolean duplicate = false;
                for (int[] existing : result) {
                    if (java.util.Arrays.equals(existing, combined)) {
                        duplicate = true; break;
                    }
                }
                if (!duplicate) result.add(combined);
            }
        }
        return result;
    }

    private static Pair<PHRepresentation, Double> phFromTraceAllOrders(
            double[] trace, int totalOrder, int maxIter, double stopCond) {
        PHRepresentation bestPH = null;
        double bestLogLi = Double.NEGATIVE_INFINITY;

        for (int br = 2; br <= totalOrder; br++) {
            List<int[]> allOrd = allOrders(br, totalOrder);
            for (int[] ordk : allOrd) {
                Pair<PHRepresentation, Double> res = phFromTraceInternal(trace, ordk, maxIter, stopCond, null);
                if (res.getRight() > bestLogLi) {
                    bestLogLi = res.getRight();
                    bestPH = res.getLeft();
                }
            }
        }
        if (bestPH == null) {
            return phFromTraceInternal(trace, new int[]{totalOrder}, maxIter, stopCond, null);
        }
        return new Pair<PHRepresentation, Double>(bestPH, bestLogLi);
    }

    private static Pair<PHRepresentation, Double> phFromTraceInternal(
            double[] trace, int[] orders, int maxIter, double stopCond, PHRepresentation initialGuess) {
        int M = orders.length;
        int K = trace.length;

        double[] alphav = new double[M];
        double[] lambd = new double[M];

        if (initialGuess != null) {
            int alphaLen = initialGuess.alpha.length();
            if (alphaLen != M) {
                throw new IllegalArgumentException("The length of the initial vector (" + alphaLen
                        + ") is not consistent with the number of branches (" + M + ")!");
            }
            for (int i = 0; i < M; i++) alphav[i] = initialGuess.alpha.get(0, i);
            for (int i = 0; i < M; i++) lambd[i] = -initialGuess.A.get(i, i);
        } else {
            double traceSum = 0.0;
            for (double tr : trace) traceSum += tr;
            double traceMean = traceSum / K;

            for (int i = 0; i < M; i++) alphav[i] = 1.0 / M;
            for (int i = 0; i < M; i++) lambd[i] = (double) orders[i] * (i + 1);

            double inim = 0.0;
            for (int i = 0; i < M; i++) inim += alphav[i] / (i + 1.0);

            double scaleFactor = inim / traceMean;
            for (int i = 0; i < M; i++) lambd[i] *= scaleFactor;
        }

        double[] W = new double[K];
        for (int i = 0; i < K; i++) W[i] = 1.0 / K;

        double[][] Q = new double[M][K];

        double[] logFactorials = new double[M];
        for (int i = 0; i < M; i++) {
            double lf = 0.0;
            for (int j = 1; j < orders[i]; j++) lf += Math.log(j);
            logFactorials[i] = lf;
        }

        double logli = 1e-14;
        double ologli = 1.0;
        int steps = 1;

        while (Math.abs((ologli - logli) / logli) > stopCond && steps <= maxIter) {
            ologli = logli;

            for (int i = 0; i < M; i++) {
                int ord = orders[i];
                double lam = lambd[i];
                double logAlpha = Math.log(alphav[i]);
                double logLam = Math.log(lam);
                double logFact = logFactorials[i];
                for (int k = 0; k < K; k++) {
                    double t = trace[k];
                    double logQ = logAlpha + (double) ord * logLam + (ord - 1) * Math.log(t) - logFact - lam * t;
                    Q[i][k] = Math.exp(logQ);
                }
            }

            double[] nor = new double[K];
            for (int k = 0; k < K; k++) {
                double s = 0.0;
                for (int i = 0; i < M; i++) s += Q[i][k];
                nor[k] = s;
            }
            for (int i = 0; i < M; i++) {
                for (int k = 0; k < K; k++) {
                    if (nor[k] > 0.0) Q[i][k] /= nor[k];
                }
            }

            logli = 0.0;
            for (int k = 0; k < K; k++) {
                if (nor[k] > 0.0) logli += Math.log(nor[k]) * W[k];
            }

            double[] v1 = new double[M];
            double[] v2 = new double[M];
            for (int i = 0; i < M; i++) {
                double s1 = 0.0, s2 = 0.0;
                for (int k = 0; k < K; k++) {
                    s1 += Q[i][k] * W[k];
                    s2 += Q[i][k] * trace[k] * W[k];
                }
                v1[i] = s1;
                v2[i] = s2;
            }
            for (int i = 0; i < M; i++) {
                alphav[i] = v1[i];
                if (v2[i] > 0.0) lambd[i] = (double) orders[i] * v1[i] / v2[i];
            }
            steps++;
        }

        int N = 0;
        for (int o : orders) N += o;
        Matrix alpha = new Matrix(1, N);
        Matrix A = new Matrix(N, N);

        int ix = 0;
        for (int i = 0; i < M; i++) {
            int ord = orders[i];
            double lam = lambd[i];
            alpha.set(0, ix, alphav[i]);
            if (ord == 1) {
                A.set(ix, ix, -lam);
            } else {
                for (int j = 0; j < ord; j++) {
                    A.set(ix + j, ix + j, -lam);
                    if (j < ord - 1) A.set(ix + j, ix + j + 1, lam);
                }
            }
            ix += ord;
        }
        return new Pair<PHRepresentation, Double>(new PHRepresentation(alpha, A), logli);
    }
}
