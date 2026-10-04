/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.fitting;

import jline.lib.butools.mc.DTMCSolve;
import jline.util.matrix.Matrix;
import jline.util.Pair;

import java.util.ArrayList;
import java.util.List;

/**
 * Top-level functions for MAP-from-trace fitting.
 */
public final class MAPFromTrace {

    private MAPFromTrace() {
    }

    /**
     * Generates all distinct sorted partitions of sumorders into branches parts,
     * each part >= 1.
     */
    private static List<List<Integer>> allOrders(int branches, int sumorders) {
        if (branches == 1) {
            List<List<Integer>> result = new ArrayList<List<Integer>>();
            List<Integer> single = new ArrayList<Integer>();
            single.add(sumorders);
            result.add(single);
            return result;
        }
        List<List<Integer>> result = new ArrayList<List<Integer>>();
        for (int i = 0; i < sumorders - branches + 1; i++) {
            List<List<Integer>> subPartitions = allOrders(branches - 1, sumorders - i - 1);
            for (List<Integer> sub : subPartitions) {
                List<Integer> candidate = new ArrayList<Integer>(sub.size() + 1);
                candidate.addAll(sub);
                candidate.add(i + 1);
                java.util.Collections.sort(candidate);
                if (!result.contains(candidate)) {
                    result.add(candidate);
                }
            }
        }
        return result;
    }

    private static double factorial(int n) {
        double f = 1.0;
        for (int i = 2; i <= n; i++) {
            f *= (double) i;
        }
        return f;
    }

    public static MAPFitResult mapFromTrace(double[] trace, int[] orders) {
        return mapFromTrace(trace, orders, 200, 1e-7, null, null);
    }

    public static MAPFitResult mapFromTrace(double[] trace, int[] orders, int maxIter) {
        return mapFromTrace(trace, orders, maxIter, 1e-7, null, null);
    }

    public static MAPFitResult mapFromTrace(double[] trace, int[] orders, int maxIter, double stopCond) {
        return mapFromTrace(trace, orders, maxIter, stopCond, null, null);
    }

    public static MAPFitResult mapFromTrace(double[] trace, int[] orders, int maxIter, double stopCond, Matrix[] initialGuess) {
        return mapFromTrace(trace, orders, maxIter, stopCond, initialGuess, null);
    }

    /**
     * Performs MAP fitting using the EM algorithm (ErCHMM).
     */
    public static MAPFitResult mapFromTrace(double[] trace, int[] orders, int maxIter, double stopCond,
                                             Matrix[] initialGuess, String resultFormat) {
        if (orders.length == 1) {
            int totalStates = orders[0];
            Triple res = mapFromTraceAllOrders(trace, totalStates, maxIter, stopCond, initialGuess);
            return new MAPFitResult(res.first, res.second, res.third);
        }

        Triple result = mapFromTraceEM(trace, orders, maxIter, stopCond, initialGuess);
        return new MAPFitResult(result.first, result.second, result.third);
    }

    private static class Triple {
        Matrix first;
        Matrix second;
        double third;
        Triple(Matrix a, Matrix b, double c) {
            this.first = a;
            this.second = b;
            this.third = c;
        }
    }

    private static Triple mapFromTraceAllOrders(double[] trace, int totalStates, int maxIter, double stopCond, Matrix[] initialGuess) {
        Matrix bestD0 = null;
        Matrix bestD1 = null;
        double bestLogli = Double.NEGATIVE_INFINITY;

        for (int br = 2; br <= totalStates; br++) {
            List<List<Integer>> allOrd = allOrders(br, totalStates);
            for (List<Integer> ordk : allOrd) {
                int[] ordArray = new int[ordk.size()];
                for (int i = 0; i < ordk.size(); i++) ordArray[i] = ordk.get(i);
                Triple result = mapFromTraceEM(trace, ordArray, maxIter, stopCond, initialGuess);
                if (result.third > bestLogli) {
                    bestLogli = result.third;
                    bestD0 = result.first;
                    bestD1 = result.second;
                }
            }
        }
        return new Triple(bestD0, bestD1, bestLogli);
    }

    private static Triple mapFromTraceEM(double[] trace, int[] orders, int maxIter, double stopCond, Matrix[] initialGuess) {
        int M = orders.length;
        int K = trace.length;

        double[] alphav = new double[M];
        double[] lambd = new double[M];
        double[][] P = new double[M][M];

        if (initialGuess == null) {
            double trm = 0.0;
            for (double t : trace) trm += t;
            trm /= (double) K;
            for (int m = 0; m < M; m++) {
                alphav[m] = 1.0 / M;
            }
            double inim = 0.0;
            for (int m = 0; m < M; m++) {
                double linVal = (double) (m + 1);
                lambd[m] = (double) orders[m] * linVal;
                inim += alphav[m] / linVal;
            }
            for (int m = 0; m < M; m++) {
                lambd[m] = lambd[m] * inim / trm;
            }
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    P[i][j] = alphav[j];
                }
            }
        } else {
            if (initialGuess.length != 2) {
                throw new IllegalArgumentException("Initial guess must contain exactly 2 matrices (D0 and D1)");
            }
            double[] initLambd = new double[M];
            double[][] initP = new double[M][M];
            extractErCHMMParams(initialGuess[0], initialGuess[1], orders, initLambd, initP);
            for (int m = 0; m < M; m++) {
                lambd[m] = initLambd[m];
            }
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    P[i][j] = initP[i][j];
                }
            }
            Matrix Pmat = new Matrix(M, M);
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    Pmat.set(i, j, P[i][j]);
                }
            }
            Matrix piVec = DTMCSolve.dtmcSolve(Pmat);
            for (int m = 0; m < M; m++) {
                alphav[m] = piVec.get(0, m);
            }
        }

        double[][] Q = new double[M][K];
        double[][] A = new double[K][M];
        double[] Ascale = new double[K];
        double[][] B = new double[K][M];
        double[] Bscale = new double[K];

        double logli = 1e-14;
        double ologli = 0.0;
        int steps = 1;

        while (Math.abs((ologli - logli) / logli) > stopCond && steps <= maxIter) {
            ologli = logli;

            for (int m = 0; m < M; m++) {
                double lm = lambd[m];
                int om = orders[m];
                double factVal = factorial(om - 1);
                for (int k = 0; k < K; k++) {
                    double lt = lm * trace[k];
                    Q[m][k] = Math.pow(lt, om - 1) / factVal * lm * Math.exp(-lt);
                }
            }

            double[] prev = alphav.clone();
            double scprev = 0.0;
            for (int k = 0; k < K; k++) {
                double[] temp = new double[M];
                for (int m = 0; m < M; m++) {
                    temp[m] = prev[m] * Q[m][k];
                }
                for (int m = 0; m < M; m++) {
                    double s = 0.0;
                    for (int j = 0; j < M; j++) {
                        s += temp[j] * P[j][m];
                    }
                    prev[m] = s;
                }
                double sumPrev = 0.0;
                for (int m = 0; m < M; m++) sumPrev += prev[m];
                double scale = log2(sumPrev);
                double invScale = Math.pow(2.0, -scale);
                for (int m = 0; m < M; m++) prev[m] *= invScale;
                Ascale[k] = scprev + scale;
                for (int m = 0; m < M; m++) A[k][m] = prev[m];
                scprev = Ascale[k];
            }

            double[] nnext = new double[M];
            for (int m = 0; m < M; m++) nnext[m] = 1.0;
            scprev = 0.0;
            for (int k = K - 1; k >= 0; k--) {
                double[] Pnnext = new double[M];
                for (int i = 0; i < M; i++) {
                    double s = 0.0;
                    for (int j = 0; j < M; j++) {
                        s += P[i][j] * nnext[j];
                    }
                    Pnnext[i] = s;
                }
                for (int m = 0; m < M; m++) {
                    nnext[m] = Q[m][k] * Pnnext[m];
                }
                double sumNext = 0.0;
                for (int m = 0; m < M; m++) sumNext += nnext[m];
                double scale = log2(sumNext);
                double invScale = Math.pow(2.0, -scale);
                for (int m = 0; m < M; m++) nnext[m] *= invScale;
                Bscale[k] = scprev + scale;
                for (int m = 0; m < M; m++) B[k][m] = nnext[m];
                scprev = Bscale[k];
            }

            double llh = 0.0;
            for (int m = 0; m < M; m++) {
                llh += alphav[m] * B[0][m];
            }
            logli = (Math.log(llh) + Bscale[0] * Math.log(2.0)) / K;
            double illh = 1.0 / llh;

            double[] oldAlphav = alphav.clone();

            double[][] AB = new double[K][M];
            for (int k = 0; k < K; k++) {
                for (int m = 0; m < M; m++) {
                    double avKM = (k == 0) ? oldAlphav[m] : A[k - 1][m];
                    AB[k][m] = avKM * B[k][m];
                }
            }
            for (int k = 0; k < K; k++) {
                double rowSum = 0.0;
                for (int m = 0; m < M; m++) rowSum += AB[k][m];
                if (rowSum > 0.0) {
                    for (int m = 0; m < M; m++) AB[k][m] /= rowSum;
                }
            }
            double[] v1 = new double[M];
            for (int m = 0; m < M; m++) {
                double s = 0.0;
                for (int k = 0; k < K; k++) s += AB[k][m];
                v1[m] = s;
            }
            double[] v2 = new double[M];
            for (int m = 0; m < M; m++) {
                double s = 0.0;
                for (int k = 0; k < K; k++) s += AB[k][m] * trace[k];
                v2[m] = s;
            }
            for (int m = 0; m < M; m++) {
                alphav[m] = v1[m] / K;
                lambd[m] = (double) orders[m] * v1[m] / v2[m];
            }

            double[][] Avv = new double[K][M];
            for (int k = 0; k < K; k++) {
                for (int m = 0; m < M; m++) {
                    double avKM = (k == 0) ? oldAlphav[m] : A[k - 1][m];
                    Avv[k][m] = avKM * Q[m][k];
                }
            }
            double[] nor = new double[K];
            for (int k = 0; k < K; k++) {
                double ascaleV = (k == 0) ? 0.0 : Ascale[k - 1];
                double bscaleV = (k < K - 1) ? Bscale[k + 1] : 0.0;
                nor[k] = illh * Math.pow(2.0, ascaleV + bscaleV - Bscale[0]);
            }
            for (int k = 0; k < K; k++) {
                for (int m = 0; m < M; m++) Avv[k][m] *= nor[k];
            }
            double[][] AvvT_BvT = new double[M][M];
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    double s = 0.0;
                    for (int k = 0; k < K; k++) {
                        double bvTkj = (k < K - 1) ? B[k + 1][j] : 1.0;
                        s += Avv[k][i] * bvTkj;
                    }
                    AvvT_BvT[i][j] = s;
                }
            }
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    P[i][j] = AvvT_BvT[i][j] * P[i][j];
                }
            }
            for (int i = 0; i < M; i++) {
                double rowSum = 0.0;
                for (int j = 0; j < M; j++) rowSum += P[i][j];
                if (rowSum > 0.0) {
                    for (int j = 0; j < M; j++) P[i][j] /= rowSum;
                }
            }

            steps++;
        }

        int N = 0;
        for (int o : orders) N += o;
        Matrix D0 = new Matrix(N, N);
        int ix = 0;
        for (int i = 0; i < M; i++) {
            if (orders[i] == 1) {
                D0.set(ix, ix, -lambd[i]);
            } else {
                for (int p = 0; p < orders[i]; p++) {
                    D0.set(ix + p, ix + p, -lambd[i]);
                }
                for (int p = 0; p < orders[i] - 1; p++) {
                    D0.set(ix + p, ix + p + 1, lambd[i]);
                }
            }
            ix += orders[i];
        }

        Matrix D1 = new Matrix(N, N);
        int[] indicesTo = new int[M];
        indicesTo[0] = 0;
        for (int j = 1; j < M; j++) {
            indicesTo[j] = indicesTo[j - 1] + orders[j - 1];
        }
        int[] indicesFrom = new int[M];
        int cumSum = 0;
        for (int i = 0; i < M; i++) {
            cumSum += orders[i];
            indicesFrom[i] = cumSum - 1;
        }
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                D1.set(indicesFrom[i], indicesTo[j], lambd[i] * P[i][j]);
            }
        }
        return new Triple(D0, D1, logli);
    }

    private static double log2(double x) {
        return Math.log(x) / Math.log(2.0);
    }

    private static void extractErCHMMParams(Matrix D0, Matrix D1, int[] orders, double[] outLambd, double[][] outP) {
        int M = orders.length;
        int ix = 0;
        for (int i = 0; i < M; i++) {
            outLambd[i] = -D0.get(ix, ix);
            ix += orders[i];
        }

        int[] indicesTo = new int[M];
        indicesTo[0] = 0;
        for (int j = 1; j < M; j++) {
            indicesTo[j] = indicesTo[j - 1] + orders[j - 1];
        }
        int[] indicesFrom = new int[M];
        int cumSum = 0;
        for (int i = 0; i < M; i++) {
            cumSum += orders[i];
            indicesFrom[i] = cumSum - 1;
        }

        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                if (outLambd[i] > 0.0) {
                    outP[i][j] = D1.get(indicesFrom[i], indicesTo[j]) / outLambd[i];
                } else {
                    outP[i][j] = 0.0;
                }
            }
        }
    }
}
