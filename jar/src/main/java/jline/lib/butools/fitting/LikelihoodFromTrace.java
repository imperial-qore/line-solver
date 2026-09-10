/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.fitting;

import java.util.Arrays;

import jline.lib.butools.mc.DTMCSolve;
import jline.util.matrix.Matrix;

public final class LikelihoodFromTrace {
    private LikelihoodFromTrace() {}

    /**
     * Evaluates the log-likelihood of a trace with the given PH distribution.
     * The result is divided by the length of the trace.
     */
    public static double likelihoodFromTracePH(double[] trace, Matrix alpha, Matrix A, double prec) {
        double[] sortedTrace = trace.clone();
        Arrays.sort(sortedTrace);
        int N = A.getNumRows();

        // Find lambda as max absolute value of diagonal
        double lambda = 0.0;
        for (int i = 0; i < N; i++) {
            double absVal = Math.abs(A.get(i, i));
            if (absVal > lambda) lambda = absVal;
        }

        // P = A/lambda + I
        Matrix P = A.scale(1.0 / lambda).add(Matrix.eye(N));

        // a = -A * ones (closing vector)
        Matrix a = new Matrix(N, 1);
        for (int i = 0; i < N; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < N; j++) {
                rowSum += A.get(i, j);
            }
            a.set(i, 0, -rowSum);
        }

        double eps = Math.max(prec, Math.pow(10.0, Math.log10(prec) + Math.log10(lambda)));
        int L = sortedTrace.length;

        // Initialize Poisson probabilities
        double[] lpoi = new double[L];
        double[] logTrace = new double[L];
        double[] poi = new double[L];
        for (int i = 0; i < L; i++) {
            lpoi[i] = -lambda * sortedTrace[i];
            logTrace[i] = Math.log(sortedTrace[i]);
            poi[i] = Math.exp(lpoi[i]);
        }
        double[] spoi = poi.clone();

        // fx = poi * (alpha * a)
        double alphaA = alpha.mult(a).get(0, 0);
        double[] fx = new double[L];
        for (int i = 0; i < L; i++) fx[i] = poi[i] * alphaA;

        int k = 1;
        int first = 0;
        Matrix coeffv = alpha.copy();
        int maxIter = 10000;

        while (first < L && k < maxIter) {
            coeffv = coeffv.mult(P);
            double coeffvA = coeffv.mult(a).get(0, 0);

            for (int i = first; i < L; i++) {
                lpoi[i] += Math.log(lambda) + logTrace[i] - Math.log(k);
                poi[i] = Math.exp(lpoi[i]);
                spoi[i] += poi[i];
                fx[i] += poi[i] * coeffvA;
            }

            k++;

            // Find new first index where spoi < 1 - eps
            int newFirst = L;
            for (int i = first; i < L; i++) {
                if (spoi[i] < 1 - eps) {
                    newFirst = i;
                    break;
                }
            }
            first = newFirst;
        }

        // Compute log-likelihood
        double logLikelihood = 0.0;
        for (int i = 0; i < L; i++) {
            if (fx[i] > 0) {
                logLikelihood += Math.log(fx[i]);
            }
        }

        return logLikelihood / L;
    }

    public static double likelihoodFromTracePH(double[] trace, Matrix alpha, Matrix A) {
        return likelihoodFromTracePH(trace, alpha, A, 1e-14);
    }

    /**
     * Evaluates the log-likelihood of a trace with the given MAP.
     * The result is divided by the length of the trace.
     */
    public static double likelihoodFromTraceMAP(double[] trace, Matrix D0, Matrix D1, double prec) {
        int N = D0.getNumRows();
        int L = trace.length;

        // Sort trace and keep original indices
        Integer[] idxArr = new Integer[L];
        for (int i = 0; i < L; i++) idxArr[i] = i;
        final double[] traceLocal = trace;
        Arrays.sort(idxArr, new java.util.Comparator<Integer>() {
            @Override
            public int compare(Integer a, Integer b) {
                return Double.compare(traceLocal[a], traceLocal[b]);
            }
        });
        double[] sortedTrace = new double[L];
        int[] ix = new int[L];
        for (int i = 0; i < L; i++) {
            sortedTrace[i] = trace[idxArr[i]];
            ix[i] = idxArr[i];
        }

        // Find lambda
        double lambda = 0.0;
        for (int i = 0; i < N; i++) {
            double absVal = Math.abs(D0.get(i, i));
            if (absVal > lambda) lambda = absVal;
        }

        // P = D0/lambda + I
        Matrix P = D0.scale(1.0 / lambda).add(Matrix.eye(N));

        double eps = Math.max(prec, Math.pow(10.0, Math.log10(prec) + Math.log10(lambda)));

        // Initialize
        double[] lpoi = new double[L];
        double[] logTrace = new double[L];
        double[] poi = new double[L];
        for (int i = 0; i < L; i++) {
            lpoi[i] = -lambda * sortedTrace[i];
            logTrace[i] = Math.log(sortedTrace[i]);
            poi[i] = Math.exp(lpoi[i]);
        }
        double[] spoi = poi.clone();

        // fx is L x N x N (stored as list of matrices)
        Matrix[] fx = new Matrix[L];
        for (int i = 0; i < L; i++) fx[i] = D1.scale(poi[i]);

        int k = 1;
        int first = 0;
        Matrix coeffv = D1.copy();
        int maxIter = 10000;

        while (first < L && k < maxIter) {
            coeffv = P.mult(coeffv);

            for (int i = first; i < L; i++) {
                lpoi[i] += Math.log(lambda) + logTrace[i] - Math.log(k);
                poi[i] = Math.exp(lpoi[i]);
                spoi[i] += poi[i];
                fx[i] = fx[i].add(coeffv.scale(poi[i]));
            }

            k++;

            int newFirst = L;
            for (int i = first; i < L; i++) {
                if (spoi[i] < 1 - eps) {
                    newFirst = i;
                    break;
                }
            }
            first = newFirst;
        }

        // Compute stationary distribution
        Matrix Pi = D0.neg().inv().mult(D1);
        Matrix alpha = DTMCSolve.dtmcSolve(Pi);

        // Compute log-likelihood by multiplying matrices in original order
        Matrix l = alpha.copy();
        int sc = 0;

        // Create reverse mapping
        int[] ixrev = new int[L];
        for (int i = 0; i < L; i++) {
            ixrev[ix[i]] = i;
        }

        for (int i = 0; i < L; i++) {
            l = l.mult(fx[ixrev[i]]);

            // Rescale periodically to avoid numerical issues
            if (i % 10 == 0) {
                double sumL = l.elementSum();
                if (sumL > 0) {
                    int scale = (int) Math.ceil(Math.log(sumL) / Math.log(2.0));
                    if (scale > 1) {
                        l = l.scale(1.0 / Math.pow(2.0, scale));
                        sc += scale;
                    }
                    if (scale < -10) {
                        int adjScale = scale + 10;
                        l = l.scale(1.0 / Math.pow(2.0, adjScale));
                        sc += adjScale;
                    }
                }
            }
        }

        return (Math.log(l.elementSum()) + sc * Math.log(2.0)) / L;
    }

    public static double likelihoodFromTraceMAP(double[] trace, Matrix D0, Matrix D1) {
        return likelihoodFromTraceMAP(trace, D0, D1, 1e-14);
    }
}
