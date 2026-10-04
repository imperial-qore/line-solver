/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import jline.lib.butools.dph.MGFromMoments.MGRepresentation;
import jline.lib.butools.reptrans.FindMarkovianRepresentation;
import jline.lib.butools.reptrans.FindMarkovianRepresentation.EvalFun;
import jline.lib.butools.reptrans.FindMarkovianRepresentation.TransFun;
import jline.util.matrix.Matrix;

public final class DPHFromMG {
    private DPHFromMG() {}

    public static MGRepresentation dphFromMG(Matrix alpha, Matrix A) {
        return dphFromMG(alpha, A, 1e-14);
    }

    /**
     * Obtains a Markovian representation of a matrix-geometric distribution
     * of the same size, if possible.
     */
    public static MGRepresentation dphFromMG(Matrix alpha, Matrix A, double precision) {
        if (!CheckMGRepresentation.checkMGRepresentation(alpha, A)) {
            throw new IllegalArgumentException("DPHFromMG: Input isn't a valid MG distribution!");
        }

        // Transform function: apply similarity transformation B to representation
        TransFun transfun = new TransFun() {
            @Override
            public List<Matrix> apply(List<Matrix> rep, Matrix B) {
                Matrix newAlpha = rep.get(0).mult(B);
                Matrix newA = B.inv().mult(rep.get(1)).mult(B);
                return Arrays.asList(newAlpha, newA);
            }
        };

        // Evaluation function: measure distance from Markovian representation
        EvalFun evalfun = new EvalFun() {
            @Override
            public double apply(List<Matrix> rep, int k) {
                Matrix ao = rep.get(0);
                Matrix Ao = rep.get(1);
                int N = Ao.getNumRows();

                // av = 1 - sum(Ao, 2) (closing vector)
                double[] av = new double[N];
                for (int i = 0; i < N; i++) {
                    double rowSum = 0.0;
                    for (int j = 0; j < N; j++) {
                        rowSum += Ao.get(i, j);
                    }
                    av[i] = 1.0 - rowSum;
                }

                // Ad = Ao - diag(diag(Ao)) (off-diagonal elements)
                Matrix Ad = Matrix.zeros(N, N);
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        if (i != j) {
                            Ad.set(i, j, Ao.get(i, j));
                        }
                    }
                }

                if (k % 2 == 0) {
                    // Return max negative value
                    double minVal = Double.MAX_VALUE;
                    for (int i = 0; i < ao.length(); i++) {
                        if (ao.get(i) < minVal) minVal = ao.get(i);
                    }
                    for (int i = 0; i < N; i++) {
                        if (av[i] < minVal) minVal = av[i];
                    }
                    for (int i = 0; i < N; i++) {
                        for (int j = 0; j < N; j++) {
                            if (Ad.get(i, j) < minVal) minVal = Ad.get(i, j);
                        }
                    }
                    return -Math.min(0.0, minVal);
                } else {
                    // Return sum of negative values
                    double negSum = 0.0;
                    for (int i = 0; i < ao.length(); i++) {
                        if (ao.get(i) < 0) negSum -= ao.get(i);
                    }
                    for (int i = 0; i < N; i++) {
                        if (av[i] < 0) negSum -= av[i];
                    }
                    for (int i = 0; i < N; i++) {
                        for (int j = 0; j < N; j++) {
                            if (Ad.get(i, j) < 0) negSum -= Ad.get(i, j);
                        }
                    }
                    return negSum;
                }
            }
        };

        List<Matrix> rep = new ArrayList<Matrix>();
        rep.add(alpha);
        rep.add(A);
        List<Matrix> nrep = FindMarkovianRepresentation.findMarkovianRepresentation(rep, transfun, evalfun, precision);

        return new MGRepresentation(nrep.get(0), nrep.get(1));
    }

    /**
     * Overload for double[] alpha.
     */
    public static MGRepresentation dphFromMG(double[] alpha, Matrix A, double precision) {
        return dphFromMG(new Matrix(alpha), A, precision);
    }

    public static MGRepresentation dphFromMG(double[] alpha, Matrix A) {
        return dphFromMG(new Matrix(alpha), A, 1e-14);
    }
}
