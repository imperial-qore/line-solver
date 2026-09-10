/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph;

import java.util.ArrayList;
import java.util.List;

import java.util.function.BiFunction;

import jline.lib.butools.reptrans.FindMarkovianRepresentation;
import jline.util.matrix.Matrix;

public final class PHFromME {
    private PHFromME() {}

    /**
     * Obtains a Markovian representation of a matrix exponential
     * distribution of the same size, if possible, using elementary
     * similarity transformations.
     */
    public static PHRepresentation phFromME(Matrix alpha, Matrix A, double prec) {
        if (!CheckMERepresentation.checkMERepresentation(alpha, A, prec)) {
            throw new IllegalArgumentException("PHFromME: Input is not a valid ME representation!");
        }

        final List<Matrix> repList = new ArrayList<Matrix>();
        repList.add(alpha);
        repList.add(A);

        BiFunction<List<Matrix>, Matrix, List<Matrix>> transfun = new BiFunction<List<Matrix>, Matrix, List<Matrix>>() {
            @Override
            public List<Matrix> apply(List<Matrix> oH, Matrix B) {
                Matrix Binv = B.inv();
                List<Matrix> out = new ArrayList<Matrix>(2);
                out.add(oH.get(0).mult(B));
                out.add(Binv.mult(oH.get(1)).mult(B));
                return out;
            }
        };

        BiFunction<List<Matrix>, Integer, Double> evalfun = new BiFunction<List<Matrix>, Integer, Double>() {
            @Override
            public Double apply(List<Matrix> oH, Integer k) {
                Matrix ao = oH.get(0);
                Matrix Ao = oH.get(1);
                int M = Ao.getNumRows();

                double[] av = new double[M];
                for (int i = 0; i < M; i++) {
                    double rowSum = 0.0;
                    for (int j = 0; j < M; j++) {
                        rowSum += Ao.get(i, j);
                    }
                    av[i] = -rowSum;
                }

                if (k % 2 == 0) {
                    double minVal = Double.MAX_VALUE;
                    for (int j = 0; j < ao.getNumCols(); j++) {
                        double v = ao.get(0, j);
                        if (v < minVal) minVal = v;
                    }
                    for (int i = 0; i < M; i++) {
                        if (av[i] < minVal) minVal = av[i];
                    }
                    for (int i = 0; i < M; i++) {
                        for (int j = 0; j < M; j++) {
                            if (i != j) {
                                double v = Ao.get(i, j);
                                if (v < minVal) minVal = v;
                            }
                        }
                    }
                    return -minVal;
                } else {
                    double sumNeg = 0.0;
                    for (int j = 0; j < ao.getNumCols(); j++) {
                        double v = ao.get(0, j);
                        if (v < 0.0) sumNeg += v;
                    }
                    for (int i = 0; i < M; i++) {
                        if (av[i] < 0.0) sumNeg += av[i];
                    }
                    for (int i = 0; i < M; i++) {
                        for (int j = 0; j < M; j++) {
                            if (i != j) {
                                double v = Ao.get(i, j);
                                if (v < 0.0) sumNeg += v;
                            }
                        }
                    }
                    return -sumNeg;
                }
            }
        };

        List<Matrix> result = FindMarkovianRepresentation.findMarkovianRepresentation(repList, transfun, evalfun, prec);
        return new PHRepresentation(result.get(0), result.get(1));
    }

    public static PHRepresentation phFromME(Matrix alpha, Matrix A) {
        return phFromME(alpha, A, 1e-14);
    }
}
