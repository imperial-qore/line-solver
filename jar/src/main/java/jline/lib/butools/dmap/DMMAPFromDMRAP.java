/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import java.util.ArrayList;
import java.util.List;

import java.util.function.BiFunction;

import jline.lib.butools.reptrans.FindMarkovianRepresentation;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class DMMAPFromDMRAP {
    private DMMAPFromDMRAP() {}

    /**
     * Obtains a Markovian representation of a discrete rational arrival process
     * of the same size, if possible.
     *
     * @param H    The H0...HN matrices of the DMRAP to transform (as MatrixCell)
     * @param prec Precision threshold
     * @return The D0...DN matrices of the DMMAP (if found)
     */
    public static MatrixCell dmmapFromDMRAP(MatrixCell H, double prec) {
        if (!CheckDMRAPRepresentation.checkDMRAPRepresentation(H, prec)) {
            throw new IllegalArgumentException("DMMAPFromDMRAP: Input isn't a valid DMRAP representation!");
        }

        final List<Matrix> HList = new ArrayList<Matrix>();
        for (int i = 0; i < H.size(); i++) {
            HList.add(H.get(i));
        }

        // Transformation function: nH = inv(B) * oH * B for each matrix
        BiFunction<List<Matrix>, Matrix, List<Matrix>> transfun = new BiFunction<List<Matrix>, Matrix, List<Matrix>>() {
            @Override
            public List<Matrix> apply(List<Matrix> oH, Matrix B) {
                Matrix Binv = B.inv();
                List<Matrix> out = new ArrayList<Matrix>(oH.size());
                for (Matrix m : oH) {
                    out.add(Binv.mult(m).mult(B));
                }
                return out;
            }
        };

        BiFunction<List<Matrix>, Integer, Double> evalfun = new BiFunction<List<Matrix>, Integer, Double>() {
            @Override
            public Double apply(List<Matrix> oH, Integer k) {
                Matrix ones = Matrix.ones(oH.get(0).getNumRows(), oH.get(0).getNumCols());
                if (k % 2 == 0) {
                    double dist = Double.POSITIVE_INFINITY;
                    for (int i = 0; i < oH.size(); i++) {
                        double minOH = oH.get(i).elementMin();
                        double minOnesMinusOH = ones.sub(oH.get(i)).elementMin();
                        dist = Math.min(dist, Math.min(minOH, minOnesMinusOH));
                    }
                    return -dist;
                } else {
                    double dist = 0.0;
                    for (int i = 0; i < oH.size(); i++) {
                        double sumNegOH = 0.0;
                        double sumNegOnesMinusOH = 0.0;
                        Matrix oHi = oH.get(i);
                        for (int r = 0; r < oHi.getNumRows(); r++) {
                            for (int c = 0; c < oHi.getNumCols(); c++) {
                                double v = oHi.get(r, c);
                                if (v < 0) sumNegOH += v;
                                double vOnes = 1.0 - v;
                                if (vOnes < 0) sumNegOnesMinusOH += vOnes;
                            }
                        }
                        dist += Math.min(sumNegOH, sumNegOnesMinusOH);
                    }
                    return -dist;
                }
            }
        };

        List<Matrix> result = FindMarkovianRepresentation.findMarkovianRepresentation(HList, transfun, evalfun, prec);

        MatrixCell resultCell = new MatrixCell(result.size());
        for (int i = 0; i < result.size(); i++) {
            resultCell.set(i, result.get(i));
        }
        return resultCell;
    }

    public static MatrixCell dmmapFromDMRAP(MatrixCell H) {
        return dmmapFromDMRAP(H, 1e-14);
    }

    public static MatrixCell dmmapFromDMRAP(Matrix[] H, double prec) {
        MatrixCell cell = new MatrixCell(H.length);
        for (int i = 0; i < H.length; i++) {
            cell.set(i, H[i]);
        }
        return dmmapFromDMRAP(cell, prec);
    }

    public static MatrixCell dmmapFromDMRAP(Matrix[] H) {
        return dmmapFromDMRAP(H, 1e-14);
    }
}
