/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map;

import java.util.ArrayList;
import java.util.List;

import jline.lib.butools.reptrans.FindMarkovianRepresentation;
import jline.lib.butools.reptrans.FindMarkovianRepresentation.EvalFun;
import jline.lib.butools.reptrans.FindMarkovianRepresentation.TransFun;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class MMAPFromMRAP {
    private MMAPFromMRAP() {}

    public static MatrixCell mmapFromMRAP(MatrixCell H) {
        return mmapFromMRAP(H, 1e-14);
    }

    /**
     * Obtains a Markovian representation of a continuous marked rational
     * arrival process of the same size, if possible.
     */
    public static MatrixCell mmapFromMRAP(MatrixCell H, double prec) {
        if (!CheckMRAPRepresentation.checkMRAPRepresentation(H, prec)) {
            throw new IllegalArgumentException("MMAPFromMRAP: Input isn't a valid MRAP representation!");
        }

        // Convert MatrixCell to List<Matrix>
        List<Matrix> HList = new ArrayList<Matrix>();
        for (int i = 0; i < H.size(); i++) {
            HList.add(H.get(i));
        }

        // Transformation function: nH = inv(B) * oH * B for each matrix
        TransFun transfun = new TransFun() {
            @Override
            public List<Matrix> apply(List<Matrix> oH, Matrix B) {
                Matrix Binv = B.inv();
                List<Matrix> out = new ArrayList<Matrix>();
                for (Matrix m : oH) {
                    out.add(Binv.mult(m).mult(B));
                }
                return out;
            }
        };

        // Evaluation function for continuous MAP
        EvalFun evalfun = new EvalFun() {
            @Override
            public double apply(List<Matrix> oH, int k) {
                if (k % 2 == 0) {
                    double dist = Double.POSITIVE_INFINITY;
                    for (int i = 0; i < oH.size(); i++) {
                        if (i == 0) {
                            int n = oH.get(0).getNumRows();
                            for (int r = 0; r < n; r++) {
                                for (int c = 0; c < n; c++) {
                                    if (r != c) {
                                        dist = Math.min(dist, oH.get(0).get(r, c));
                                    }
                                }
                            }
                        } else {
                            dist = Math.min(dist, oH.get(i).elementMin());
                        }
                    }
                    return -dist;
                } else {
                    double dist = 0.0;
                    for (int i = 0; i < oH.size(); i++) {
                        if (i == 0) {
                            int n = oH.get(0).getNumRows();
                            for (int r = 0; r < n; r++) {
                                for (int c = 0; c < n; c++) {
                                    if (r != c) {
                                        double v = oH.get(0).get(r, c);
                                        if (v < 0) dist += v;
                                    }
                                }
                            }
                        } else {
                            for (int r = 0; r < oH.get(i).getNumRows(); r++) {
                                for (int c = 0; c < oH.get(i).getNumCols(); c++) {
                                    double v = oH.get(i).get(r, c);
                                    if (v < 0) dist += v;
                                }
                            }
                        }
                    }
                    return -dist;
                }
            }
        };

        List<Matrix> result = FindMarkovianRepresentation.findMarkovianRepresentation(HList, transfun, evalfun, prec);

        // Convert back to MatrixCell
        MatrixCell resultCell = new MatrixCell(result.size());
        for (int i = 0; i < result.size(); i++) {
            resultCell.set(i, result.get(i));
        }
        return resultCell;
    }

    /**
     * Overload for Matrix[].
     */
    public static MatrixCell mmapFromMRAP(Matrix[] H, double prec) {
        MatrixCell cell = new MatrixCell(H.length);
        for (int i = 0; i < H.length; i++) {
            cell.set(i, H[i]);
        }
        return mmapFromMRAP(cell, prec);
    }

    public static MatrixCell mmapFromMRAP(Matrix[] H) {
        return mmapFromMRAP(H, 1e-14);
    }
}
