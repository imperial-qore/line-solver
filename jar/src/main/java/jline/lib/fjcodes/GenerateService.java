/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.lib.fjcodes;

import jline.util.matrix.Matrix;

public final class GenerateService {
    private GenerateService() {}

    /**
     * Generate Phase-Type representation for service time.
     */
    public static GenerateServiceResult generateService(FJService services, FJServiceH service_h, int C, Matrix S) {
        int dim = service_h.getBeta().length();
        int m = services.getTau_st().length();

        Matrix indexes_notbusy = FjCodesUtilsKt.build_index(m, 1);
        int dim_NB = indexes_notbusy.getNumRows();

        int dim_C = C + 1;

        int newdim = dim_C * dim;
        int dim_notbusy = dim_C * dim_NB;

        Matrix T = new Matrix(newdim + dim_notbusy, newdim + dim_notbusy);
        Matrix t = new Matrix(newdim + dim_notbusy, 1);

        for (int i = 0; i < dim_NB; i++) {
            t.set(newdim + dim_notbusy - dim_NB + i, 0, services.getSt().get(i, 0));
        }

        FjCodesUtilsKt.setSubMatrix(T, 0, 0, S);

        Matrix S_long = new Matrix(dim, dim_NB);

        for (int row = 0; row < dim; row++) {
            double[] countvect = FjCodesUtilsKt.getRowAsArray(service_h.getService_phases(), row);

            for (int i = m; i < 2 * m; i++) {
                if (countvect[i] > 0.0) {
                    double[] tovect = new double[m];
                    for (int k = 0; k < m; k++) {
                        tovect[k] = countvect[k];
                    }

                    int col = FjCodesUtilsKt.vectmatch(tovect, indexes_notbusy) - 1;

                    if (col >= 0) {
                        double currentVal = S_long.get(row, col);
                        S_long.set(row, col, currentVal + countvect[i] * services.getSt().get(i - m, 0));
                    }
                }
            }
        }

        for (int row = 0; row < dim_C - 1; row++) {
            FjCodesUtilsKt.setSubMatrix(T, row * dim, newdim + row * dim_NB, S_long);
        }

        Matrix S_last = new Matrix(dim, dim_NB);

        for (int row = 0; row < dim; row++) {
            double[] countvect = FjCodesUtilsKt.getRowAsArray(service_h.getService_phases(), row);

            for (int k = 0; k <= 1; k++) {
                int startIdx = k * m;
                int endIdx = (k + 1) * m;

                for (int i = startIdx; i < endIdx; i++) {
                    if (countvect[i] > 0.0) {
                        int otherStart = (1 - k) * m;
                        double[] tovect = new double[m];
                        for (int j = 0; j < m; j++) {
                            tovect[j] = countvect[otherStart + j];
                        }

                        int col = FjCodesUtilsKt.vectmatch(tovect, indexes_notbusy) - 1;

                        if (col >= 0) {
                            double currentVal = S_last.get(row, col);
                            S_last.set(row, col, currentVal + countvect[i] * services.getSt().get(i - startIdx, 0));
                        }
                    }
                }
            }
        }

        FjCodesUtilsKt.setSubMatrix(T, (dim_C - 1) * dim, newdim + (dim_C - 1) * dim_NB, S_last);

        for (int row = 0; row < dim_C; row++) {
            FjCodesUtilsKt.setSubMatrix(T, newdim + row * dim_NB, newdim + row * dim_NB, services.getST());
        }

        Matrix A = new Matrix(m, m);
        for (int i = 0; i < m; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < m; j++) {
                rowSum += services.getST().get(i, j);
            }
            for (int j = 0; j < m; j++) {
                A.set(i, j, -rowSum * services.getTau_st().get(j, 0));
            }
        }

        for (int row = 0; row < dim_C - 1; row++) {
            FjCodesUtilsKt.setSubMatrix(T, newdim + row * dim_NB, newdim + (row + 1) * dim_NB, A);
        }

        for (int row = newdim; row < newdim + dim_notbusy; row++) {
            T.set(row, row, 0.0);
            double rowSum = 0.0;
            for (int col = 0; col < newdim + dim_notbusy; col++) {
                rowSum += T.get(row, col);
            }
            rowSum += t.get(row, 0);
            T.set(row, row, -rowSum);
        }

        return new GenerateServiceResult(T, newdim, dim_notbusy);
    }
}
