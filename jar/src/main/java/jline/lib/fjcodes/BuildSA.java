/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.lib.fjcodes;

import jline.util.matrix.Matrix;

public final class BuildSA {
    private BuildSA() {}

    /**
     * Build state transition matrices S and A_jump.
     */
    public static SAResult build_SA(FJService services, FJServiceH service_h, int C) {
        int dim = service_h.getBeta().length();
        int m = services.getTau_st().length();
        int dim_C = C + 1;
        int newdim = dim_C * dim;

        Matrix S = new Matrix(newdim, newdim);
        Matrix A_jump = new Matrix(newdim, newdim);

        for (int row = 0; row < dim_C; row++) {
            FjCodesUtils.setSubMatrix(S, row * dim, row * dim, service_h.getS());
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

        Matrix S_Cminus1 = new Matrix(dim, dim);

        for (int row = 0; row < dim; row++) {
            double[] countvect = FjCodesUtils.getRowAsArray(service_h.getService_phases(), row);

            for (int i = 0; i < m; i++) {
                if (countvect[i] > 0.0) {
                    for (int j = 0; j < m; j++) {
                        double[] tovect = countvect.clone();
                        tovect[i] = tovect[i] - 1.0;
                        tovect[j] = tovect[j] + 1.0;

                        int col = FjCodesUtils.vectmatch(tovect, service_h.getService_phases()) - 1;

                        if (col >= 0) {
                            double currentVal = S_Cminus1.get(row, col);
                            S_Cminus1.set(row, col, currentVal + countvect[i] * A.get(i, j));
                        }
                    }
                }
            }
        }

        for (int c = C; c >= 1; c--) {
            int sourceRow = (C - c) * dim;
            int targetCol = (C - c + 1) * dim;
            FjCodesUtils.setSubMatrix(S, sourceRow, targetCol, S_Cminus1);
        }

        Matrix A_Cplus1 = new Matrix(dim, dim);

        for (int row = 0; row < dim; row++) {
            double[] countvect = FjCodesUtils.getRowAsArray(service_h.getService_phases(), row);

            for (int i = m; i < 2 * m; i++) {
                if (countvect[i] > 0.0) {
                    for (int j = m; j < 2 * m; j++) {
                        double[] tovect = countvect.clone();
                        tovect[i] = tovect[i] - 1.0;
                        tovect[j] = tovect[j] + 1.0;

                        int col = FjCodesUtils.vectmatch(tovect, service_h.getService_phases()) - 1;

                        if (col >= 0) {
                            double currentVal = A_Cplus1.get(row, col);
                            A_Cplus1.set(row, col, currentVal + A.get(i - m, j - m));
                        }
                    }
                }
            }
        }

        for (int c = C - 1; c >= 1; c--) {
            int sourceRow = (C - c) * dim;
            int targetCol = (C - c - 1) * dim;
            FjCodesUtils.setSubMatrix(A_jump, sourceRow, targetCol, A_Cplus1);
        }

        FjCodesUtils.setSubMatrix(A_jump, 0, 0, A_Cplus1);

        Matrix A_last = new Matrix(dim, dim);

        for (int row = 0; row < dim; row++) {
            double[] countvect = FjCodesUtils.getRowAsArray(service_h.getService_phases(), row);

            for (int k = 0; k <= 1; k++) {
                int startIdx = k * m;
                int endIdx = (k + 1) * m;

                for (int i = startIdx; i < endIdx; i++) {
                    if (countvect[i] > 0.0) {
                        int otherStart = (1 - k) * m;
                        double[] tovect = new double[2 * m];
                        for (int idx = 0; idx < m; idx++) {
                            tovect[idx] = countvect[otherStart + idx];
                        }

                        for (int j = 0; j < m; j++) {
                            double[] temp_tovector = tovect.clone();
                            temp_tovector[m + j] = 1.0;

                            int col = FjCodesUtils.vectmatch(temp_tovector, service_h.getService_phases()) - 1;

                            if (col >= 0) {
                                double currentVal = A_last.get(row, col);
                                A_last.set(row, col, currentVal + countvect[i] * A.get(i - startIdx, j));
                            }
                        }
                    }
                }
            }
        }

        FjCodesUtils.setSubMatrix(A_jump, C * dim, (C - 1) * dim, A_last);

        return new SAResult(S, A_jump);
    }
}
