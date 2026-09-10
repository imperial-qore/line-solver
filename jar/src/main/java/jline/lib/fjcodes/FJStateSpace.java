/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.lib.fjcodes;

import jline.util.matrix.Matrix;

public final class FJStateSpace {
    private FJStateSpace() {}

    /**
     * Construct state space matrix for not-all-busy states.
     *
     * @param C Capacity parameter
     * @param services Service process
     * @param service_h 2-node service representation
     * @return State transition matrix for not-all-busy states
     */
    public static Matrix constructNotAllBusy(int C, FJService services, FJServiceH service_h) {
        int dim = service_h.getBeta().length();
        int m = services.getTau_st().length();

        Matrix indexes_notbusy = FJUtils.build_index(m, 1);
        int dim_NB = indexes_notbusy.getNumRows();
        int dim_C = C + 1;
        int dim_notbusy = (dim_C - 1) * dim_NB + 1;

        Matrix S_notallbusy = new Matrix(dim_notbusy, dim_notbusy);

        // From not-busy to not-busy
        for (int row = 0; row < dim_C - 1; row++) {
            int startRow = row * dim_NB;
            FJUtils.setSubMatrix(S_notallbusy, startRow, startRow, services.getST());
        }

        // A = -sum(ST, 2) * tau_st  (row sum times initial prob)
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

        // Transitions between levels
        for (int row = 0; row < dim_C - 2; row++) {
            int startRow = row * dim_NB;
            int startCol = (row + 1) * dim_NB;
            FJUtils.setSubMatrix(S_notallbusy, startRow, startCol, A);
        }

        // Last transition to absorbing state
        int startRow = (dim_C - 2) * dim_NB;
        int endRow = (dim_C - 1) * dim_NB;
        int startCol = (dim_C - 1) * dim_NB;
        for (int i = startRow; i < endRow; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < m; j++) {
                rowSum += services.getST().get(i - startRow, j);
            }
            S_notallbusy.set(i, startCol, -rowSum);
        }

        // Normalize diagonal
        for (int row = 0; row < dim_notbusy; row++) {
            double rowSum = 0.0;
            for (int col = 0; col < dim_notbusy; col++) {
                if (row != col) {
                    rowSum += S_notallbusy.get(row, col);
                }
            }
            S_notallbusy.set(row, row, -rowSum);
        }

        return S_notallbusy;
    }

    /**
     * Construct extended state space matrices (Se, Sestar, R0, Ke, Kc).
     *
     * @param C Capacity parameter
     * @param services Service process
     * @param service_h 2-node service representation
     * @param S Busy state transition matrix
     * @return SRKResult with all constructed matrices
     */
    public static SRKResult constructSRK(int C, FJService services, FJServiceH service_h, Matrix S) {
        int dim = service_h.getBeta().length();
        int m = services.getTau_st().length();

        Matrix indexes_notbusy = FJUtils.build_index(m, 1);
        int dim_NB = indexes_notbusy.getNumRows();
        int dim_C = C + 1;
        int newdim = dim_C * dim;
        int dim_notbusy = (dim_C - 1) * dim_NB + 1;

        Matrix Se = new Matrix(newdim + dim_notbusy, newdim + dim_notbusy);
        Matrix Sestar = new Matrix(newdim + dim_notbusy, newdim + dim_notbusy);

        // Copy S into upper-left block
        FJUtils.setSubMatrix(Se, 0, 0, S);

        // From busy to not-busy: when the job in the shorter queue completes service
        Matrix S_long = new Matrix(dim, dim_NB);

        for (int row = 0; row < dim; row++) {
            double[] countvect = FJUtils.getRowAsArray(service_h.getService_phases(), row);

            for (int i = m; i < 2 * m; i++) {  // Starting phase
                if (countvect[i] > 0) {
                    // Extract tovect (first m elements)
                    double[] tovect = new double[m];
                    for (int k = 0; k < m; k++) {
                        tovect[k] = countvect[k];
                    }

                    int col = FJUtils.vectmatch(tovect, indexes_notbusy) - 1;

                    if (col >= 0) {
                        double currentVal = S_long.get(row, col);
                        S_long.set(row, col, currentVal + countvect[i] * services.getSt().get(i - m, 0));
                    }
                }
            }
        }

        // Set S_long blocks in Se
        FJUtils.setSubMatrix(Se, 0, newdim, S_long);
        for (int row = 1; row < dim_C - 1; row++) {
            FJUtils.setSubMatrix(Se, row * dim, newdim + (row - 1) * dim_NB, S_long);
        }

        // When the 2 queues have the same length
        Matrix S_last = new Matrix(dim, dim_NB);

        for (int row = 0; row < dim; row++) {
            double[] countvect = FJUtils.getRowAsArray(service_h.getService_phases(), row);

            for (int k = 0; k <= 1; k++) {
                int startIdx = k * m;
                int endIdx = (k + 1) * m;

                for (int i = startIdx; i < endIdx; i++) {
                    if (countvect[i] > 0) {
                        // Extract tovect from the other queue
                        int otherStart = (1 - k) * m;
                        double[] tovect = new double[m];
                        for (int j = 0; j < m; j++) {
                            tovect[j] = countvect[otherStart + j];
                        }

                        int col = FJUtils.vectmatch(tovect, indexes_notbusy) - 1;

                        if (col >= 0) {
                            double currentVal = S_last.get(row, col);
                            S_last.set(row, col, currentVal + countvect[i] * services.getSt().get(i - startIdx, 0));
                        }
                    }
                }
            }
        }

        FJUtils.setSubMatrix(Se, (dim_C - 1) * dim, newdim + (dim_C - 2) * dim_NB, S_last);

        // Sestar gets transitions from busy to not-busy
        FJUtils.setSubMatrix(Sestar, 0, newdim,
                Matrix.getSubMatrix(Se, 0, newdim, newdim, newdim + dim_notbusy));

        // From not-busy to not-busy (same as constructNotAllBusy logic)
        for (int row = 0; row < dim_C - 1; row++) {
            int start = newdim + row * dim_NB;
            FJUtils.setSubMatrix(Se, start, start, services.getST());
        }

        // A = -sum(ST, 2) * tau_st
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

        for (int row = 0; row < dim_C - 2; row++) {
            int startRow = newdim + row * dim_NB;
            int startCol = newdim + (row + 1) * dim_NB;
            FJUtils.setSubMatrix(Se, startRow, startCol, A);
        }

        // Last block
        int lastBlockRow = newdim + (dim_C - 2) * dim_NB;
        int lastBlockCol = newdim + (dim_C - 1) * dim_NB;
        for (int i = 0; i < dim_NB; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < m; j++) {
                rowSum += services.getST().get(i, j);
            }
            Se.set(lastBlockRow + i, lastBlockCol, -rowSum);
        }

        // Normalize diagonal for not-busy part
        for (int row = newdim; row < newdim + dim_notbusy; row++) {
            double rowSum = 0.0;
            for (int col = 0; col < newdim + dim_notbusy; col++) {
                if (row != col) {
                    rowSum += Se.get(row, col);
                }
            }
            Se.set(row, row, -rowSum);
        }

        // Construct R0 matrix: from not-busy to busy
        Matrix R0 = new Matrix(newdim + dim_notbusy, newdim + dim_notbusy);
        Matrix R_NB = new Matrix(dim_NB, dim);

        for (int row = 0; row < dim_NB; row++) {
            double[] countvect = FJUtils.getRowAsArray(indexes_notbusy, row);

            for (int i = 0; i < m; i++) {
                // tovect_short has a 1 at position i
                double[] tovect = new double[2 * m];
                for (int j = 0; j < m; j++) {
                    tovect[j] = countvect[j];
                }
                tovect[m + i] = 1.0;

                int col = FJUtils.vectmatch(tovect, service_h.getService_phases()) - 1;

                if (col >= 0) {
                    double currentVal = R_NB.get(row, col);
                    R_NB.set(row, col, currentVal + services.getTau_st().get(i, 0));
                }
            }
        }

        for (int row = 0; row < dim_C - 1; row++) {
            FJUtils.setSubMatrix(R0, newdim + row * dim_NB, row * dim, R_NB);
        }

        // Last row of R0
        for (int i = 0; i < dim; i++) {
            R0.set(newdim + dim_notbusy - 1, newdim - dim + i, service_h.getBeta().get(0, i));
        }

        // Construct Ke and Kc
        Matrix Ke = new Matrix(newdim, newdim + dim_notbusy);
        for (int i = 0; i < newdim; i++) {
            Ke.set(i, i, 1.0);
        }

        Matrix Kc = new Matrix(newdim + dim_notbusy, newdim);
        for (int i = 0; i < newdim; i++) {
            Kc.set(i, i, 1.0);
        }

        return new SRKResult(Se, Sestar, R0, Ke, Kc, S_long, S_last);
    }
}
