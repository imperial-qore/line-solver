/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.lib.fjcodes;

import jline.util.matrix.Matrix;

public final class ReturnRT2 {
    private ReturnRT2() {}

    /**
     * Compute response time percentiles for K=2 Fork-Join queue.
     *
     * @param arrival Arrival process
     * @param service Service process for single subtask
     * @param pers Array of percentile levels to compute
     * @param C Capacity parameter
     * @param tMode T-matrix computation method ("NARE" or "Sylvest")
     * @return Matrix with percentile results
     */
    public static Matrix returnRT2(FJArrival arrival, FJService service, double[] pers, int C, String tMode) {
        FJServiceH service_h = BuildServiceH.build_Service_h(service);

        ComputeT.ComputeTResult computeTResult = ComputeT.computeT(arrival, service, service_h, C, tMode);
        Matrix T = computeTResult.T;
        Matrix S = computeTResult.S;
        Matrix A_jump = computeTResult.A_jump;
        Matrix S_Arr = computeTResult.S_Arr;
        Matrix sum_Ajump = computeTResult.sum_Ajump;



        int mWait = A_jump.getNumRows();

        // phi = sum(T - S_Arr, 2)
        Matrix phi = new Matrix(T.getNumRows(), 1);
        for (int i = 0; i < T.getNumRows(); i++) {
            double sum = 0.0;
            for (int j = 0; j < T.getNumCols(); j++) {
                sum += T.get(i, j) - S_Arr.get(i, j);
            }
            phi.set(i, 0, sum);
        }

        PiResult piResult = ComputePi.computePi(T, arrival, service, service_h, C, S, A_jump);
        Matrix pi0 = piResult.pi0;
        double En1 = piResult.En1;

        ReturnWaitResult waitResult = ReturnWait.returnWait(En1, pi0, T, phi, sum_Ajump);
        Matrix wait_alpha = waitResult.wait_alpha;
        Matrix wait_Smat = waitResult.wait_Smat;
        double prob_wait = waitResult.prob_wait;
        Matrix alfa = waitResult.alfa;

        GenerateServiceResult serviceResult = GenerateService.generateService(service, service_h, C, S);
        Matrix ST = serviceResult.getT();
        int dim = serviceResult.getNewdim();
        int dim_notbusy = serviceResult.getDim_notbusy();

        ST = ST.kron(Matrix.eye(arrival.getLambda0().getNumRows()));

        // Normalize pi0
        double pi0Sum = elementSumLocal(pi0);
        pi0 = pi0.scale(1.0 / pi0Sum);

        int ma = arrival.getLambda0().getNumRows();
        int dim_ma = ma * dim;
        int dim_notbusy_ma = ma * dim_notbusy;
        int dim_service = dim_ma + dim_notbusy_ma;

        // Starting state of service for a job in not-all-busy period
        Matrix notbusy_start = new Matrix(1, dim_service);
        for (int i = 0; i < dim_ma; i++) {
            notbusy_start.set(0, i, (1.0 - prob_wait) * pi0.get(0, i));
        }

        // Starting state of service for a job in all-busy
        Matrix busy_start = new Matrix(1, dim_service);

        // TS = T - S_Arr
        Matrix TS = T.add(-1.0, S_Arr);

        Matrix alfaTS = alfa.mult(TS);
        double sumAlfaTS = elementSumLocal(alfaTS);

        for (int i = 0; i < dim_ma; i++) {
            busy_start.set(0, i, prob_wait * alfaTS.get(0, i) / sumAlfaTS);
        }


        // PH representation of response time
        int Tr = ST.getNumRows();
        int Sc = wait_Smat.getNumCols();

        // TS_sum
        Matrix TS_sum = new Matrix(TS.getNumRows(), ma * mWait);
        for (int i = 0; i < TS.getNumRows(); i++) {
            double rowSum = 0.0;
            for (int j = 0; j < TS.getNumCols(); j++) {
                rowSum += TS.get(i, j);
            }
            for (int j = 0; j < ma * mWait; j++) {
                TS_sum.set(i, j, rowSum);
            }
        }

        // TS_normalized
        Matrix TS_normalized = new Matrix(TS.getNumRows(), TS.getNumCols());
        for (int i = 0; i < TS.getNumRows(); i++) {
            for (int j = 0; j < TS.getNumCols(); j++) {
                if (TS_sum.get(i, 0) > 1e-12) {
                    TS_normalized.set(i, j, TS.get(i, j) / TS_sum.get(i, 0));
                }
            }
        }

        // stat_service_phase = busy_start / (-ST)
        Matrix negST = ST.scale(-1.0);
        Matrix stat_service_phase = busy_start.rightMatrixDivide(negST);


        // Select the service phases with meaningful stationary mass. MATLAB uses
        // busy_nz = stat_service_phase > 0, but the sub-generator balance of C_res
        // is sensitive to phases whose stat is only floating-point noise around 0
        // (~1e-17): their sign differs between MATLAB's LAPACK and this linear
        // algebra, changing which phases are kept and corrupting the balance. Use a
        // relative threshold so only phases with genuine mass are kept, which is
        // robust to that noise and keeps C_res close to a sub-generator.
        double statMax = 0.0;
        for (int i = 0; i < stat_service_phase.length(); i++) {
            statMax = Math.max(statMax, stat_service_phase.get(0, i));
        }
        double statTol = 1e-10 * statMax;

        boolean[] busy_nz = new boolean[Tr];
        int m_tr_ST = 0;

        for (int i = 0; i < stat_service_phase.length(); i++) {
            double val_i = stat_service_phase.get(0, i);
            if (val_i > statTol) {
                busy_nz[i] = true;
                m_tr_ST++;
            }
        }

        // tr_start_state_full
        Matrix tr_start_state_full = new Matrix(1, Tr);
        for (int i = 0; i < Tr; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < ST.getNumCols(); j++) {
                rowSum += ST.get(i, j);
            }
            tr_start_state_full.set(0, i, -rowSum * stat_service_phase.get(0, i));
        }

        // Build tr_ST (time-reversed service time matrix). MATLAB restricts to
        // tr_ST(busy_nz, busy_nz): both the row AND the column must be a busy_nz
        // phase. The earlier port only guarded the row, leaving entries in
        // non-busy_nz columns that are dropped in the reduction below, which broke
        // the row-sum balance of C_res (making it a non-sub-generator).
        Matrix tr_ST_full = new Matrix(Tr, Tr);
        for (int i = 0; i < Tr; i++) {
            for (int j = 0; j < Tr; j++) {
                if (busy_nz[i] && busy_nz[j]) {
                    tr_ST_full.set(i, j,
                            ST.get(j, i) * stat_service_phase.get(0, j) / stat_service_phase.get(0, i));
                }
            }
        }

        // tr_ST_exit = sum(-tr_ST, 2)
        Matrix tr_ST_exit_full = new Matrix(Tr, 1);
        for (int i = 0; i < Tr; i++) {
            double sum = 0.0;
            for (int j = 0; j < Tr; j++) {
                sum += tr_ST_full.get(i, j);
            }
            tr_ST_exit_full.set(i, 0, -sum);
        }

        // tr_ST_exit_mat
        Matrix tr_ST_exit_mat = new Matrix(Sc, Sc);
        for (int i = 0; i < Sc; i++) {
            for (int j = 0; j < Sc; j++) {
                tr_ST_exit_mat.set(i, j, tr_ST_exit_full.get(i, 0));
            }
        }

        // TS2 = (TS')*diag(alfa): TS2[a][b] = TS_normalized[b][a] * alfa[b].
        // (a indexes TS columns, b indexes the ma*mWait phases.)
        int nColsTS = TS_normalized.getNumCols();
        int nRowsTS = TS_normalized.getNumRows();  // = ma*mWait
        Matrix TS2 = new Matrix(nColsTS, nRowsTS);
        for (int a = 0; a < nColsTS; a++) {
            for (int b = 0; b < nRowsTS; b++) {
                TS2.set(a, b, TS_normalized.get(b, a) * alfa.get(0, b));
            }
        }

        // Normalize each row of TS2 by its row sum (MATLAB: TS2 ./ (sum(TS2,2)*ones)).
        Matrix tr_TS_ind2_normalized = new Matrix(nColsTS, nRowsTS);
        for (int a = 0; a < nColsTS; a++) {
            double rs = 0.0;
            for (int b = 0; b < nRowsTS; b++) {
                rs += TS2.get(a, b);
            }
            if (Math.abs(rs) < 1e-12) {
                rs = 1.0;
            }
            for (int b = 0; b < nRowsTS; b++) {
                tr_TS_ind2_normalized.set(a, b, TS2.get(a, b) / rs);
            }
        }

        // tildeP = [tr_ST_exit_mat .* tr_TS_ind2_normalized(1:Sc,1:Sc); zeros(Tr-Sc, Sc)]
        Matrix tildeP = new Matrix(Tr, Sc);
        for (int i = 0; i < Sc; i++) {
            for (int j = 0; j < Sc; j++) {
                tildeP.set(i, j, tr_ST_exit_mat.get(i, j) * tr_TS_ind2_normalized.get(i, j));
            }
        }

        // Extract non-zero phases
        Matrix tr_ST = new Matrix(m_tr_ST, m_tr_ST);
        Matrix tr_start_state = new Matrix(1, m_tr_ST);
        Matrix tildeP_reduced = new Matrix(m_tr_ST, Sc);

        int idx = 0;
        int[] indexMap = new int[Tr];
        for (int i = 0; i < Tr; i++) {
            if (busy_nz[i]) {
                indexMap[i] = idx;
                tr_start_state.set(0, idx, tr_start_state_full.get(0, i));
                idx++;
            }
        }

        idx = 0;
        for (int i = 0; i < Tr; i++) {
            if (busy_nz[i]) {
                int idx2 = 0;
                for (int j = 0; j < Tr; j++) {
                    if (busy_nz[j]) {
                        tr_ST.set(idx, idx2, tr_ST_full.get(i, j));
                        idx2++;
                    }
                }
                for (int j = 0; j < Sc; j++) {
                    tildeP_reduced.set(idx, j, tildeP.get(i, j));
                }
                idx++;
            }
        }

        // gamma_res
        Matrix gamma_res = new Matrix(1, dim_service + m_tr_ST + Sc);
        for (int i = 0; i < dim_service; i++) {
            gamma_res.set(0, i, notbusy_start.get(0, i));
        }
        for (int i = 0; i < m_tr_ST; i++) {
            gamma_res.set(0, dim_service + i, tr_start_state.get(0, i));
        }

        // C_res block matrix
        Matrix C_res = new Matrix(Tr + m_tr_ST + Sc, Tr + m_tr_ST + Sc);
        FJUtils.setSubMatrix(C_res, 0, 0, ST);
        FJUtils.setSubMatrix(C_res, Tr, Tr, tr_ST);
        FJUtils.setSubMatrix(C_res, Tr, Tr + m_tr_ST, tildeP_reduced);
        FJUtils.setSubMatrix(C_res, Tr + m_tr_ST, Tr + m_tr_ST, wait_Smat);

        return ReturnPer.returnPer(gamma_res, C_res, pers);
    }

    public static Matrix returnRT2(FJArrival arrival, FJService service, double[] pers, int C) {
        return returnRT2(arrival, service, pers, C, "NARE");
    }

    private static double elementSumLocal(Matrix m) {
        double sum = 0.0;
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                sum += m.get(i, j);
            }
        }
        return sum;
    }
}
