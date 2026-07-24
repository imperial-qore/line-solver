/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mam.handlers;

import java.util.HashMap;
import java.util.Map;

import jline.api.mc.Ctmc_solve;
import jline.lib.butools.MMAPPH1FCFS;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Finite-buffer marginal for MMAP[K]/PH[K]/1/FCFS.
 *
 * <p>Port of matlab/src/solvers/MAM/mam_truncate_renorm.m. Solves the
 * infinite-buffer MMAP[K]/PH[K]/1/FCFS queue via BuTools MMAPPH1FCFS, truncates
 * the marginal queue length distribution at the buffer capacity, and
 * renormalizes. For an M/M/1 input the renormalized distribution coincides
 * exactly with the M/M/1/K marginal; for general MMAP/PH it is an ASTA-style
 * approximation.</p>
 */
public final class Mam_truncate_renorm {
    private Mam_truncate_renorm() {}

    /**
     * Outcome of the finite-buffer truncation.
     */
    public static final class Result {
        /** mean number of jobs in system, clipped to [0, capK]. */
        public final double meanQ;
        /** blocking probability, the renormalized boundary mass p(N=capK). */
        public final double lossProb;
        /** renormalized truncated marginal, 1x(capK+1). */
        public final Matrix p_norm;

        Result(double meanQ, double lossProb, Matrix p_norm) {
            this.meanQ = meanQ;
            this.lossProb = lossProb;
            this.p_norm = p_norm;
        }
    }

    /**
     * @param D_arr    {D0, D_class1, ..., D_classK} arrival MMAP
     * @param pie_cell per-class PH initial distributions
     * @param D0_cell  per-class PH subgenerators
     * @param capK     buffer capacity (max jobs in system)
     */
    public static Result mam_truncate_renorm(MatrixCell D_arr,
                                             Map<Integer, Matrix> pie_cell,
                                             Map<Integer, Matrix> D0_cell,
                                             int capK) {
        int nClasses = pie_cell.size();

        MatrixCell D_call;
        Map<Integer, Matrix> pie_call;
        Map<Integer, Matrix> D0_call;

        if (nClasses > 1) {
            // see _kb/06-solver-catalog.md for rationale
            Matrix D0 = D_arr.get(0);
            Matrix Dsum = new Matrix(D0.getNumRows(), D0.getNumCols());
            for (int k = 0; k < nClasses; k++) {
                Dsum = Dsum.add(1.0, D_arr.get(k + 1));
            }
            Matrix e_arr = Matrix.ones(D0.getNumRows(), 1);
            Matrix theta = Ctmc_solve.ctmc_solve(D0.add(1.0, Dsum));

            double[] lambda_k = new double[nClasses];
            double sumL = 0.0;
            for (int k = 0; k < nClasses; k++) {
                Matrix v = theta.mult(D_arr.get(k + 1)).mult(e_arr);
                lambda_k[k] = v.get(0, 0);
                sumL += lambda_k[k];
            }
            double[] w = new double[nClasses];
            for (int k = 0; k < nClasses; k++) {
                w[k] = (sumL > 0) ? lambda_k[k] / sumL : 1.0 / nClasses;
            }

            int n_total = 0;
            int[] n_k = new int[nClasses];
            for (int k = 0; k < nClasses; k++) {
                n_k[k] = pie_cell.get(k).length();
                n_total += n_k[k];
            }
            Matrix alpha_mix = new Matrix(1, n_total);
            Matrix T_mix = new Matrix(n_total, n_total);
            int offset = 0;
            for (int k = 0; k < nClasses; k++) {
                Matrix pk = pie_cell.get(k);
                for (int i = 0; i < n_k[k]; i++) {
                    alpha_mix.set(0, offset + i, w[k] * pk.get(i));
                }
                Matrix Tk = D0_cell.get(k);
                for (int i = 0; i < n_k[k]; i++) {
                    for (int j = 0; j < n_k[k]; j++) {
                        T_mix.set(offset + i, offset + j, Tk.get(i, j));
                    }
                }
                offset += n_k[k];
            }

            D_call = new MatrixCell(2);
            D_call.set(0, D0);
            D_call.set(1, Dsum);
            pie_call = new HashMap<Integer, Matrix>();
            pie_call.put(0, alpha_mix);
            D0_call = new HashMap<Integer, Matrix>();
            D0_call.put(0, T_mix);
        } else {
            D_call = D_arr;
            pie_call = pie_cell;
            D0_call = D0_cell;
        }

        int nLevels = capK + 1;
        Map<String, Map<Integer, Matrix>> res = MMAPPH1FCFS.MMAPPH1FCFS(
                D_call, pie_call, D0_call, null, Integer.valueOf(nLevels),
                null, null, false, false, null, null);

        Matrix pdistrM = null;
        if (res.containsKey("ncDistr") && res.get("ncDistr").containsKey(0)) {
            pdistrM = res.get("ncDistr").get(0);
        }

        double[] p_in = new double[nLevels];
        double massIn = 0.0;
        if (pdistrM != null) {
            for (int i = 0; i < nLevels && i < pdistrM.length(); i++) {
                p_in[i] = Math.abs(pdistrM.get(i));
                massIn += p_in[i];
            }
        }

        Matrix p_norm = new Matrix(1, nLevels);
        if (massIn <= 0) {
            p_norm.set(0, 0, 1.0);
        } else {
            for (int i = 0; i < nLevels; i++) {
                p_norm.set(0, i, p_in[i] / massIn);
            }
        }

        double meanQ = 0.0;
        for (int i = 0; i <= capK; i++) {
            meanQ += i * p_norm.get(0, i);
        }
        meanQ = Math.max(0.0, Math.min((double) capK, meanQ));
        double lossProb = p_norm.get(0, nLevels - 1);

        return new Result(meanQ, lossProb, p_norm);
    }
}
