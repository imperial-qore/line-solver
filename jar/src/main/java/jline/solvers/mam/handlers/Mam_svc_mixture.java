/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mam.handlers;

import java.util.Map;

import jline.api.mc.Ctmc_solve;
import jline.api.qsys.QsysServiceLaw;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Arrival-weighted phase-type mixture of the per-class service laws at a
 * station, as the service descriptor accepted by
 * {@link jline.api.qsys.Qsys_mapg1k} and {@link jline.api.qsys.Qsys_mmapg1k}.
 *
 * <p>Port of {@code matlab/src/solvers/MAM/mam_svc_mixture.m}. The mixture is
 * PH(alpha_mix, T_mix) with alpha_mix = [w_1 pie_1, ...],
 * T_mix = blkdiag(D0_1, ...), and w_k = lambda_k/sum_j lambda_j the fraction of
 * arrivals belonging to class k. It is therefore the service law of an
 * arbitrary packet, and it reduces to the common law exactly, as a
 * distribution, when every class shares one.
 *
 * <p>This is the same construction {@link Mam_truncate_renorm} already applied
 * inline when a station carries more than one class, factored out so that the
 * exact finite-buffer branch and the truncate-and-renormalize fallback rest on
 * identical service assumptions and stay comparable.
 */
public final class Mam_svc_mixture {

    private Mam_svc_mixture() {
    }

    /** The mixture, kept as its two blocks plus the aggregated arrival process. */
    public static final class Result {
        /** 1 x n, the arrival-weighted initial probability row. */
        public final Matrix alpha;
        /** n x n, the block-diagonal subgenerator. */
        public final Matrix T;
        /** sum_k D_arr(k+1), the aggregate arrival matrix. */
        public final Matrix Dsum;
        /** Per-class arrival rates lambda_k. */
        public final double[] lambdaClass;

        Result(Matrix alpha, Matrix T, Matrix Dsum, double[] lambdaClass) {
            this.alpha = alpha;
            this.T = T;
            this.Dsum = Dsum;
            this.lambdaClass = lambdaClass;
        }

        /** The mixture as a service law descriptor. */
        public QsysServiceLaw toServiceLaw() {
            return QsysServiceLaw.phaseType(alpha, T);
        }
    }

    /**
     * @param D_arr    {D0, D_class1, ..., D_classR} arrival MMAP
     * @param pie_cell per-class PH initial distributions, keyed 0..R-1
     * @param D0_cell  per-class PH subgenerators, keyed 0..R-1
     */
    public static Result mam_svc_mixture(MatrixCell D_arr, Map<Integer, Matrix> pie_cell,
                                         Map<Integer, Matrix> D0_cell) {
        int nClasses = pie_cell.size();
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
            lambda_k[k] = theta.mult(D_arr.get(k + 1)).mult(e_arr).get(0, 0);
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
        return new Result(alpha_mix, T_mix, Dsum, lambda_k);
    }
}
