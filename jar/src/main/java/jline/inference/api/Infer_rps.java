/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference.api;

import jline.inference.util.NnlsSolver;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.List;

public final class Infer_rps {
    private Infer_rps() {}

    /**
     * Regression for Processor Sharing (RPS) demand estimation.
     *
     * Based on mean-value analysis for PS stations:
     *   E[R_r] = E[D_r] * E[Q_bar_A] / V
     *
     * where Q_bar_A is the total number of jobs seen upon admission
     * (including the arriving job) and V is the number of servers.
     *
     * @param rt response time samples (column vector)
     * @param classVec class of each request sample (0-based)
     * @param ql queue length samples (n x R matrix)
     * @param V number of cores/servers
     * @return estimated demands (1 x R array)
     */
    public static double[] infer_rps(double[] rt, int[] classVec, Matrix ql, int V) {
        int R = ql.getNumCols();
        double[] demandEst = new double[R];

        for (int r = 0; r < R; r++) {
            List<Integer> indices = new ArrayList<Integer>();
            for (int k = 0; k < classVec.length; k++) {
                if (classVec[k] == r) indices.add(k);
            }
            if (indices.isEmpty()) continue;

            double[] respTimes = new double[indices.size()];
            for (int k = 0; k < indices.size(); k++) {
                respTimes[k] = rt[indices.get(k)];
            }
            double[] qBarA = new double[indices.size()];
            for (int k = 0; k < indices.size(); k++) {
                int idx = indices.get(k);
                double totalQL = 0.0;
                for (int c = 0; c < R; c++) {
                    totalQL += ql.get(idx, c);
                }
                qBarA[k] = (totalQL + 1.0) / V;
            }

            // NNLS: min ||qBarA * d - respTimes||^2 s.t. d >= 0
            double[][] A = new double[respTimes.length][1];
            for (int k = 0; k < respTimes.length; k++) {
                A[k][0] = qBarA[k];
            }
            double[] result = NnlsSolver.lsqnonneg(A, respTimes);
            demandEst[r] = result[0];
        }

        return demandEst;
    }
}
