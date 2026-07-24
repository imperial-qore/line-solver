/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference.api;

import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.List;

public final class Infer_get_qlen_arrival {
    private Infer_get_qlen_arrival() {}

    /**
     * Compute queue lengths at arrival from legacy cell data format.
     *
     * Wrapper around infer_compute_ql_at_arrival for the legacy cell-based
     * data format. data[3][k] contains arrival times (in ms) and data[4][k]
     * contains response times for class k.
     *
     * @param data cell array in standard format: data[row][class],
     *             data[3][k] = arrival times (ms), data[4][k] = response times (s)
     * @param K number of classes
     * @return list of K matrices, each numSamples(k) x K
     */
    public static List<Matrix> infer_get_qlen_arrival(double[][][] data, int K) {
        List<Double> atList = new ArrayList<Double>();
        List<Double> rtList = new ArrayList<Double>();
        List<Integer> classList = new ArrayList<Integer>();
        int[] numObs = new int[K];

        for (int k = 0; k < K; k++) {
            double[] arvTimes = data[3][k];
            double[] respTimes = data[4][k];
            if (arvTimes == null || respTimes == null) continue;
            numObs[k] = arvTimes.length;
            for (int i = 0; i < arvTimes.length; i++) {
                atList.add(arvTimes[i] / 1000.0); // convert ms to secs
                rtList.add(respTimes[i]);
                classList.add(k);
            }
        }

        double[] at = new double[atList.size()];
        for (int i = 0; i < atList.size(); i++) at[i] = atList.get(i).doubleValue();
        double[] rt = new double[rtList.size()];
        for (int i = 0; i < rtList.size(); i++) rt[i] = rtList.get(i).doubleValue();
        int[] classVec = new int[classList.size()];
        for (int i = 0; i < classList.size(); i++) classVec[i] = classList.get(i).intValue();
        int n = at.length;
        int[] jobid = new int[n];
        for (int i = 0; i < n; i++) jobid[i] = i;

        Matrix qlUnsorted = Infer_compute_ql_at_arrival.infer_compute_ql_at_arrival(at, jobid, rt, jobid, classVec, K);

        // Split into per-class matrices in original order
        List<Matrix> ql = new ArrayList<Matrix>(K);
        int counter = 0;
        for (int k = 0; k < K; k++) {
            Matrix mat = new Matrix(numObs[k], K);
            for (int i = 0; i < numObs[k]; i++) {
                for (int r = 0; r < K; r++) {
                    mat.set(i, r, qlUnsorted.get(counter + i, r));
                }
            }
            ql.add(mat);
            counter += numObs[k];
        }

        return ql;
    }
}
