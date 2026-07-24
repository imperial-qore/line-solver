/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference.api;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import jline.lang.Network;
import jline.lang.nodes.Queue;
import jline.util.matrix.Matrix;

public final class Infer_minps_setup {
    private Infer_minps_setup() {}

    /**
     * Setup function for MINPS with legacy cell-based data format.
     */
    public static double[] infer_minps_setup(double[][][] data,
                                             int initSample,
                                             int sampleSize,
                                             int V,
                                             Network model,
                                             Queue node) {
        int R = data[0].length - 1;

        // Count samples per class
        int[] sampleNumber = new int[R];
        for (int k = 0; k < R; k++) {
            sampleNumber[k] = (data[3][k] != null) ? data[3][k].length : 0;
        }

        // Remove classes without samples
        int newR = 0;
        for (int k = 0; k < R; k++) {
            if (sampleNumber[k] > 0) newR++;
        }
        double[][][] data2 = new double[7][newR + 1][];
        int rIdx = 0;
        for (int k = 0; k <= R; k++) {
            if ((k < R && sampleNumber[k] > 0) || k == R) {
                for (int j = 0; j < 7; j++) {
                    data2[j][rIdx] = data[j][k];
                }
                rIdx++;
            }
        }
        R = newR;

        // Get queue length at arrival times
        List<Matrix> qls = Infer_get_qlen_arrival.infer_get_qlen_arrival(data2, R);

        List<Double> rtList = new ArrayList<Double>();
        List<Integer> classList = new ArrayList<Integer>();
        List<double[]> qlList = new ArrayList<double[]>();
        List<Double> atList = new ArrayList<Double>();

        for (int k = 0; k < R; k++) {
            double[] respTimes = data2[4][k];
            double[] arvTimes = data2[3][k];
            if (respTimes == null || arvTimes == null) continue;
            Matrix qlk = qls.get(k);
            for (int i = 0; i < respTimes.length; i++) {
                rtList.add(respTimes[i]);
                classList.add(k);
                double[] qlRow = new double[R];
                for (int r = 0; r < R; r++) {
                    qlRow[r] = qlk.get(i, r);
                }
                qlList.add(qlRow);
                atList.add(arvTimes[i] / 1000.0); // ms to seconds
            }
        }

        // Sort by arrival time
        Integer[] indicesObj = new Integer[atList.size()];
        for (int i = 0; i < atList.size(); i++) indicesObj[i] = i;
        final List<Double> atListFinal = atList;
        java.util.Arrays.sort(indicesObj, new Comparator<Integer>() {
            @Override
            public int compare(Integer a, Integer b) {
                return Double.compare(atListFinal.get(a), atListFinal.get(b));
            }
        });
        int n = indicesObj.length;
        int[] indices = new int[n];
        for (int i = 0; i < n; i++) indices[i] = indicesObj[i];

        double[] rt = new double[n];
        int[] classArr = new int[n];
        Matrix qlMatrix = new Matrix(n, R);
        for (int i = 0; i < n; i++) {
            rt[i] = rtList.get(indices[i]);
            classArr[i] = classList.get(indices[i]);
            for (int r = 0; r < R; r++) {
                qlMatrix.set(i, r, qlList.get(indices[i])[r]);
            }
        }

        int actualSampleSize = (sampleSize == 0) ? qlMatrix.getNumRows() : sampleSize;
        int firstSample = initSample;
        int finalSample = Math.min(firstSample + actualSampleSize, qlMatrix.getNumRows());

        int sliceLen = finalSample - firstSample;
        double[] rtExp = new double[sliceLen];
        int[] classExp = new int[sliceLen];
        Matrix qlExp = new Matrix(sliceLen, R);
        for (int i = 0; i < sliceLen; i++) {
            rtExp[i] = rt[firstSample + i];
            classExp[i] = classArr[firstSample + i];
            for (int r = 0; r < R; r++) {
                qlExp.set(i, r, qlMatrix.get(firstSample + i, r));
            }
        }

        // Remove zero response times
        List<Integer> valid = new ArrayList<Integer>();
        for (int i = 0; i < rtExp.length; i++) {
            if (rtExp[i] > 0) valid.add(i);
        }
        double[] rtValid = new double[valid.size()];
        int[] classValid = new int[valid.size()];
        Matrix qlValid = new Matrix(valid.size(), R);
        for (int i = 0; i < valid.size(); i++) {
            int vi = valid.get(i);
            rtValid[i] = rtExp[vi];
            classValid[i] = classExp[vi];
            for (int r = 0; r < R; r++) {
                qlValid.set(i, r, qlExp.get(vi, r));
            }
        }
        // suppress unused warnings
        Collections.<Integer>emptyList();
        return Infer_minps.infer_minps(model, node, rtValid, classValid, qlValid);
    }
}
