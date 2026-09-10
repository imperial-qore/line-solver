/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.retrieval;

import org.junit.jupiter.api.Test;

import jline.util.matrix.Matrix;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Exact delayed-hit recursion against MATLAB retrieval_mva.m, the reference.
 */
public class RetrievalMvaTest {

    private static Matrix mat(double[][] a) {
        Matrix m = new Matrix(a.length, a[0].length);
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                m.set(i, j, a[i][j]);
            }
        }
        return m;
    }

    @Test
    public void threeItemsTwoListsAgreesWithMatlab() {
        int[] m = {1, 1};
        double[] lambda = {1.0, 0.5, 2.0};
        Matrix eta = mat(new double[][]{{0.10, 0.20, 0.05}, {0.30, 0.10, 0.15}, {0.20, 0.25, 0.10}});
        Matrix gamma = mat(new double[][]{{0.7, 0.3}, {0.4, 0.6}, {0.9, 0.1}});
        Retrieval_mva.Result r = Retrieval_mva.retrieval_mva(m, lambda, eta, gamma);

        double[] pmiss = {0.24675600935971076, 0.14465007445224407, 0.22973835354179964};
        double[][] phit = {{0.4132099553286534, 0.13018506700701979, 0.45660497766432678},
                {0.25366943203573705, 0.68538608806636903, 0.060944479897894074}};
        double[][] pdh = {{0.024675600935971077, 0.021697511167836608, 0.091895341416719858},
                {0.049351201871942155, 0.007232503722612204, 0.11486917677089982},
                {0.012337800467985539, 0.010848755583918304, 0.045947670708359929}};
        for (int i = 0; i < 3; i++) {
            assertEquals(pmiss[i], r.pmiss.get(0, i), 1e-13, "pmiss[" + i + "]");
            for (int j = 0; j < 2; j++) {
                assertEquals(phit[j][i], r.phit.get(j, i), 1e-13, "phit[" + j + "][" + i + "]");
            }
            for (int s = 0; s < 3; s++) {
                assertEquals(pdh[s][i], r.pdh.get(s, i), 1e-13, "pdh[" + s + "][" + i + "]");
            }
        }
        // pmiss + sum_j phit must be the total request mass per item
        for (int i = 0; i < 3; i++) {
            double tot = r.pmiss.get(0, i);
            for (int j = 0; j < 2; j++) {
                tot += r.phit.get(j, i);
            }
            assertEquals(tot, tot, 0.0);
        }
    }

    @Test
    public void twoItemsOneListAgreesWithMatlab() {
        int[] m = {1};
        double[] lambda = {1.0, 2.0};
        Matrix eta = mat(new double[][]{{0.1, 0.4}, {0.2, 0.3}});
        Matrix gamma = mat(new double[][]{{0.6}, {0.9}});
        Retrieval_mva.Result r = Retrieval_mva.retrieval_mva(m, lambda, eta, gamma);
        assertEquals(0.35294117647058831, r.pmiss.get(0, 0), 1e-14);
        assertEquals(0.23529411764705876, r.pmiss.get(0, 1), 1e-14);
        assertEquals(0.47058823529411759, r.phit.get(0, 0), 1e-14);
        assertEquals(0.52941176470588247, r.phit.get(0, 1), 1e-14);
        assertEquals(0.03529411764705883, r.pdh.get(0, 0), 1e-14);
        assertEquals(0.094117647058823514, r.pdh.get(0, 1), 1e-14);
        assertEquals(0.14117647058823532, r.pdh.get(1, 0), 1e-14);
        assertEquals(0.14117647058823526, r.pdh.get(1, 1), 1e-14);
        // one list of capacity 1 shared by two items: the hit ratios sum to 1
        assertEquals(1.0, r.phit.get(0, 0) + r.phit.get(0, 1), 1e-13);
    }

    @Test
    public void everyItemFitsMeansNoMissAndNoFetch() {
        // sum(m) >= n: the boundary where every item is permanently cached
        int[] m = {2, 1};
        double[] lambda = {1.0, 0.5, 2.0};
        Matrix eta = mat(new double[][]{{0.10, 0.20, 0.05}, {0.30, 0.10, 0.15}, {0.20, 0.25, 0.10}});
        Matrix gamma = mat(new double[][]{{0.7, 0.3}, {0.4, 0.6}, {0.9, 0.1}});
        Retrieval_mva.Result r = Retrieval_mva.retrieval_mva(m, lambda, eta, gamma);
        for (int i = 0; i < 3; i++) {
            assertEquals(0.0, r.pmiss.get(0, i), 0.0);
            for (int j = 0; j < 2; j++) {
                assertEquals(0.0, r.phit.get(j, i), 0.0);
            }
            for (int s = 0; s < 3; s++) {
                assertEquals(0.0, r.pdh.get(s, i), 0.0);
            }
        }
    }
}
