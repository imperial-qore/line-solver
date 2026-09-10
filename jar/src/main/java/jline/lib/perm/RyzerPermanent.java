/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.perm;

import java.util.ArrayList;
import java.util.List;

import jline.util.matrix.Matrix;

/**
 * Implementation of Ryzer's algorithm to calculate the permanent.
 */
public class RyzerPermanent extends PermSolver {
    private final String mode;

    public RyzerPermanent(Matrix matrix) {
        this(matrix, "graycode", false);
    }

    public RyzerPermanent(Matrix matrix, String mode) {
        this(matrix, mode, false);
    }

    public RyzerPermanent(Matrix matrix, String mode, boolean solve) {
        super(matrix);
        this.mode = mode;
        if (solve) {
            solve();
        }
    }

    @Override
    public void compute() {
        if ("graycode".equals(mode)) {
            value = ryzerAlgorithmGraycode();
        } else {
            value = ryzerAlgorithmNaive();
        }
    }

    /**
     * Gray code version of Ryzer's algorithm.
     */
    private double ryzerAlgorithmGraycode() {
        double permanent = 0.0;

        double[] rowSum = new double[n];

        List<Integer> bitToModifyList = bitToModify(n);

        boolean[] currentBit = new boolean[n];

        for (int bitIndex : bitToModifyList) {
            currentBit[bitIndex] = !currentBit[bitIndex];

            double multiplier = currentBit[bitIndex] ? 1.0 : -1.0;
            for (int i = 0; i < n; i++) {
                rowSum[i] += multiplier * matrix.get(i, bitIndex);
            }

            int nbCol = 0;
            for (boolean b : currentBit) {
                if (b) nbCol++;
            }

            double product = 1.0;
            for (int i = 0; i < n; i++) {
                product *= rowSum[i];
            }

            permanent += Math.pow(-1.0, nbCol) * product;
        }

        permanent *= Math.pow(-1.0, n);
        return permanent;
    }

    /**
     * Naive version of Ryzer's algorithm.
     */
    private double ryzerAlgorithmNaive() {
        double permanent = 0.0;

        for (int i = 0; i <= n; i++) {
            List<List<Integer>> combinations = generateCombinations(n, n - i);

            for (List<Integer> combination : combinations) {
                double[] rowSums = new double[n];
                for (int row = 0; row < n; row++) {
                    double s = 0.0;
                    for (int col : combination) {
                        s += matrix.get(row, col);
                    }
                    rowSums[row] = s;
                }

                double product = 1.0;
                for (double sum : rowSums) {
                    product *= sum;
                }

                permanent += Math.pow(-1.0, n - i) * product;
            }
        }

        permanent *= Math.pow(-1.0, n);
        return permanent;
    }

    /**
     * Create the order of bits to modify to satisfy Gray code ordering.
     */
    private List<Integer> bitToModify(int m) {
        List<Integer> result = new ArrayList<Integer>();
        if (m == 1) {
            result.add(0);
            return result;
        }
        List<Integer> sub = bitToModify(m - 1);
        result.addAll(sub);
        result.add(m - 1);
        result.addAll(sub);
        return result;
    }

    /**
     * Generate all combinations of k elements from n elements.
     */
    private List<List<Integer>> generateCombinations(int n, int k) {
        List<List<Integer>> result = new ArrayList<List<Integer>>();
        if (k == 0) {
            result.add(new ArrayList<Integer>());
            return result;
        }
        if (k > n) return result;

        List<Integer> current = new ArrayList<Integer>();
        backtrack(result, current, 0, n, k);
        return result;
    }

    private void backtrack(List<List<Integer>> result, List<Integer> current, int start, int n, int k) {
        if (current.size() == k) {
            result.add(new ArrayList<Integer>(current));
            return;
        }
        for (int i = start; i < n; i++) {
            current.add(i);
            backtrack(result, current, i + 1, n, k);
            current.remove(current.size() - 1);
        }
    }
}
