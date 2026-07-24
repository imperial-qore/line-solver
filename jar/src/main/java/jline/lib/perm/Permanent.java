/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.perm;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.util.Pair;
import jline.util.matrix.Matrix;

/**
 * Implementation of the MATLAB perm.m permanent computation algorithm.
 */
public class Permanent extends PermSolver {

    public Permanent(Matrix matrix) {
        this(matrix, false);
    }

    public Permanent(Matrix matrix, boolean solve) {
        super(matrix);
        if (solve) {
            solve();
        }
    }

    @Override
    public void compute() {
        value = computeWithMultiplicities();
    }

    /**
     * Computes the permanent using the MATLAB perm.m algorithm with multiplicities.
     */
    private double computeWithMultiplicities() {
        Pair<Matrix, int[]> uniqueResult = findUniqueColumnsWithMultiplicities();
        Matrix uniqueMatrix = uniqueResult.getLeft();
        int[] multiplicities = uniqueResult.getRight();

        int R = multiplicities.length;
        int n = 0;
        for (int m : multiplicities) n += m;
        double result = 0.0;

        int[] f = new int[R];

        boolean done = false;
        while (!done) {
            int fSum = 0;
            for (int v : f) fSum += v;
            double term = Math.pow(-1.0, fSum);

            for (int j = 0; j < R; j++) {
                term *= binomialCoefficient(multiplicities[j], f[j]);
            }

            for (int i = 0; i < n; i++) {
                double sumTerm = 0.0;
                for (int k = 0; k < R; k++) {
                    sumTerm += f[k] * uniqueMatrix.get(i, k);
                }
                term *= sumTerm;
            }

            result += term;
            f = pprodNext(f, multiplicities);
            if (f.length == 0) {
                done = true;
            } else {
                boolean allMinus1 = true;
                for (int v : f) {
                    if (v != -1) { allMinus1 = false; break; }
                }
                if (allMinus1) done = true;
            }
        }

        return Math.pow(-1.0, n) * result;
    }

    /**
     * Finds unique columns in the matrix and their multiplicities.
     */
    private Pair<Matrix, int[]> findUniqueColumnsWithMultiplicities() {
        Map<List<Double>, Integer> columnMap = new HashMap<List<Double>, Integer>();
        List<List<Double>> uniqueColumns = new ArrayList<List<Double>>();
        List<Integer> multiplicities = new ArrayList<Integer>();

        for (int j = 0; j < matrix.getNumCols(); j++) {
            List<Double> column = new ArrayList<Double>(matrix.getNumRows());
            for (int i = 0; i < matrix.getNumRows(); i++) {
                column.add(matrix.get(i, j));
            }
            Integer index = columnMap.get(column);
            if (index == null) {
                columnMap.put(column, uniqueColumns.size());
                uniqueColumns.add(column);
                multiplicities.add(1);
            } else {
                multiplicities.set(index, multiplicities.get(index) + 1);
            }
        }

        boolean allOne = true;
        for (int m : multiplicities) {
            if (m != 1) { allOne = false; break; }
        }

        if (allOne) {
            Map<List<Double>, Integer> rowMap = new HashMap<List<Double>, Integer>();
            List<List<Double>> uniqueRows = new ArrayList<List<Double>>();
            List<Integer> rowMultiplicities = new ArrayList<Integer>();

            for (int i = 0; i < matrix.getNumRows(); i++) {
                List<Double> row = new ArrayList<Double>(matrix.getNumCols());
                for (int j = 0; j < matrix.getNumCols(); j++) {
                    row.add(matrix.get(i, j));
                }
                Integer index = rowMap.get(row);
                if (index == null) {
                    rowMap.put(row, uniqueRows.size());
                    uniqueRows.add(row);
                    rowMultiplicities.add(1);
                } else {
                    rowMultiplicities.set(index, rowMultiplicities.get(index) + 1);
                }
            }

            boolean anyRepeated = false;
            for (int m : rowMultiplicities) {
                if (m > 1) { anyRepeated = true; break; }
            }

            if (anyRepeated) {
                Matrix transposedMatrix = new Matrix(uniqueRows.get(0).size(), uniqueRows.size());
                for (int i = 0; i < uniqueRows.size(); i++) {
                    List<Double> row = uniqueRows.get(i);
                    for (int j = 0; j < row.size(); j++) {
                        transposedMatrix.set(j, i, row.get(j));
                    }
                }
                int[] arr = new int[rowMultiplicities.size()];
                for (int k = 0; k < arr.length; k++) arr[k] = rowMultiplicities.get(k);
                return new Pair<Matrix, int[]>(transposedMatrix, arr);
            } else {
                Matrix uniqueMatrix = new Matrix(matrix.getNumRows(), uniqueColumns.size());
                for (int j = 0; j < uniqueColumns.size(); j++) {
                    List<Double> col = uniqueColumns.get(j);
                    for (int i = 0; i < col.size(); i++) {
                        uniqueMatrix.set(i, j, col.get(i));
                    }
                }
                int[] arr = new int[multiplicities.size()];
                for (int k = 0; k < arr.length; k++) arr[k] = multiplicities.get(k);
                return new Pair<Matrix, int[]>(uniqueMatrix, arr);
            }
        } else {
            Matrix uniqueMatrix = new Matrix(matrix.getNumRows(), uniqueColumns.size());
            for (int j = 0; j < uniqueColumns.size(); j++) {
                List<Double> col = uniqueColumns.get(j);
                for (int i = 0; i < col.size(); i++) {
                    uniqueMatrix.set(i, j, col.get(i));
                }
            }
            int[] arr = new int[multiplicities.size()];
            for (int k = 0; k < arr.length; k++) arr[k] = multiplicities.get(k);
            return new Pair<Matrix, int[]>(uniqueMatrix, arr);
        }
    }

    /**
     * Computes binomial coefficient C(n, k).
     */
    private double binomialCoefficient(int n, int k) {
        if (k > n || k < 0) return 0.0;
        if (k == 0 || k == n) return 1.0;

        double result = 1.0;
        int upper = Math.min(k, n - k);
        for (int i = 1; i <= upper; i++) {
            result = result * (n - i + 1) / i;
        }
        return result;
    }

    /**
     * MATLAB pprod iterator - generates the next state in the sequence.
     */
    private int[] pprodNext(int[] current, int[] bounds) {
        int[] n = current.clone();
        int[] N = bounds;
        int R = N.length;

        boolean atMax = true;
        for (int i = 0; i < R; i++) {
            if (n[i] != N[i]) { atMax = false; break; }
        }
        if (atMax) {
            return new int[] { -1 };
        }

        int s = R - 1;
        while (s >= 0 && n[s] == N[s]) {
            n[s] = 0;
            s--;
        }

        if (s < 0) {
            return new int[] { -1 };
        }

        n[s]++;
        return n;
    }
}
