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
     * Log10 magnitude of the largest Ryser term for this orientation.
     *
     * The inclusion-exclusion expansion is largest when every column is selected,
     * giving prod_i (sum_j a_ij). Since per(A) equals per(A transposed), the two
     * orientations return the same value but not the same cancellation, so this
     * is the quantity to minimize when choosing between them. Returns positive
     * infinity when a row sum vanishes, so such an orientation is never chosen.
     */
    private static double ryserConditioning(Matrix m) {
        double total = 0.0;
        for (int i = 0; i < m.getNumRows(); i++) {
            double lineSum = 0.0;
            for (int j = 0; j < m.getNumCols(); j++) {
                lineSum += Math.abs(m.get(i, j));
            }
            if (lineSum <= 0.0) {
                return Double.POSITIVE_INFINITY;
            }
            total += Math.log10(lineSum);
        }
        return total;
    }

    /**
     * Groups repeated columns of the better-conditioned orientation of the matrix.
     *
     * Exploiting repeated rows means expanding the transpose, which leaves the
     * permanent unchanged but can raise the largest intermediate term by many
     * orders of magnitude. Vandermonde-like matrices, such as the A_x of
     * pfqn_lcfsqn_nc whose rows are geometric in the column index, lose every
     * significant digit that way. Orientation is chosen by conditioning first
     * and grouping applied second, even when that forgoes the grouping.
     */
    private Pair<Matrix, int[]> findUniqueColumnsWithMultiplicities() {
        Matrix source = matrix;
        Matrix transposed = matrix.transpose();
        if (ryserConditioning(transposed) < ryserConditioning(matrix)) {
            source = transposed;
        }

        Map<List<Double>, Integer> columnMap = new HashMap<List<Double>, Integer>();
        List<List<Double>> uniqueColumns = new ArrayList<List<Double>>();
        List<Integer> multiplicities = new ArrayList<Integer>();

        for (int j = 0; j < source.getNumCols(); j++) {
            List<Double> column = new ArrayList<Double>(source.getNumRows());
            for (int i = 0; i < source.getNumRows(); i++) {
                column.add(source.get(i, j));
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

        Matrix uniqueMatrix = new Matrix(source.getNumRows(), uniqueColumns.size());
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
