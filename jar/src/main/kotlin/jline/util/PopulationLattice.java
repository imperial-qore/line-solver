/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

// Copyright (c) 2012-2026, Imperial College London
// All rights reserved.

package jline.util;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import static jline.util.Maths.multichoose;

/**
 * Data structure modeling a lattice used to describe a combination of job populations.
 */
public class PopulationLattice {

    /**
     * Computes a hash index for a population vector within a lattice.
     * Maps a population vector to a unique index in the lattice space.
     * 
     * @param n the current population vector
     * @param N the maximum population vector (lattice bounds)
     * @return unique hash index for the population state
     */
    public static int hashpop(Matrix n, Matrix N) {
        int idx = 0;
        int R = N.length();
        for (int r = 0; r < R; r++) {
            double prod = 1;
            for (int j = 0; j < r; j++) {
                prod *= (N.get(j) + 1);
            }
            idx += prod * n.get(r);
        }
        return idx;
    }

    /**
     * Computes a hash index for a population vector using precomputed products.
     * Optimized version that uses precomputed product values for efficiency.
     * 
     * @param n the current population vector
     * @param N the maximum population vector (unused in this overload)
     * @param R the number of job classes
     * @param prods precomputed product values for hash computation
     * @return unique hash index for the population state
     */
    public static int hashpop(Matrix n, Matrix N, int R, Matrix prods) {
        int idx = 0;
        for (int r = 0; r < R; r++) {
            idx += prods.get(r) * n.get(r);
        }
        return idx;
    }

    /**
     * Initializes a population product iterator.
     * Returns a zero matrix to start iterating through population states.
     * 
     * @param n the maximum population vector
     * @return initial state (zero matrix) for population iteration
     */
    public static Matrix pprod(Matrix n) {
        return new Matrix(n.getNumRows(), n.getNumCols());
    }

    /**
     * Advances to the next population state in lexicographic order.
     * Generates the next non-negative vector less than or equal to N.
     * 
     * @param n the current population state
     * @param N the maximum population bounds
     * @return next population state, or null if enumeration is complete
     */
    public static Matrix pprod(Matrix n, Matrix N) {
        if (N.isEmpty()) {
            N = n.copy();
            n = new Matrix(N.getNumRows(), N.getNumCols());
            return n;
        }

        int R = N.length();
        int countEqual = 0;
        // Count equal elements treating input as a vector
        for (int i = 0; i < R; i++) {
            if (n.get(i) == N.get(i)) {
                countEqual++;
            }
        }
        if (countEqual == R) {
            n = new Matrix(1, 1);
            n.set(0, 0, -1);
            return n;
        }

        int s = R;
        while (s > 0 && n.get(s - 1) == N.get(s - 1)) {
            n.set(s - 1, 0);
            s--;
        }

        if (s == 0) {
            return n;
        }

        n.set(s - 1, n.get(s - 1) + 1);
        return n;
    }

    public static Matrix pprodcon(Matrix n, Matrix lb, Matrix ub) {
        if (n == null) {
            lb = new Matrix(lb.getNumRows(), lb.getNumCols());
            n = new Matrix(lb.getNumRows(), lb.getNumCols());
            return n;
        }
        int R = ub.length();
        int sum = 0;
        for (int row = 0; row < n.getNumRows(); row++) {
            for (int col = 0; col < n.getNumCols(); col++) {
                if (n.get(row, col) == ub.get(row, col)) {
                    sum += 1;
                }
            }
        }
        if (sum == R) {
            for (int row = 0; row < n.getNumRows(); row++) {
                for (int col = 0; col < n.getNumCols(); col++) {
                    n.set(row, col, -1);
                }
            }
            return n;
        }
        int s = R;
        while (s > 0 && n.get(s - 1) == ub.get(s - 1)) {
            n.set(s - 1, lb.get(s - 1));
            s = s - 1;
        }
        if (s == 0) {
            return n;
        }
        n.set(s - 1, n.get(s - 1) + 1);
        return n;
    }

    public static sprodResult sprod(int M, Matrix N) {
        int R = N.length();
        Matrix S = new Matrix(1, R);
        MatrixCell D = new MatrixCell(R);
        for (int r = 0; r < R; r++) {
            D.set(r, multichoose(M, N.get(r)));
            S.set(r, D.get(r).getNumRows() - 1);
        }
        Matrix s = pprod(S);
        Matrix n = new Matrix(D.get(0).getNumCols(), R);
        for (int r = 0; r < R; r++) {
            n.setColumn(r, D.get(r).getRow((int) (s.get(r))).transpose());
        }
        return new sprodResult(s, n, S, D);
    }

    public static sprodResult sprod(Matrix s, Matrix S, MatrixCell D) {
        int M = D.get(0).getNumCols();
        int R = D.size();
        Matrix n = Matrix.zeros(M, R);

        s = pprod(s, S);
        if (s.value() == -1) {
            n = n.sub(Matrix.createLike(n));
            return new sprodResult(s, n, S, D);
        }

        for (int r = 0; r < R; r++) {
            n.setColumn(r, D.get(r).getRow((int) (s.get(r))).transpose());
        }
        return new sprodResult(s, n, S, D);
    }

    public static class sprodResult {
        public Matrix s;
        public Matrix n;
        public Matrix S;
        public MatrixCell D;

        public sprodResult(Matrix s, Matrix n, Matrix S, MatrixCell D) {
            this.s = s;
            this.n = n;
            this.S = S;
            this.D = D;
        }
    }
}
