/**
 * @file Argument normalisation shared by the shortest-job-next solvers
 *
 * Ported at parity from MATLAB matlab/src/api/pfqn/private/sjn_args.m.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.io.InputOutput;
import jline.util.matrix.Matrix;

final class SjnArgs {
    private SjnArgs() {}

    /** Read a Matrix laid out as a row or a column into a plain vector. */
    static double[] vector(Matrix v, int R, String what) {
        int n = v.getNumRows() * v.getNumCols();
        if (n != R) {
            throw new IllegalArgumentException("the " + what + " has " + n + " entries but the"
                    + " demand matrix has " + R + " classes");
        }
        double[] out = new double[R];
        int i = 0;
        for (int a = 0; a < v.getNumRows(); a++) {
            for (int b = 0; b < v.getNumCols(); b++) {
                out[i++] = v.get(a, b);
            }
        }
        return out;
    }

    static double[][] matrix(Matrix m, int M, int R) {
        if (m.getNumRows() != M || m.getNumCols() != R) {
            throw new IllegalArgumentException("expected a " + M + " by " + R + " matrix");
        }
        double[][] out = new double[M][R];
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < R; j++) {
                out[i][j] = m.get(i, j);
            }
        }
        return out;
    }

    static double[][] ones(int M, int R) {
        double[][] out = new double[M][R];
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < R; j++) {
                out[i][j] = 1;
            }
        }
        return out;
    }

    static double[] onesVector(int R) {
        double[] out = new double[R];
        for (int i = 0; i < R; i++) {
            out[i] = 1;
        }
        return out;
    }

    /**
     * Per-visit service times. The job size the discipline compares is one visit's service time,
     * not the demand accumulated over all visits.
     */
    static double[][] perVisit(double[][] L, double[][] V) {
        int M = L.length;
        int R = L[0].length;
        double[][] S = new double[M][R];
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < R; j++) {
                S[i][j] = V[i][j] > 0 ? L[i][j] / V[i][j] : 0;
            }
        }
        return S;
    }

    static int[] stationSet(int[] sjnset, int M) {
        if (sjnset == null) {
            return new int[0];
        }
        for (int i = 0; i < sjnset.length; i++) {
            if (sjnset[i] < 0 || sjnset[i] >= M) {
                throw new IllegalArgumentException("sjnset contains a station index outside the"
                        + " demand matrix");
            }
            for (int j = i + 1; j < sjnset.length; j++) {
                if (sjnset[i] == sjnset[j]) {
                    throw new IllegalArgumentException("sjnset repeats a station index");
                }
            }
        }
        return sjnset.clone();
    }

    static int indexOf(int[] set, int value) {
        for (int i = 0; i < set.length; i++) {
            if (set[i] == value) {
                return i;
            }
        }
        return -1;
    }

    /** Mixed-radix decoding of a linear lattice index into a population vector. */
    static double[] decode(int idx, int[] stride, double[] N) {
        int R = N.length;
        double[] n = new double[R];
        for (int r = 0; r < R; r++) {
            n[r] = (idx / stride[r]) % ((int) N[r] + 1);
        }
        return n;
    }

    static Matrix rowMatrix(double[] v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) {
            m.set(0, i, v[i]);
        }
        return m;
    }

    static Matrix toMatrix(double[][] a) {
        Matrix m = new Matrix(a.length, a[0].length);
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                m.set(i, j, a[i][j]);
            }
        }
        return m;
    }

    static double[][][] copy3(double[][][] a) {
        double[][][] out = new double[a.length][][];
        for (int i = 0; i < a.length; i++) {
            out[i] = copy2(a[i]);
        }
        return out;
    }

    static double[][] copy2(double[][] a) {
        double[][] out = new double[a.length][];
        for (int i = 0; i < a.length; i++) {
            out[i] = a[i].clone();
        }
        return out;
    }

    static void warnCapped(double umax) {
        InputOutput.line_warning(SjnArgs.class.getName(), "the utilization cap of " + umax
                + " was binding at an SJN station: the station is in the starvation regime, where"
                + " long jobs are held back and the arrival theorem is badly violated. The results"
                + " are stable but their accuracy is not warranted, use SolverCTMC or SolverLDES"
                + " there.");
    }
}
