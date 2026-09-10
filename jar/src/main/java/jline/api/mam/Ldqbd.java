/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.mam;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import jline.util.matrix.Matrix;

/**
 * Level-Dependent Quasi-Birth-Death process solver.
 */
public final class Ldqbd {
    private Ldqbd() {}

    public static LdqbdResult ldqbd(List<Matrix> Q0, List<Matrix> Q1, List<Matrix> Q2) {
        return ldqbd(Q0, Q1, Q2, new LdqbdOptions());
    }

    public static LdqbdResult ldqbd(List<Matrix> Q0, List<Matrix> Q1, List<Matrix> Q2, LdqbdOptions options) {
        int N = Q1.size() - 1;
        if (options.getVerbose()) {
            System.out.println("LD-QBD Solver: N=" + N + " levels, epsilon=" + options.getEpsilon());
        }
        List<Matrix> R = computeAllRateMatrices(N, Q0, Q1, Q2, options);
        if (options.getVerbose()) {
            for (int n = 0; n < N; n++) {
                System.out.println("  R^(" + (n + 1) + ") computed (" + R.get(n).getNumRows() + "x" + R.get(n).getNumCols() + ")");
            }
        }
        List<Matrix> piCells = new ArrayList<Matrix>();
        Matrix pi = computeStationaryDist(R, Q0, Q1, Q2, N, options, piCells);
        return new LdqbdResult(R, pi, piCells.isEmpty() ? null : piCells);
    }

    private static List<Matrix> computeAllRateMatrices(int N, List<Matrix> Q0, List<Matrix> Q1, List<Matrix> Q2,
                                                       LdqbdOptions options) {
        Matrix[] R = new Matrix[N];
        Matrix Q0_Nm1 = Q0.get(N - 1);
        Matrix Q1_N = Q1.get(N);
        Matrix U_N = Q1_N.scale(-1.0);
        R[N - 1] = safeMatrixDivide(Q0_Nm1, U_N);

        for (int n = N - 1; n >= 1; n--) {
            Matrix Q0_nm1 = Q0.get(n - 1);
            Matrix Q1_n = Q1.get(n);
            Matrix Q2_np1 = Q2.get(n);
            Matrix R_np1 = R[n];

            Matrix RQ2 = R_np1.mult(Q2_np1);
            Matrix U = Q1_n.scale(-1.0).sub(RQ2);
            R[n - 1] = safeMatrixDivide(Q0_nm1, U);
        }

        List<Matrix> result = new ArrayList<Matrix>(N);
        for (int i = 0; i < N; i++) result.add(R[i]);
        return result;
    }

    private static Matrix safeMatrixDivide(Matrix A, Matrix U) {
        if (U.length() == 1) {
            double uVal = U.get(0, 0);
            if (Math.abs(uVal) > 1e-14) {
                return A.scale(1.0 / uVal);
            } else {
                return new Matrix(A.getNumRows(), A.getNumCols());
            }
        } else {
            double det = U.det();
            if (Math.abs(det) > 1e-14) {
                return A.mult(U.inv());
            } else {
                return A.mult(pinv(U));
            }
        }
    }

    private static Matrix pinv(Matrix A) {
        jline.io.Ret.SVD svd = A.svd();
        Matrix U = svd.u;
        Matrix S = svd.s;
        Matrix V = svd.v;

        int m = A.getNumRows();
        int n = A.getNumCols();
        Matrix Splus = new Matrix(n, m);
        double tol = 1e-10 * Math.max(m, n) * (S.length() > 0 ? S.get(0, 0) : 0.0);
        int rank = S.getNumRows();
        for (int i = 0; i < rank; i++) {
            double sigma = S.get(i, 0);
            if (Math.abs(sigma) > tol) {
                Splus.set(i, i, 1.0 / sigma);
            }
        }
        return V.mult(Splus).mult(U.transpose());
    }

    private static Matrix computeStationaryDist(List<Matrix> R, List<Matrix> Q0, List<Matrix> Q1, List<Matrix> Q2,
                                                int N, LdqbdOptions options, List<Matrix> piCellsOut) {
        boolean isScalar = true;
        for (Matrix m : Q1) {
            if (m.length() != 1) { isScalar = false; break; }
        }
        Set<Integer> dimSet = new HashSet<Integer>();
        for (Matrix m : Q1) dimSet.add(m.getNumRows());
        boolean isHomogeneous = dimSet.size() == 1;

        if (isScalar && isHomogeneous) {
            Matrix pi = new Matrix(1, N + 1);
            pi.set(0, 0, 1.0);
            for (int n = 1; n <= N; n++) {
                if (R.get(n - 1).length() == 1) {
                    pi.set(0, n, pi.get(0, n - 1) * R.get(n - 1).get(0, 0));
                } else {
                    pi.set(0, n, 0.0);
                }
            }
            double sum = pi.sumRows(0);
            if (sum > 0) pi.scaleEq(1.0 / sum);
            for (int n = 0; n <= N; n++) {
                Matrix cell = new Matrix(1, 1);
                cell.set(0, 0, pi.get(0, n));
                piCellsOut.add(cell);
            }
            return pi;
        } else {
            List<Matrix> piCells = new ArrayList<Matrix>(N + 1);
            for (int i = 0; i <= N; i++) piCells.add(new Matrix(1, 1));

            if (Q1.get(0).length() == 1) {
                Matrix pi0 = new Matrix(1, 1);
                pi0.set(0, 0, 1.0);
                piCells.set(0, pi0);
            } else {
                Matrix Q1_0 = Q1.get(0);
                Matrix Q2_1 = Q2.get(0);
                Matrix A = Q1_0.add(1.0, R.get(0).mult(Q2_1));
                Matrix pi0 = solveLeftNullSpace(A);
                piCells.set(0, pi0);
            }

            for (int n = 1; n <= N; n++) {
                piCells.set(n, piCells.get(n - 1).mult(R.get(n - 1)));
            }

            double total = 0.0;
            for (int n = 0; n <= N; n++) total += piCells.get(n).sumRows(0);

            Matrix pi = new Matrix(1, N + 1);
            for (int n = 0; n <= N; n++) {
                if (total > 0) piCells.get(n).scaleEq(1.0 / total);
                pi.set(0, n, piCells.get(n).sumRows(0));
            }
            piCellsOut.addAll(piCells);
            return pi;
        }
    }

    /**
     * Boundary vector of the level-0 block: the solution of
     * pi_0 * (Q1^(0) + R^(1) Q2^(1)) = 0.
     *
     * That matrix is the generator of the process censored on level 0, so its
     * stationary distribution IS the boundary vector and CTMC_SOLVE is the right
     * instrument. A power iteration on A' is not: it converges to the DOMINANT
     * left direction of A, which is not the null one, and returned a plausible
     * but wrong vector for every block bigger than 1x1. No caller reached this
     * branch before the bgchain method, whose level 0 carries the arrival and
     * environment phases.
     */
    private static Matrix solveLeftNullSpace(Matrix A) {
        return jline.api.mc.Ctmc_solve.ctmc_solve(jline.api.mc.Ctmc_makeinfgen.ctmc_makeinfgen(A));
    }
}
