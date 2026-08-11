/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.lib.fjcodes;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class ComputePi {
    private ComputePi() {}

    /**
     * Compute pseudo-inverse of a matrix using SVD.
     * Handles rank-deficient matrices by zeroing out small singular values.
     */
    private static Matrix pseudoInverse(Matrix A, double tol) {
        Ret.SVD svdResult = A.svd();
        Matrix U = svdResult.u;
        Matrix S = svdResult.s;
        Matrix V = svdResult.v;

        int m = A.getNumRows();
        int n = A.getNumCols();
        int k = Math.min(m, n);
        Matrix Splus = new Matrix(k, k);

        int numSingularValues = S.getNumRows();
        for (int i = 0; i < numSingularValues; i++) {
            double sigma = S.get(i, 0);
            if (Math.abs(sigma) > tol) {
                Splus.set(i, i, 1.0 / sigma);
            }
        }

        return V.mult(Splus).mult(U.transpose());
    }

    private static Matrix pseudoInverse(Matrix A) {
        return pseudoInverse(A, 1e-10);
    }

    /**
     * Solve overdetermined/rank-deficient system A * x = b using pseudo-inverse.
     */
    private static Matrix solveLeastSquares(Matrix A, Matrix b) {
        Matrix Aplus = pseudoInverse(A);
        return Aplus.mult(b);
    }

    /**
     * Compute steady-state distribution for the case with 2 replicas.
     *
     * Solves for the steady-state probability pi0 and expected number En1.
     *
     * @param T T-matrix from NARE solution
     * @param arrival Arrival process
     * @param services Service process for single subtask
     * @param service_h Service representation for 2-node job
     * @param C Capacity parameter
     * @param S State transition matrix
     * @param A_jump Jump matrix
     * @return PiResult with pi0 and En1
     */
    public static PiResult computePi(
            Matrix T,
            FJArrival arrival,
            FJService services,
            FJServiceH service_h,
            int C,
            Matrix S,
            Matrix A_jump) {

        if (services.getSerChoice() == 1) {
            // Exponential service case
            int ms = S.getNumCols();
            Matrix S_notallbusy = FJStateSpace.constructNotAllBusy(C, services, service_h);
            Matrix Q0 = FJUtils.kronsum(S_notallbusy, arrival.getLambda0());

            int da = arrival.getLambda0().getNumRows();

            // Solve Lyapunov equation: T*Igral + Igral*(kron(I, lambda0)) = -I
            Matrix B = Matrix.eye(ms).kron(arrival.getLambda0());
            Matrix C_lyap = Matrix.eye(da * ms).scale(-1.0);
            Matrix Igral = Matrix.lyap(T, B, C_lyap, null);

            // pi0mat = Igral * kron(A_jump, I) * inv(Q0) * kron(I, lambda1)
            Matrix Q0inv = pseudoInverse(Q0);
            Matrix pi0mat = Igral
                    .mult(A_jump.kron(Matrix.eye(da)))
                    .mult(Q0inv)
                    .mult(Matrix.eye(ms).kron(arrival.getLambda1()));

            // Solve for pi0: pi0 * (pi0mat - I) = 0, pi0 * 1 = 1
            Matrix rhsVec = new Matrix(1, ms * da + 1);
            for (int i = 0; i < ms * da; i++) {
                rhsVec.set(0, i, 0.0);
            }
            rhsVec.set(0, ms * da, -1.0);

            Matrix lhsMat = new Matrix(ms * da, ms * da + 1);

            for (int i = 0; i < ms * da; i++) {
                for (int j = 0; j < ms * da; j++) {
                    double value = (i == j) ? pi0mat.get(i, j) - 1.0 : pi0mat.get(i, j);
                    lhsMat.set(i, j, value);
                }
            }

            Matrix Tinv = pseudoInverse(T);
            for (int i = 0; i < ms * da; i++) {
                double rowSum = 0.0;
                for (int j = 0; j < ms * da; j++) {
                    rowSum += Tinv.get(i, j);
                }
                lhsMat.set(i, ms * da, rowSum);
            }

            Matrix pi0_T = solveLeastSquares(lhsMat.transpose(), rhsVec.transpose());
            Matrix pi0 = pi0_T.transpose();

            double sumPi0 = 0.0;
            for (int i = 0; i < ms * da; i++) {
                sumPi0 += pi0.get(0, i);
            }

            Matrix pi0_times_pi0mat = pi0.mult(pi0mat);
            double sumResult = 0.0;
            for (int i = 0; i < ms * da; i++) {
                sumResult += pi0_times_pi0mat.get(0, i);
            }

            double En1 = (1.0 / sumPi0) * sumResult;

            return new PiResult(pi0, En1);

        } else {
            // General PH service case
            SRKResult result = FJStateSpace.constructSRK(C, services, service_h, S);
            Matrix Se = result.Se;
            Matrix Sestar = result.Sestar;
            Matrix R0 = result.R0;
            Matrix Ke = result.Ke;
            Matrix Kc = result.Kc;

            int dtmat = T.getNumRows() / arrival.getLambda0().getNumCols();
            int dsexp = Se.getNumRows();

            // Extract submatrices
            Matrix Sedash = Matrix.getSubMatrix(Se, dtmat, dsexp, dtmat, dsexp);
            Matrix Rbusy = Matrix.getSubMatrix(R0, dtmat, dsexp, 0, R0.getNumCols());

            int da = arrival.getLambda0().length();
            int Iidle_rows = dsexp - dtmat;
            Matrix Iidle = new Matrix(Iidle_rows, da);
            for (int i = 0; i < Math.min(Iidle_rows, da); i++) {
                Iidle.set(i, i, 1.0);
            }

            Matrix Q0 = FJUtils.kronsum(Sedash, arrival.getLambda0());
            Matrix Qbusy = Rbusy.kron(arrival.getLambda1());

            Matrix Qidle = Q0;
            Matrix Qidle_inv = pseudoInverse(Qidle);
            Matrix Bmap = Iidle.mult(Qidle_inv).scale(-1.0).mult(Qbusy);

            Matrix Kemap = Ke.kron(Matrix.eye(da));
            Matrix Kcmap = Kc.kron(Matrix.eye(da));

            // Solve Lyapunov equation for Igral
            Matrix BB = Matrix.eye(Sestar.getNumCols()).kron(arrival.getLambda0());
            Matrix C_lyap = Kemap.mult(Sestar.kron(Matrix.eye(da)));
            Matrix Igral = Matrix.lyap(T, BB, C_lyap, null);

            Matrix pi0mat = Igral.mult(Bmap).mult(Kcmap);

            int ms_da = dtmat * da;

            Matrix rhsVec = new Matrix(1, ms_da + 1);
            for (int i = 0; i < ms_da; i++) {
                rhsVec.set(0, i, 0.0);
            }
            rhsVec.set(0, ms_da, -1.0);

            Matrix lhsMat = new Matrix(ms_da, ms_da + 1);

            for (int i = 0; i < ms_da; i++) {
                for (int j = 0; j < ms_da; j++) {
                    double value = (i == j) ? pi0mat.get(i, j) - 1.0 : pi0mat.get(i, j);
                    lhsMat.set(i, j, value);
                }
            }

            Matrix Tinv = pseudoInverse(T);
            for (int i = 0; i < ms_da; i++) {
                double rowSum = 0.0;
                for (int j = 0; j < ms_da; j++) {
                    rowSum += Tinv.get(i, j);
                }
                lhsMat.set(i, ms_da, rowSum);
            }

            Matrix pi0_T = solveLeastSquares(lhsMat.transpose(), rhsVec.transpose());
            Matrix pi0 = pi0_T.transpose();

            double sumPi0 = 0.0;
            for (int i = 0; i < ms_da; i++) {
                sumPi0 += pi0.get(0, i);
            }

            Matrix term = pi0.scale(-1.0)
                    .mult(Igral)
                    .mult(Iidle)
                    .mult(Qidle_inv)
                    .mult(Matrix.eye(Sedash.getNumRows()).kron(arrival.getLambda1()));

            double sumTerm = 0.0;
            for (int i = 0; i < term.getNumCols(); i++) {
                sumTerm += term.get(0, i);
            }

            double En1 = (1.0 / sumPi0) * sumTerm;

            return new PiResult(pi0, En1);
        }
    }
}
