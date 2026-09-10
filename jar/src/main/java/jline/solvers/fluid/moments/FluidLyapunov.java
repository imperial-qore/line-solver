/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.moments;

import jline.util.matrix.Matrix;
import jline.io.Ret;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.util.FastMath;

import java.util.List;

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.mfilename;

/**
 * Stationary covariance of the linear noise approximation. Java twin of the
 * MATLAB {@code fluid_lyapunov}.
 *
 * <p>Around a fixed point x* of the fluid drift, the fluctuation process
 * Z = (X - x*) obeys the linear stochastic differential equation
 * {@code dZ = A*Z*dt + sqrt(Qdiff)*dW}, whose stationary covariance solves the
 * Lyapunov equation</p>
 *
 * <pre>
 *   A*Sigma + Sigma*A' + Qdiff = 0,   Qdiff = D*diag(r(x*))*D'
 * </pre>
 *
 * <p>A is singular whenever the model conserves population: every closed class
 * contributes a left null vector, so the equation has no unique solution on the
 * full state space. It does have one on the reachable subspace, which is exactly
 * range(D): the state can only move along jump directions, so the fluctuation
 * lives there and nowhere else. Both {@code A = D*diag(rateBase)*G} and Qdiff
 * map into range(D) as well, so restricting to an orthonormal basis V of
 * range(D) is an exact reduction, not an approximation, and the reduced Lyapunov
 * equation is nonsingular whenever the fixed point is stable.</p>
 *
 * @see FluidMomentTerms
 */
public final class FluidLyapunov {

    private static final double SQRT_EPS = FastMath.sqrt(2.220446049250313e-16);

    private FluidLyapunov() {
    }

    /** The stationary covariance together with the diagnostics of the reduced solve. */
    public static final class Result {
        public final Matrix sigma;
        public final int rank;
        public final double maxRealEig;
        public final boolean stable;

        public Result(Matrix sigma, int rank, double maxRealEig, boolean stable) {
            this.sigma = sigma;
            this.rank = rank;
            this.maxRealEig = maxRealEig;
            this.stable = stable;
        }
    }

    /**
     * Solves the Lyapunov equation on the reachable subspace.
     *
     * @param A     drift Jacobian at the fixed point (n x n)
     * @param Qdiff diffusion matrix D*diag(r)*D' (n x n)
     * @param D     jump matrix (n x nevents), spanning the reachable subspace
     * @param tol   stability margin; eigenvalues of the reduced A with real part
     *              above -tol are reported as non-hyperbolic
     * @return the stationary covariance, supported on range(D)
     */
    public static Result solve(Matrix A, Matrix Qdiff, Matrix D, double tol) {
        if (!(tol > 0)) {
            tol = SQRT_EPS;
        }
        int n = A.getNumRows();
        Matrix V = orth(D);
        if (V.getNumCols() == 0) {
            return new Result(new Matrix(n, n), 0, Double.NEGATIVE_INFINITY, true);
        }

        Matrix Vt = V.transpose();
        Matrix Ar = Vt.mult(A).mult(V);
        Matrix Qr = Vt.mult(Qdiff).mult(V);
        symmetrise(Qr);

        double maxRe = Double.NEGATIVE_INFINITY;
        List<Complex> ev = Ar.eig();
        for (int i = 0; i < ev.size(); i++) {
            maxRe = FastMath.max(maxRe, ev.get(i).getReal());
        }
        boolean stable = maxRe < -tol;
        if (!stable) {
            throw new FluidNonHyperbolicException(String.format(
                    "The fluid fixed point is not exponentially stable on the reachable subspace "
                    + "(largest Jacobian eigenvalue has real part %g), so the linear noise approximation has no "
                    + "stationary covariance. This happens at an unstable model or at a drift kink; use "
                    + "options.method=\"closing\" for the mean only.", maxRe));
        }

        // Bartels-Stewart: Matrix.sylv(A,B,C) solves A*X + X*B = -C
        Matrix W = Matrix.sylv(Ar, Ar.transpose(), Qr);
        symmetrise(W);

        Matrix sigma = V.mult(W).mult(Vt);
        symmetrise(sigma);
        return new Result(sigma, V.getNumCols(), maxRe, stable);
    }

    /**
     * Orthonormal basis of the column space, mirroring MATLAB {@code orth}: the
     * left singular vectors whose singular value exceeds
     * {@code max(size(A))*eps*max(s)}.
     *
     * @param A any matrix
     * @return an (rows x rank) matrix with orthonormal columns
     */
    public static Matrix orth(Matrix A) {
        if (A.getNumRows() == 0 || A.getNumCols() == 0) {
            return new Matrix(A.getNumRows(), 0);
        }
        Ret.SVD svd = A.svd();
        Matrix U = svd.u;
        Matrix S = svd.s;
        double smax = 0;
        for (int i = 0; i < S.getNumRows(); i++) {
            smax = FastMath.max(smax, S.get(i, 0));
        }
        double tol = FastMath.max(A.getNumRows(), A.getNumCols()) * 2.220446049250313e-16 * smax;
        int rank = 0;
        for (int i = 0; i < S.getNumRows() && i < U.getNumCols(); i++) {
            if (S.get(i, 0) > tol) {
                rank++;
            }
        }
        Matrix V = new Matrix(U.getNumRows(), rank);
        for (int i = 0; i < U.getNumRows(); i++) {
            for (int j = 0; j < rank; j++) {
                V.set(i, j, U.get(i, j));
            }
        }
        return V;
    }

    private static void symmetrise(Matrix m) {
        int n = m.getNumRows();
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double v = 0.5 * (m.get(i, j) + m.get(j, i));
                m.set(i, j, v);
                m.set(j, i, v);
            }
        }
    }
}
