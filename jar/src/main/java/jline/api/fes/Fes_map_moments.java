/**
 * @file Moments and index of dispersion of an inter-departure MAP
 *
 * @since LINE 3.0
 */
package jline.api.fes;

import jline.api.mc.Ctmc_solve;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Computes the first three moments, the lag-1 joint moment and the index of dispersion of
 * a MAP without forming (-T0)^-1.
 *
 * Evaluates equations (4), (5) and (7) of Casale, Mi, Cherkasova and Smirni, IEEE Trans.
 * Soft. Eng. 37(5), 2011. The inverse (-T0)^-1 is dense even when T0 is sparse, so it is
 * never formed: the moments follow from the vector recursion v_{k+1} = v_k (-T0)^-1, each
 * step being a linear solve. Method "euler" replaces the solve by the quadrature of
 * v*int_0^inf exp(T0 t) dt integrated by the trapezoid rule with the Euler approximation
 * exp(T0 dt) ~ I + T0 dt. Method "ssolve" is the default because it is exact and faster.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public final class Fes_map_moments {
    private Fes_map_moments() {}

    /** Relative mass left when the Euler quadrature stops. */
    public static final double TOL = 1e-12;
    /** Fraction of the uniformization bound used as integration step. */
    public static final double STEP_SAFETY = 0.1;
    /** Maximum number of Euler integration steps. */
    public static final int ITER_MAX = 1000000;

    /**
     * Computes the descriptors of an inter-departure MAP by sparse linear solves.
     *
     * @param MAP the pair (T0,T1)
     * @return the four descriptors and the index of dispersion
     */
    public static FesMapMomentsResult fes_map_moments(MatrixCell MAP) {
        return fes_map_moments(MAP.get(0), MAP.get(1), "ssolve");
    }

    /**
     * Computes the descriptors of an inter-departure MAP.
     *
     * @param T0     hidden transitions of the MAP
     * @param T1     marked transitions of the MAP
     * @param method "ssolve" for the linear solve, "euler" for the quadrature
     * @return the four descriptors and the index of dispersion
     */
    public static FesMapMomentsResult fes_map_moments(Matrix T0, Matrix T1, String method) {
        return fes_map_moments(T0, T1, method, STEP_SAFETY);
    }

    /**
     * Computes the descriptors of an inter-departure MAP.
     *
     * @param T0          hidden transitions of the MAP
     * @param T1          marked transitions of the MAP
     * @param method      "ssolve" for the linear solve, "euler" for the quadrature
     * @param stepSafety  fraction of the uniformization bound used as integration step
     * @return the four descriptors and the index of dispersion
     */
    public static FesMapMomentsResult fes_map_moments(Matrix T0, Matrix T1, String method, double stepSafety) {
        int dim = T0.getNumRows();
        Matrix e = new Matrix(dim, 1);
        for (int i = 0; i < dim; i++) {
            e.set(i, 0, 1.0);
        }

        Matrix Q = T0.add(1.0, T1);
        Matrix phi = Ctmc_solve.ctmc_solve(Q);
        Matrix pie = phi.mult(T1);
        double lambda = pie.elementSum();
        pie = pie.scale(1.0 / lambda);

        boolean euler = "euler".equalsIgnoreCase(method);
        double dt = 0;
        if (euler) {
            double dmax = 0;
            for (int i = 0; i < dim; i++) {
                dmax = Math.max(dmax, Math.abs(T0.get(i, i)));
            }
            dt = stepSafety / dmax;
        }
        Matrix negT0 = T0.scale(-1.0);
        Matrix negT0t = negT0.transpose();

        Matrix v1 = euler ? Fes_map_euler.fes_map_euler(pie, T0, dt, TOL, ITER_MAX) : solveRow(negT0t, pie);
        Matrix v2 = euler ? Fes_map_euler.fes_map_euler(v1, T0, dt, TOL, ITER_MAX) : solveRow(negT0t, v1);
        Matrix v3 = euler ? Fes_map_euler.fes_map_euler(v2, T0, dt, TOL, ITER_MAX) : solveRow(negT0t, v2);
        Matrix v2T1 = v2.mult(T1);
        Matrix v4 = euler ? Fes_map_euler.fes_map_euler(v2T1, T0, dt, TOL, ITER_MAX) : solveRow(negT0t, v2T1);

        double e1 = v1.elementSum();
        double e2 = 2 * v2.elementSum();
        double e3 = 6 * v3.elementSum();
        double e11 = v4.elementSum();

        // equation (7), with pie*inv(Q+e*phi) from the rank-one update y*Q = pie-phi
        // under the normalization y*e = 1
        Matrix A = Q.copy();
        for (int i = 0; i < dim; i++) {
            A.set(i, dim - 1, 1.0);
        }
        Matrix rhs = pie.add(-1.0, phi);
        rhs.set(0, dim - 1, 1.0);
        Matrix y = solveRow(A.transpose(), rhs);
        double idc = 1 + 2 * (lambda - y.mult(T1).mult(e).get(0, 0));

        return new FesMapMomentsResult(e1, e2, e3, e11, idc);
    }

    /** Solves the row system x*A = b given the transpose of A. */
    private static Matrix solveRow(Matrix At, Matrix b) {
        return At.leftMatrixDivide(b.transpose()).transpose();
    }
}
