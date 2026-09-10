/**
 * @file Quadrature of v*(-T0)^-1 by uniformized Euler integration
 *
 * @since LINE 3.0
 */
package jline.api.fes;

import jline.util.matrix.Matrix;

/**
 * Integrates v*int_0^inf exp(T0 t) dt by the trapezoid rule.
 *
 * Evaluates the product v*(-T0)^-1 without any factorization, as done in Section 5.2.2 of
 * Casale, Mi, Cherkasova and Smirni, IEEE Trans. Soft. Eng. 37(5), 2011. The propagated
 * vector uses the Euler approximation exp(T0 dt) ~ I + T0 dt, so only vector-matrix
 * products are performed and the sparsity of T0 is preserved throughout.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public final class Fes_map_euler {
    private Fes_map_euler() {}

    /**
     * Approximates v*(-T0)^-1 by quadrature.
     *
     * @param v       row vector to be multiplied by (-T0)^-1
     * @param T0      hidden transitions of the MAP, a stable matrix
     * @param dt      integration step, below 1/max(abs(diag(T0)))
     * @param tol     relative mass left when the integration stops
     * @param iterMax maximum number of integration steps
     * @return row vector approximating v*(-T0)^-1
     */
    public static Matrix fes_map_euler(Matrix v, Matrix T0, double dt, double tol, int iterMax) {
        Matrix y = new Matrix(v.getNumRows(), v.getNumCols());
        Matrix z = v.copy();
        double nrm0 = norm1(v);

        for (int it = 0; it < iterMax; it++) {
            Matrix znext = z.add(dt, z.mult(T0));
            y = y.add(dt / 2, z).add(dt / 2, znext);
            z = znext;
            if (norm1(z) <= tol * nrm0) {
                break;
            }
        }
        return y;
    }

    private static double norm1(Matrix v) {
        double s = 0;
        for (int i = 0; i < v.getNumCols(); i++) {
            s += Math.abs(v.get(0, i));
        }
        return s;
    }
}
