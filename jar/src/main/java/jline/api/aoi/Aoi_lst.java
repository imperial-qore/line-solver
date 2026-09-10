/**
 * @file Laplace-Stieltjes Transform functions for AoI analysis
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

import jline.util.matrix.Matrix;

public final class Aoi_lst {
    private Aoi_lst() {}

    /**
     * LST for exponential distribution: H*(s) = mu / (mu + s)
     */
    public static LstFunction aoi_lst_exp(final double mu) {
        if (mu <= 0) {
            throw new IllegalArgumentException("Rate mu must be positive");
        }
        return new LstFunction() {
            @Override
            public double evaluate(double s) {
                return mu / (mu + s);
            }
        };
    }

    /**
     * LST for deterministic (constant) distribution: H*(s) = exp(-s * d)
     */
    public static LstFunction aoi_lst_det(final double d) {
        if (d <= 0) {
            throw new IllegalArgumentException("Constant d must be positive");
        }
        return new LstFunction() {
            @Override
            public double evaluate(double s) {
                return Math.exp(-s * d);
            }
        };
    }

    /**
     * LST for Erlang-k distribution: H*(s) = (mu / (mu + s))^k
     */
    public static LstFunction aoi_lst_erlang(final int k, final double mu) {
        if (k < 1) {
            throw new IllegalArgumentException("Shape k must be a positive integer");
        }
        if (mu <= 0) {
            throw new IllegalArgumentException("Rate mu must be positive");
        }
        return new LstFunction() {
            @Override
            public double evaluate(double s) {
                return Math.pow(mu / (mu + s), (double) k);
            }
        };
    }

    /**
     * LST for phase-type distribution: H*(s) = alpha * (s*I - T)^{-1} * (-T * e)
     */
    public static LstFunction aoi_lst_ph(final Matrix alpha, final Matrix T) {
        final int n = alpha.getNumCols();
        if (T.getNumRows() != n || T.getNumCols() != n) {
            throw new IllegalArgumentException("T must be " + n + " x " + n + " to match alpha");
        }
        Matrix ones = Matrix.ones(n, 1);
        final Matrix exitRates = T.mult(ones).scale(-1.0);

        return new LstFunction() {
            @Override
            public double evaluate(double s) {
                Matrix sIminusT = Matrix.eye(n).scale(s).sub(T);
                Matrix x = new Matrix(n, 1);
                Matrix.solveSafe(sIminusT, exitRates, x);
                return alpha.mult(x).get(0, 0);
            }
        };
    }
}
