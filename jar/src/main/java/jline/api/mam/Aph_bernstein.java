/**
 * @file Absorbing Phase-type distribution Bernstein polynomial approximation
 *
 * Constructs APH distributions using Bernstein exponential approximation methods.
 * Advanced technique for approximating arbitrary density functions with phase-type distributions.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.function.DoubleUnaryOperator;

import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class Aph_bernstein {
    private Aph_bernstein() {}

    /**
     * Fits an Acyclic Phase-type distribution using Bernstein's approximation.
     *
     * Based on: Andras Horvath, Enrico Vicario: Construction of Phase Type Distributions
     * by Bernstein Exponentials. EPEW 2023: 201-215.
     *
     * @param f Function for the density f(x) to approximate
     * @param order Approximation order
     * @return Pair containing (D0, D1) matrices representing the APH distribution
     */
    public static Pair<Matrix, Matrix> aph_bernstein(DoubleUnaryOperator f, int order) {
        int n = order;

        // Compute normalization constant c
        double c = 0.0;
        for (int i = 1; i <= n; i++) {
            c += f.applyAsDouble(-Math.log((double) i / n)) / i;
        }

        // Create the generator matrix T (diagonal with -[1:n] and super-diagonal with [1:(n-1)])
        Matrix T = Matrix.zeros(n, n);
        for (int i = 0; i < n; i++) {
            T.set(i, i, -(double) (i + 1));
            if (i < n - 1) {
                T.set(i, i + 1, (double) (i + 1));
            }
        }

        // Compute initial probability vector alpha
        Matrix alpha = Matrix.zeros(1, n);
        for (int i = 1; i <= n; i++) {
            alpha.set(0, i - 1, f.applyAsDouble(-Math.log((double) i / n)) / (i * c));
        }

        // Create P matrix (replicated alpha vector)
        Matrix P = Matrix.zeros(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                P.set(i, j, alpha.get(0, j));
            }
        }

        // Compute D0 and D1 matrices
        Matrix D0 = T;
        Matrix D1 = T.scale(-1.0).mult(P);

        return new Pair<Matrix, Matrix>(D0, D1);
    }
}
