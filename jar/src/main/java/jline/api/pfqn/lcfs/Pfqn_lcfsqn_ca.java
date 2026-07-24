/**
 * Convolution Algorithm for LCFS Queueing Networks
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.lcfs;

import jline.util.PopulationLattice;
import jline.util.matrix.Matrix;

public final class Pfqn_lcfsqn_ca {
    private Pfqn_lcfsqn_ca() {}

    /**
     * Convolution algorithm for multiclass LCFS queueing networks.
     */
    public static LcfsqnCaResult pfqn_lcfsqn_ca(Matrix alpha, Matrix beta, Matrix N) {
        int R = alpha.length();
        Matrix populationVector = (N != null) ? N : Matrix.ones(1, R);

        int K = (int) populationVector.elementSum();

        if (K == 0) {
            return new LcfsqnCaResult(1.0, 1.0);
        }

        Matrix prods = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            double prod = 1.0;
            for (int j = 0; j < r; j++) {
                prod *= (populationVector.get(j) + 1);
            }
            prods.set(r, prod);
        }

        int totalSize = 1;
        for (int r = 0; r < R; r++) {
            totalSize *= ((int) populationVector.get(r) + 1);
        }

        double[] G = new double[totalSize];
        double[] V = new double[totalSize];

        Matrix n = PopulationLattice.pprod(populationVector);

        while (n.get(0) >= 0) {
            int idx = hashpopLcfs(n, populationVector, R, prods);
            int sumN = (int) n.elementSum();

            if (sumN == 0) {
                G[idx] = 1.0;
                V[idx] = 1.0;
            } else {
                V[idx] = 0.0;
                G[idx] = 0.0;

                for (int r = 0; r < R; r++) {
                    if (n.get(r) > 0) {
                        Matrix n_minus_r = oner(n, r);
                        int idx_r = hashpopLcfs(n_minus_r, populationVector, R, prods);

                        V[idx] = V[idx] + V[idx_r];

                        G[idx] = G[idx] + Math.pow(alpha.get(r), (double) (sumN - 1)) * beta.get(r) * G[idx_r];
                    }
                }

                double prodAlphaN = 1.0;
                for (int r = 0; r < R; r++) {
                    prodAlphaN *= Math.pow(alpha.get(r), n.get(r));
                }
                V[idx] = prodAlphaN * V[idx];
                G[idx] = G[idx] + V[idx];
            }

            n = PopulationLattice.pprod(n, populationVector);
        }

        return new LcfsqnCaResult(G[totalSize - 1], V[totalSize - 1]);
    }

    public static LcfsqnCaResult pfqn_lcfsqn_ca(Matrix alpha, Matrix beta) {
        return pfqn_lcfsqn_ca(alpha, beta, null);
    }

    /**
     * Computes hash index for LCFS population lattice.
     */
    private static int hashpopLcfs(Matrix n, Matrix N, int R, Matrix prods) {
        int idx = 0;
        for (int r = 0; r < R; r++) {
            idx += (int) (prods.get(r) * n.get(r));
        }
        return idx;
    }

    /**
     * Decrements the r-th component of vector n by 1.
     */
    private static Matrix oner(Matrix n, int r) {
        Matrix result = n.copy();
        result.set(r, result.get(r) - 1);
        return result;
    }
}
