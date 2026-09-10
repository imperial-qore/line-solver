/**
 * Normalizing Constant for LCFS Queueing Networks
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.lcfs;

import jline.lib.perm.Permanent;
import jline.util.matrix.Matrix;

public final class Pfqn_lcfsqn_nc {
    private Pfqn_lcfsqn_nc() {}

    /**
     * Result for LCFS NC algorithm.
     */
    public static final class LcfsqnNcResult {
        public final double G;
        public final Matrix[] Ax;

        public LcfsqnNcResult(double G, Matrix[] Ax) {
            this.G = G;
            this.Ax = Ax;
        }
    }

    /**
     * Normalizing constant for multiclass LCFS queueing networks.
     */
    public static LcfsqnNcResult pfqn_lcfsqn_nc(Matrix alpha, Matrix beta, Matrix N) {
        int K = (int) N.elementSum();
        int R = N.length();
        double G = 0.0;

        Matrix[] Ax = new Matrix[K + 1];
        for (int x = 0; x <= K; x++) {
            Ax[x] = makeA(alpha, beta, x, K, R);
        }

        for (int x = 0; x <= K; x++) {
            G += perm(Ax[x], N);
        }

        // see _kb/03-api-layer.md for rationale
        double normalizer = 1.0;
        for (int r = 0; r < R; r++) {
            int Nr = (int) N.get(r);
            for (int k = 2; k <= Nr; k++) {
                normalizer *= k;
            }
        }
        G = G / normalizer;

        return new LcfsqnNcResult(G, Ax);
    }

    private static Matrix makeA(Matrix alpha, Matrix beta, int x, int K, int R) {
        Matrix A = new Matrix(R, K);
        for (int i = 0; i < R; i++) {
            for (int j = 0; j < x; j++) {
                A.set(i, j, Math.pow(alpha.get(i), (double) (j + 1)));
            }
            for (int j = 0; j < (K - x); j++) {
                A.set(i, x + j, Math.pow(alpha.get(i), (double) (x + j)) * beta.get(i));
            }
        }
        return A;
    }

    private static double perm(Matrix A, Matrix N) {
        int K = (int) N.elementSum();
        int R = N.length();

        Matrix expandedMatrix = new Matrix(K, K);
        int rowIdx = 0;
        for (int r = 0; r < R; r++) {
            int count = (int) N.get(r);
            for (int rep = 0; rep < count; rep++) {
                for (int col = 0; col < K; col++) {
                    expandedMatrix.set(rowIdx, col, A.get(r, col));
                }
                rowIdx++;
            }
        }

        Permanent permSolver = new Permanent(expandedMatrix, true);
        return permSolver.getValue();
    }
}
