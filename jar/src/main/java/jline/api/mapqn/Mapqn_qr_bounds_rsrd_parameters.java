package jline.api.mapqn;

import jline.util.matrix.Matrix;

/**
 * Parameters for QR Bounds RSRD (Repetitive Service Random Destination) model.
 */
public class Mapqn_qr_bounds_rsrd_parameters extends Mapqn_parameters {
    public final int M;
    public final int N;
    public final int[] F;
    public final int[] K;
    public final Matrix[] mu;
    public final Matrix[] v;
    public final double[][] alpha;
    public final Matrix r;

    public Mapqn_qr_bounds_rsrd_parameters(int M, int N, int[] F, int[] K, Matrix[] mu,
                                            Matrix[] v, double[][] alpha, Matrix r) {
        this.M = M;
        this.N = N;
        this.F = F;
        this.K = K;
        this.mu = mu;
        this.v = v;
        this.alpha = alpha;
        this.r = r;
    }

    @Override public int getM() { return M; }
    @Override public int getN() { return N; }

    @Override
    public void validate() {
        super.validate();
        if (F.length != M) throw new IllegalArgumentException("F array size must equal M");
        if (K.length != M) throw new IllegalArgumentException("K array size must equal M");
        if (mu.length != M) throw new IllegalArgumentException("mu array size must equal M");
        if (v.length != M) throw new IllegalArgumentException("v array size must equal M");
        if (alpha.length != M) throw new IllegalArgumentException("alpha array size must equal M");
        if (r.getNumRows() != M || r.getNumCols() != M) throw new IllegalArgumentException("r must be MxM matrix");

        for (int i = 0; i < M; i++) {
            if (F[i] <= 0) throw new IllegalArgumentException("F[" + i + "] must be positive");
            if (K[i] <= 0) throw new IllegalArgumentException("K[" + i + "] must be positive");
            if (alpha[i].length != N) throw new IllegalArgumentException("alpha[" + i + "] array size must equal N");
            if (mu[i].getNumRows() != K[i] || mu[i].getNumCols() != K[i])
                throw new IllegalArgumentException("mu[" + i + "] must be " + K[i] + "x" + K[i] + " matrix");
            if (v[i].getNumRows() != K[i] || v[i].getNumCols() != K[i])
                throw new IllegalArgumentException("v[" + i + "] must be " + K[i] + "x" + K[i] + " matrix");
            for (double a : alpha[i]) {
                if (a < 0) throw new IllegalArgumentException("alpha[" + i + "] values must be non-negative");
            }
        }
    }

    public double q(int i, int j, int k, int h, int n) {
        if (n == 0) return 0.0;
        if (j != i) {
            return r.get(i, j) * mu[i].get(k, h) * alpha[i][n - 1];
        } else {
            return v[i].get(k, h) * alpha[i][n - 1] + r.get(i, i) * mu[i].get(k, h) * alpha[i][n - 1];
        }
    }
}
