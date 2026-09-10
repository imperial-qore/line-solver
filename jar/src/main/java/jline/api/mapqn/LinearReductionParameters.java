package jline.api.mapqn;

import jline.util.matrix.Matrix;

/**
 * Parameters for linear reduction models with phases.
 */
public class LinearReductionParameters extends Mapqn_parameters {
    public final int M;
    public final int N;
    public final int[] K;
    public final Matrix[] mu;
    public final Matrix r;
    public final Matrix[] v;

    public LinearReductionParameters(int M, int N, int[] K, Matrix[] mu, Matrix r, Matrix[] v) {
        this.M = M;
        this.N = N;
        this.K = K;
        this.mu = mu;
        this.r = r;
        this.v = v;
    }

    @Override public int getM() { return M; }
    @Override public int getN() { return N; }

    @Override
    public void validate() {
        super.validate();
        if (K.length != M) throw new IllegalArgumentException("K array size must equal M");
        if (mu.length != M) throw new IllegalArgumentException("mu array size must equal M");
        if (r.getNumRows() != M || r.getNumCols() != M) throw new IllegalArgumentException("r must be MxM matrix");
        if (v.length != M) throw new IllegalArgumentException("v array size must equal M");
        for (int i = 0; i < M; i++) {
            if (K[i] <= 0) throw new IllegalArgumentException("K[" + i + "] must be positive");
            if (mu[i].getNumRows() != K[i] || mu[i].getNumCols() != K[i])
                throw new IllegalArgumentException("mu[" + i + "] must be " + K[i] + "x" + K[i] + " matrix");
            if (v[i].getNumRows() != K[i] || v[i].getNumCols() != K[i])
                throw new IllegalArgumentException("v[" + i + "] must be " + K[i] + "x" + K[i] + " matrix");
            for (int row = 0; row < r.getNumRows(); row++) {
                for (int col = 0; col < r.getNumCols(); col++) {
                    if (r.get(row, col) < 0) throw new IllegalArgumentException("Routing probabilities must be non-negative");
                }
            }
            for (int row = 0; row < mu[i].getNumRows(); row++) {
                for (int col = 0; col < mu[i].getNumCols(); col++) {
                    if (mu[i].get(row, col) < 0) throw new IllegalArgumentException("Service rates must be non-negative");
                }
            }
            for (int row = 0; row < v[i].getNumRows(); row++) {
                for (int col = 0; col < v[i].getNumCols(); col++) {
                    if (v[i].get(row, col) < 0) throw new IllegalArgumentException("Background rates must be non-negative");
                }
            }
        }
    }

    public double q(int i, int j, int k, int h) {
        if (j != i) {
            return r.get(i, j) * mu[i].get(k, h);
        } else {
            return v[i].get(k, h) + r.get(i, i) * mu[i].get(k, h);
        }
    }
}
