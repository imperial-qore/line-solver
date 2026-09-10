package jline.api.mapqn;

import jline.util.matrix.Matrix;

/**
 * Parameters for MVA version models.
 */
public class MVAVersionParameters extends Mapqn_parameters {
    public final int M;
    public final int N;
    public final int K;
    public final double[] muM;
    public final Matrix muMAP;
    public final Matrix r;
    public final Matrix v;

    public MVAVersionParameters(int M, int N, int K, double[] muM, Matrix muMAP, Matrix r, Matrix v) {
        this.M = M;
        this.N = N;
        this.K = K;
        this.muM = muM;
        this.muMAP = muMAP;
        this.r = r;
        this.v = v;
    }

    @Override public int getM() { return M; }
    @Override public int getN() { return N; }

    @Override
    public void validate() {
        super.validate();
        if (K <= 0) throw new IllegalArgumentException("K must be positive");
        if (muM.length != M - 1) throw new IllegalArgumentException("muM array size must equal M-1");
        if (muMAP.getNumRows() != K || muMAP.getNumCols() != K) throw new IllegalArgumentException("muMAP must be KxK matrix");
        if (r.getNumRows() != M || r.getNumCols() != M) throw new IllegalArgumentException("r must be MxM matrix");
        if (v.getNumRows() != K || v.getNumCols() != K) throw new IllegalArgumentException("v must be KxK matrix");

        for (double mv : muM) {
            if (mv < 0) throw new IllegalArgumentException("Service rates must be non-negative");
        }
        for (int row = 0; row < muMAP.getNumRows(); row++) {
            for (int col = 0; col < muMAP.getNumCols(); col++) {
                if (muMAP.get(row, col) < 0) throw new IllegalArgumentException("MAP service rates must be non-negative");
            }
        }
        for (int row = 0; row < r.getNumRows(); row++) {
            for (int col = 0; col < r.getNumCols(); col++) {
                if (r.get(row, col) < 0) throw new IllegalArgumentException("Routing probabilities must be non-negative");
            }
        }
        for (int row = 0; row < v.getNumRows(); row++) {
            for (int col = 0; col < v.getNumCols(); col++) {
                if (v.get(row, col) < 0) throw new IllegalArgumentException("Level change rates must be non-negative");
            }
        }
    }

    public double q(int i, int j, int k, int h) {
        if (i < M - 1) {
            return (k == h) ? r.get(i, j) * muM[i] : 0.0;
        } else {
            if (j < M - 1) {
                return r.get(M - 1, j) * muMAP.get(k, h);
            } else {
                if (k != h) {
                    return v.get(k, h) + r.get(M - 1, M - 1) * muMAP.get(k, h);
                } else {
                    return 0.0;
                }
            }
        }
    }
}
