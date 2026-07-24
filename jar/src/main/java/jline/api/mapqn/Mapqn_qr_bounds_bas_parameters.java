package jline.api.mapqn;

import jline.util.matrix.Matrix;

/**
 * Parameters for QR Bounds BAS (Blocking After Service) model.
 */
public class Mapqn_qr_bounds_bas_parameters extends Mapqn_parameters {
    public final int M;
    public final int N;
    public final int MR;
    public final int f;
    public final int[] K;
    public final int[] F;
    public final Matrix MM;
    public final Matrix MM1;
    public final int[] ZZ;
    public final Matrix BB;
    public final Matrix[] mu;
    public final Matrix[] v;
    public final Matrix r;
    public final int ZM;

    public Mapqn_qr_bounds_bas_parameters(int M, int N, int MR, int f, int[] K, int[] F,
                                           Matrix MM, Matrix MM1, int[] ZZ, Matrix BB,
                                           Matrix[] mu, Matrix[] v, Matrix r) {
        this.M = M;
        this.N = N;
        this.MR = MR;
        this.f = f;
        this.K = K;
        this.F = F;
        this.MM = MM;
        this.MM1 = MM1;
        this.ZZ = ZZ;
        this.BB = BB;
        this.mu = mu;
        this.v = v;
        this.r = r;
        int max = 0;
        for (int z : ZZ) if (z > max) max = z;
        this.ZM = max;
    }

    @Override public int getM() { return M; }
    @Override public int getN() { return N; }

    @Override
    public void validate() {
        super.validate();
        if (MR <= 0) throw new IllegalArgumentException("MR must be positive");
        if (f < 1 || f > M) throw new IllegalArgumentException("f must be between 1 and M");
        if (K.length != M) throw new IllegalArgumentException("K array size must equal M");
        if (F.length != M) throw new IllegalArgumentException("F array size must equal M");
        if (MM.getNumRows() != MR || MM.getNumCols() != 2)
            throw new IllegalArgumentException("MM must be MRx2 matrix (blocking order pairs)");
        if (MM1.getNumRows() != MR || MM1.getNumCols() != M)
            throw new IllegalArgumentException("MM1 must be MRxM matrix");
        if (ZZ.length != MR) throw new IllegalArgumentException("ZZ array size must equal MR");
        if (BB.getNumRows() != MR || BB.getNumCols() != M) throw new IllegalArgumentException("BB must be MRxM matrix");
        if (mu.length != M) throw new IllegalArgumentException("mu array size must equal M");
        if (v.length != M) throw new IllegalArgumentException("v array size must equal M");
        if (r.getNumRows() != M || r.getNumCols() != M) throw new IllegalArgumentException("r must be MxM matrix");
    }

    public double q(int i, int j, int k, int h) {
        if (j != i) {
            return r.get(i, j) * mu[i].get(k, h);
        } else {
            return v[i].get(k, h) + r.get(i, i) * mu[i].get(k, h);
        }
    }

    public Matrix[][] q() {
        Matrix[][] result = new Matrix[M][M];
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                Matrix mat = new Matrix(K[i], K[i]);
                for (int ki = 0; ki < K[i]; ki++) {
                    for (int hi = 0; hi < K[i]; hi++) {
                        mat.set(ki, hi, q(i, j, ki, hi));
                    }
                }
                result[i][j] = mat;
            }
        }
        return result;
    }
}
