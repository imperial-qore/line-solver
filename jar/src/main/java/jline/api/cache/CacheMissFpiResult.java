package jline.api.cache;

import jline.util.matrix.Matrix;

/**
 * Result of cache_miss_fpi: global miss rate, per-user miss rates, per-item miss rates,
 * and per-item miss probabilities.
 */
public final class CacheMissFpiResult {
    public final double M;
    public final Matrix MU;
    public final Matrix MI;
    public final Matrix pi0;

    public CacheMissFpiResult(double M, Matrix MU, Matrix MI, Matrix pi0) {
        this.M = M;
        this.MU = MU;
        this.MI = MI;
        this.pi0 = pi0;
    }

    public CacheMissFpiResult(double M) {
        this(M, null, null, null);
    }

    public double getM() { return M; }
    public Matrix getMU() { return MU; }
    public Matrix getMI() { return MI; }
    public Matrix getPi0() { return pi0; }
}
