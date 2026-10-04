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
    /**
     * Converged mean-field state, item-major over lists 0..h as a column of
     * n*(h+1) entries, or null when the solver produced no list-resolved state.
     * List 0 is outside the cache, so pi0 above is its first n entries; the hit
     * probability of an item RESOLVED BY LIST reads lists 1..h.
     */
    public final Matrix xss;

    public CacheMissFpiResult(double M, Matrix MU, Matrix MI, Matrix pi0, Matrix xss) {
        this.M = M;
        this.MU = MU;
        this.MI = MI;
        this.pi0 = pi0;
        this.xss = xss;
    }

    public CacheMissFpiResult(double M, Matrix MU, Matrix MI, Matrix pi0) {
        this(M, MU, MI, pi0, null);
    }

    public CacheMissFpiResult(double M) {
        this(M, null, null, null, null);
    }

    public double getM() { return M; }
    public Matrix getMU() { return MU; }
    public Matrix getMI() { return MI; }
    public Matrix getPi0() { return pi0; }
    public Matrix getXss() { return xss; }
}
