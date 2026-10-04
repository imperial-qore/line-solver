package jline.api.cache;

import jline.util.matrix.Matrix;

/**
 * Result of cache_rrm_meanfield: the terminal occupancy of the RANDOM(m)
 * mean field, the global miss rate and the miss ratio.
 */
public final class CacheRrmMeanfieldResult {
    /** (n x (1+h)) occupancy; column 0 is the probability of being out of cache. */
    public final Matrix x;
    /** lambda' * x(:,0). */
    public final double missrate;
    /** missrate / sum(lambda). */
    public final double missratio;

    public CacheRrmMeanfieldResult(Matrix x, double missrate, double missratio) {
        this.x = x;
        this.missrate = missrate;
        this.missratio = missratio;
    }

    public Matrix getX() { return x; }
    public double getMissrate() { return missrate; }
    public double getMissratio() { return missratio; }
}
