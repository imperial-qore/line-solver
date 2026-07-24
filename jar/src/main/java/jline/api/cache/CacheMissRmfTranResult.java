/**
 * @file Transient RMF cache-miss result.
 *
 * @since LINE 3.0
 */
package jline.api.cache;

/**
 * Transient refined-mean-field cache trajectory returned by
 * {@link Cache_miss_rmf#cache_miss_rmf_tran}.
 *
 * <p>Java port of the transient outputs of the MATLAB {@code cache_miss_rmf.m}
 * tspan path: the time grid, the per-item list-0 occupancy trajectory, the
 * per-user miss-rate trajectory, and the full DDPP occupancy trajectory used to
 * carry the cache mean occupancy across environment switches.</p>
 */
public final class CacheMissRmfTranResult {
    /** Time grid (nt,). */
    public final double[] t;
    /** Per-item list-0 occupancy (miss probability) trajectory (n x nt). */
    public final double[][] pi0_t;
    /** Per-user miss-rate trajectory (u x nt). */
    public final double[][] MU_t;
    /** Full DDPP occupancy trajectory (dim x nt). */
    public final double[][] xtraj;

    public CacheMissRmfTranResult(double[] t, double[][] pi0_t, double[][] MU_t, double[][] xtraj) {
        this.t = t;
        this.pi0_t = pi0_t;
        this.MU_t = MU_t;
        this.xtraj = xtraj;
    }
}
