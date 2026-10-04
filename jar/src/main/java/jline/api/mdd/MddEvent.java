package jline.api.mdd;

/**
 * One event of the Kronecker rate descriptor.
 *
 * <p>The rate matrix is R = sum_e (kron_k W_k^e) restricted to the reachable
 * set. An event touches only the levels it names; every other level carries the
 * identity, which {@link Mdd_mcd} supplies rather than storing.</p>
 */
public class MddEvent {

    /** Station (or place) the event departs from, 0-based. */
    public final int a;
    /** Station (or place) the event arrives at, 0-based; equals a for an internal event. */
    public final int b;
    /** Levels the event touches, as 0-based station indices, aligned with W. */
    public final int[] lev;
    /** Local matrices at the levels named by lev. */
    public final MddLocalMatrix[] W;

    public MddEvent(int a, int b, int[] lev, MddLocalMatrix[] W) {
        this.a = a;
        this.b = b;
        this.lev = lev;
        this.W = W;
    }
}
