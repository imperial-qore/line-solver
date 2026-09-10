package jline.solvers;

import jline.util.matrix.Matrix;

/**
 * Blocking parameters for the QRF (Quadratic Reduction Framework) bounds
 * {@code qrf.bas} and {@code qrf.rsrd}.
 *
 * An explicit OVERRIDE of what {@link jline.api.sn.SnToQrfBlocking} derives
 * from the model. Leaving it null is the normal path: the blocking structure
 * fixes every field below, so the analyzer builds them rather than demanding
 * them. What it must never do is INVENT them -- assuming no blocking (MR = 1,
 * unbounded buffers) answers a different model and lands roughly 31x farther
 * from exact. Matches MATLAB's {@code options.config.qrf_params} and the
 * native-Python {@code config['qrf_params']} contract.
 *
 * <p>{@code qrf.rsrd} reads NONE of these: RS-RD carries no blocking tables at
 * all, only F, which it takes from {@link jline.api.sn.SnToQrfCapacity}.
 *
 * @since LINE 3.0
 */
public class QrfParams {
    /** The ONE finite-capacity queue, 1-based. The formulation has no room for a second. */
    public int f = 1;
    /** Number of blocking configurations. Configuration 1 MUST be the empty one. */
    public int MR = 1;
    /** Maximum blocking depth, which must equal max(ZZ) and must be the REACHABLE maximum. */
    public int ZM = 0;
    /** (M) occupancy bound of each station: its buffer where that binds, else the population. */
    public int[] F;
    /** (MR) blocking depth of each configuration, i.e. how many queues it holds blocked. */
    public int[] ZZ;
    /** (MR x .) blocking order, 1-based; only column 0 is read, the head that takes f's freed slot. */
    public Matrix MM;
    /** (MR x M) successor map: MM1(m,j) is the configuration reached when queue j becomes blocked. */
    public Matrix MM1;
    /** (MR x M) blocking state: BB(m,i) is 1 iff queue i is blocked in configuration m. */
    public Matrix BB;

    public QrfParams() {
    }

    public QrfParams copy() {
        QrfParams c = new QrfParams();
        c.f = this.f;
        c.MR = this.MR;
        c.ZM = this.ZM;
        c.F = (this.F == null) ? null : this.F.clone();
        c.ZZ = (this.ZZ == null) ? null : this.ZZ.clone();
        c.MM = (this.MM == null) ? null : new Matrix(this.MM);
        c.MM1 = (this.MM1 == null) ? null : new Matrix(this.MM1);
        c.BB = (this.BB == null) ? null : new Matrix(this.BB);
        return c;
    }
}
