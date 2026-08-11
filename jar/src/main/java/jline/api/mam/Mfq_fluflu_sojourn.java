/**
 * @file Sojourn-time distribution of a fluid/fluid queue
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.lib.butools.queues.FluidQueueSTD;
import jline.util.matrix.Matrix;

public final class Mfq_fluflu_sojourn {
    private Mfq_fluflu_sojourn() {}

    /** Sojourn-time distribution of a fluid queue with fluid-modulated service, as {alpha, A}. */
    public static Matrix[] mfq_fluflu_sojourn(Matrix Qin, Matrix Rin, Matrix Qout, Matrix Rout, boolean srv0stop, boolean transToPH) {
        return FluidQueueSTD.fluFluSTD(Qin, Rin, Qout, Rout, srv0stop, transToPH);
    }
}
