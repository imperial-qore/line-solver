/**
 * @file Sojourn-time distribution of a Markov-modulated fluid queue
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.lib.butools.queues.FluidQueueSTD;
import jline.util.matrix.Matrix;

public final class Mfq_sojourn {
    private Mfq_sojourn() {}

    /** Sojourn-time distribution as {alpha, A} (ME if transToPH=false, PH otherwise). */
    public static Matrix[] mfq_sojourn(Matrix Q, Matrix Rin, Matrix Rout, Matrix Q0, boolean transToPH) {
        return FluidQueueSTD.fluidQueueSTD(Q, Rin, Rout, Q0, transToPH);
    }
}
