/**
 * @file Continuous-time fluid priority queue (Horvath 2015)
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.List;

import jline.lib.butools.queues.FluidPrioQueue;
import jline.util.matrix.Matrix;

public final class Mfq_prio_queue {
    private Mfq_prio_queue() {}

    /** Performance measures of a continuous-time fluid priority queue (delegates to BUTools FluidPrioQueue). */
    public static List<double[]> mfq_prio_queue(Matrix Q, Matrix R, double d, int[] classes,
                                                double prec, int erlMaxOrder, Object... measures) {
        return FluidPrioQueue.fluidPrioQueue(Q, R, d, classes, prec, erlMaxOrder, measures);
    }
}
