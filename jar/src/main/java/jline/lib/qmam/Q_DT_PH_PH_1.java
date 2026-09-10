/**
 * @file Q_DT_PH_PH_1 - Discrete-Time DPH/DPH/1 Queue Analyzer
 *
 * Computes the queue length distribution of a discrete-time DPH/DPH/1/FCFS
 * queue. Port of Q_DT_PH_PH_1.m of the Q-MAM library by J. F. Perez,
 * J. Van Velthoven and B. Van Houdt (VALUETOOLS 2008).
 *
 * @since LINE 3.1.0
 */
package jline.lib.qmam;

import jline.util.matrix.Matrix;

public final class Q_DT_PH_PH_1 {
    private Q_DT_PH_PH_1() {}

    /**
     * Queue length distribution of a discrete-time DPH/DPH/1/FCFS queue.
     *
     * <p>A discrete phase-type law (alpha, T) with absorption vector
     * {@code t = e - T*e} renews at every event, so the arrival and service
     * streams are the DMAPs {@code (T, t*alpha)} and {@code (S, s*beta)} and
     * the queue is the DMAP/DMAP/1 one. Routing through
     * {@link Q_DT_MAP_MAP_1} keeps a single set of QBD blocks, and hence a
     * single statement of the late-arrival-with-delayed-access convention, in
     * the codebase.
     *
     * @param alpha 1 x ma initial phase vector of the interarrival DPH
     * @param T ma x ma transient matrix of the interarrival DPH
     * @param beta 1 x ms initial phase vector of the service DPH
     * @param S ms x ms transient matrix of the service DPH
     * @param options mode, maximum number of components and verbosity
     * @return the queue length distribution, entry i being Prob[i in system]
     */
    public static DTQueueResult qDtPhPh1(Matrix alpha, Matrix T, Matrix beta, Matrix S,
                                         MAPMAP1Options options) {
        Matrix C0 = T.copy();
        Matrix C1 = absorptionBlock(T, alpha);
        Matrix D0 = S.copy();
        Matrix D1 = absorptionBlock(S, beta);
        return Q_DT_MAP_MAP_1.qDtMapMap1(C0, C1, D0, D1, options);
    }

    /** Builds {@code (e - A*e) * alpha}, the event block of the renewal DMAP. */
    private static Matrix absorptionBlock(Matrix A, Matrix alpha) {
        int m = A.getNumRows();
        if (alpha.getNumCols() != m || alpha.getNumRows() != 1) {
            throw new IllegalArgumentException("The initial vector must be 1 x " + m);
        }
        Matrix rowSums = A.mult(Matrix.ones(m, 1));
        Matrix block = new Matrix(m, m);
        for (int i = 0; i < m; i++) {
            double a = 1.0 - rowSums.get(i, 0);
            for (int j = 0; j < m; j++) {
                block.set(i, j, a * alpha.get(0, j));
            }
        }
        return block;
    }
}
