/**
 * @file MAP/M/1 queueing system analysis
 *
 * Implements analysis of MAP/M/1 queues as a convenience wrapper
 * for qsys_mapmc with a single server (c=1).
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Qsys_mapm1 {
    private Qsys_mapm1() {}

    /**
     * Analyzes a MAP/M/1 queue (single server with exponential service).
     *
     * @param D0 MAP hidden transition matrix (n x n)
     * @param D1 MAP arrival transition matrix (n x n)
     * @param mu Exponential service rate
     * @param maxNumComp Maximum number of queue length probabilities
     * @return QsysMapPhResult with performance metrics and analyzer set to "Q-MAM:MAP/M/1"
     */
    public static QsysMapPhResult qsys_mapm1(Matrix D0, Matrix D1, double mu, int maxNumComp) {
        QsysMapPhResult result = Qsys_mapmc.qsys_mapmc(D0, D1, mu, 1, "SylvesCR", maxNumComp);
        return new QsysMapPhResult(
                result.getMeanQueueLength(),
                result.getMeanWaitingTime(),
                result.getMeanSojournTime(),
                result.getUtilization(),
                result.getQueueLengthDist(),
                result.getQueueLengthMoments(),
                result.getSojournTimeMoments(),
                "Q-MAM:MAP/M/1"
        );
    }

    public static QsysMapPhResult qsys_mapm1(Matrix D0, Matrix D1, double mu) {
        return qsys_mapm1(D0, D1, mu, 1000);
    }

    /**
     * Simplified MAP/M/1 analysis using MatrixCell input for arrival.
     *
     * @param arrival MAP arrival process as MatrixCell [D0, D1]
     * @param mu Exponential service rate
     * @return QsysMapPhResult with performance metrics
     */
    public static QsysMapPhResult qsys_mapm1(MatrixCell arrival, double mu) {
        if (arrival.size() < 2) {
            throw new IllegalArgumentException("Arrival MAP must have at least 2 matrices [D0, D1]");
        }
        return qsys_mapm1(arrival.get(0), arrival.get(1), mu);
    }
}
