/**
 * @file MAP/D/1 queueing system analysis
 *
 * Implements analysis of MAP/D/1 queues as a convenience wrapper for qsys_mapdc with c=1.
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Qsys_mapd1 {
    private Qsys_mapd1() {}

    /**
     * Analyzes a MAP/D/1 queue (single server with deterministic service).
     *
     * @param D0 MAP hidden transition matrix (n x n)
     * @param D1 MAP arrival transition matrix (n x n)
     * @param s Deterministic service time (positive scalar)
     * @param maxNumComp Maximum number of queue length components
     * @param numSteps Number of waiting time distribution points per service interval
     * @return QsysMapDcResult with performance metrics
     */
    public static QsysMapDcResult qsys_mapd1(
            Matrix D0,
            Matrix D1,
            double s,
            int maxNumComp,
            int numSteps) {
        return Qsys_mapdc.qsys_mapdc(D0, D1, s, 1, maxNumComp, numSteps);
    }

    public static QsysMapDcResult qsys_mapd1(Matrix D0, Matrix D1, double s, int maxNumComp) {
        return Qsys_mapdc.qsys_mapdc(D0, D1, s, 1, maxNumComp, 1);
    }

    public static QsysMapDcResult qsys_mapd1(Matrix D0, Matrix D1, double s) {
        return Qsys_mapdc.qsys_mapdc(D0, D1, s, 1, 1000, 1);
    }

    /**
     * Simplified MAP/D/1 analysis using MatrixCell input for arrival.
     *
     * @param arrival MAP arrival process as MatrixCell [D0, D1]
     * @param s Deterministic service time
     * @return QsysMapDcResult with performance metrics
     */
    public static QsysMapDcResult qsys_mapd1(MatrixCell arrival, double s) {
        if (!(arrival.size() >= 2)) {
            throw new IllegalArgumentException("Arrival MAP must have at least 2 matrices [D0, D1]");
        }
        return Qsys_mapdc.qsys_mapdc(arrival.get(0), arrival.get(1), s, 1);
    }
}
