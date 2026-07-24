/**
 * @file MAP/D/c queueing system analysis
 *
 * Implements analysis of MAP/D/c queues using Q-MAM solver.
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys;

import jline.api.mc.Ctmc_solve;
import jline.lib.qmam.MAPDcOptions;
import jline.lib.qmam.MAPDcResult;
import jline.lib.qmam.Q_CT_MAP_D_C;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Qsys_mapdc {
    private Qsys_mapdc() {}

    /**
     * Analyzes a MAP/D/c queue (multi-server with deterministic service).
     */
    public static QsysMapDcResult qsys_mapdc(Matrix D0, Matrix D1, double s, int c, int maxNumComp,
                                              int numSteps, int verbose) {
        MAPDcResult result = Q_CT_MAP_D_C.qCtMapDC(D0, D1, s, c, new MAPDcOptions(maxNumComp, verbose, numSteps));

        Matrix theta = Ctmc_solve.ctmc_solve(D0.add(D1));
        double lambda = theta.mult(D1).elementSum();
        double rho = lambda * s / c;

        Matrix ql = result.getQueueLength();
        double meanQL = 0.0;
        for (int i = 0; i < ql.getNumCols(); i++) {
            meanQL += i * ql.get(0, i);
        }

        Matrix w = result.getWaitingTime();

        // see _kb/03-api-layer.md for rationale
        double meanWT = Math.max(0.0, meanQL / lambda - s);

        double meanST = meanWT + s;

        return new QsysMapDcResult(meanQL, meanWT, meanST, rho, ql, w, "Q-MAM:MAP/D/" + c);
    }

    public static QsysMapDcResult qsys_mapdc(Matrix D0, Matrix D1, double s, int c, int maxNumComp, int numSteps) {
        return qsys_mapdc(D0, D1, s, c, maxNumComp, numSteps, 0);
    }

    public static QsysMapDcResult qsys_mapdc(Matrix D0, Matrix D1, double s, int c, int maxNumComp) {
        return qsys_mapdc(D0, D1, s, c, maxNumComp, 1, 0);
    }

    public static QsysMapDcResult qsys_mapdc(Matrix D0, Matrix D1, double s, int c) {
        return qsys_mapdc(D0, D1, s, c, 1000, 1, 0);
    }

    /**
     * Simplified MAP/D/c analysis using MatrixCell input for arrival.
     */
    public static QsysMapDcResult qsys_mapdc(MatrixCell arrival, double s, int c) {
        if (!(arrival.size() >= 2)) {
            throw new IllegalArgumentException("Arrival MAP must have at least 2 matrices [D0, D1]");
        }
        return qsys_mapdc(arrival.get(0), arrival.get(1), s, c, 1000, 1, 0);
    }

    /**
     * Analyzes a PH/D/c queue.
     */
    public static QsysMapDcResult qsys_phdc(Matrix alpha, Matrix T, double s, int c) {
        Matrix ones = Matrix.ones(T.getNumRows(), 1);
        Matrix exitRates = T.mult(ones).scale(-1.0);
        Matrix D1 = exitRates.mult(alpha);

        return qsys_mapdc(T, D1, s, c, 1000, 1, 0);
    }
}
