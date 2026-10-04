/**
 * @file MAP/M/c queueing system analysis
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys;

import jline.api.mc.Ctmc_solve;
import jline.lib.qmam.MAPMcOptions;
import jline.lib.qmam.Q_CT_MAP_M_C;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Qsys_mapmc {
    private Qsys_mapmc() {}

    /**
     * Analyzes a MAP/M/c queue (multi-server with exponential service).
     */
    public static QsysMapPhResult qsys_mapmc(Matrix D0, Matrix D1, double mu, int c,
                                             String mode, int maxNumComp, int verbose) {
        jline.lib.qmam.MAPMcResult result = Q_CT_MAP_M_C.qCtMapMC(D0, D1, mu, c,
                new MAPMcOptions(mode, maxNumComp, verbose));

        Matrix theta = Ctmc_solve.ctmc_solve(D0.add(D1));
        double lambda = theta.mult(D1).elementSum();
        double rho = lambda / (mu * c);

        Matrix ql = result.getQueueLength();
        double meanQL = 0.0;
        for (int i = 0; i < ql.getNumCols(); i++) {
            meanQL += i * ql.get(0, i);
        }

        double meanWT = 0.0;
        if (result.getWaitAlpha() != null && result.getSmat() != null) {
            Matrix negSinv = result.getSmat().scale(-1.0).inv();
            Matrix ones = Matrix.ones(result.getSmat().getNumRows(), 1);
            meanWT = result.getWaitAlpha().mult(negSinv).mult(ones).get(0, 0);
        }

        double meanService = 1.0 / mu;
        double meanST = meanWT + meanService;

        return new QsysMapPhResult(meanQL, meanWT, meanST, rho, ql, null, null, "Q-MAM:MAP/M/" + c);
    }

    public static QsysMapPhResult qsys_mapmc(Matrix D0, Matrix D1, double mu, int c,
                                             String mode, int maxNumComp) {
        return qsys_mapmc(D0, D1, mu, c, mode, maxNumComp, 0);
    }

    public static QsysMapPhResult qsys_mapmc(Matrix D0, Matrix D1, double mu, int c) {
        return qsys_mapmc(D0, D1, mu, c, "SylvesCR", 1000, 0);
    }

    public static QsysMapPhResult qsys_mapmc(MatrixCell arrival, double mu, int c) {
        if (arrival.size() < 2) {
            throw new IllegalArgumentException("Arrival MAP must have at least 2 matrices [D0, D1]");
        }
        return qsys_mapmc(arrival.get(0), arrival.get(1), mu, c);
    }

    /**
     * Analyzes a PH/M/c queue.
     */
    public static QsysMapPhResult qsys_phmc(Matrix alpha, Matrix T, double mu, int c) {
        Matrix ones = Matrix.ones(T.getNumRows(), 1);
        Matrix exitRates = T.mult(ones).scale(-1.0);
        Matrix D1 = exitRates.mult(alpha);
        return qsys_mapmc(T, D1, mu, c);
    }
}
