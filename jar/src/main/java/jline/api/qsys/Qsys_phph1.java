/**
 * @file PH/PH/1 queueing system analysis
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.lib.butools.MMAPPH1FCFS;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Qsys_phph1 {
    private Qsys_phph1() {}

    /**
     * Analyzes a PH/PH/1 queue.
     */
    public static QsysMapPhResult qsys_phph1(Matrix alpha, Matrix T, Matrix beta, Matrix S,
                                             int numQLMoms, int numQLProbs, int numSTMoms) {
        Pair<Matrix, Matrix> map = phToMap(alpha, T);
        Matrix D0 = map.getLeft();
        Matrix D1 = map.getRight();

        MatrixCell D = new MatrixCell(2);
        D.set(0, D0);
        D.set(1, D1);

        Map<Integer, Matrix> sigmaMap = new HashMap<Integer, Matrix>();
        sigmaMap.put(0, beta);
        Map<Integer, Matrix> sMap = new HashMap<Integer, Matrix>();
        sMap.put(0, S);

        Map<String, Map<Integer, Matrix>> result = MMAPPH1FCFS.MMAPPH1FCFS(
                D, sigmaMap, sMap,
                numQLMoms, numQLProbs, numSTMoms,
                null, false, false, null, null);

        Matrix negTinv = T.scale(-1.0).inv();
        Matrix ones = Matrix.ones(T.getNumRows(), 1);
        double meanInterarrival = alpha.mult(negTinv).mult(ones).get(0, 0);
        double lambda = 1.0 / meanInterarrival;

        Matrix negSinv = S.scale(-1.0).inv();
        Matrix onesS = Matrix.ones(S.getNumRows(), 1);
        double meanService = beta.mult(negSinv).mult(onesS).get(0, 0);
        double mu = 1.0 / meanService;
        double rho = lambda / mu;

        Matrix ncMoms = null;
        if (result.get("ncMoms") != null) {
            ncMoms = result.get("ncMoms").get(0);
        }
        double meanQL = (ncMoms != null) ? ncMoms.get(0, 0) : 0.0;

        Matrix stMoms = null;
        if (result.get("stNoms") != null) {
            stMoms = result.get("stNoms").get(0);
        }
        double meanST = (stMoms != null) ? stMoms.get(0, 0) : 0.0;

        double meanWT = Math.max(0.0, meanST - meanService);

        Matrix ncDistr = null;
        if (result.get("ncDistr") != null) {
            ncDistr = result.get("ncDistr").get(0);
        }

        return new QsysMapPhResult(
                meanQL,
                meanWT,
                meanST,
                rho,
                ncDistr,
                ncMoms,
                stMoms,
                "BUTools:MMAPPH1FCFS");
    }

    public static QsysMapPhResult qsys_phph1(Matrix alpha, Matrix T, Matrix beta, Matrix S,
                                             int numQLMoms, int numQLProbs) {
        return qsys_phph1(alpha, T, beta, S, numQLMoms, numQLProbs, 3);
    }

    public static QsysMapPhResult qsys_phph1(Matrix alpha, Matrix T, Matrix beta, Matrix S, int numQLMoms) {
        return qsys_phph1(alpha, T, beta, S, numQLMoms, 100, 3);
    }

    public static QsysMapPhResult qsys_phph1(Matrix alpha, Matrix T, Matrix beta, Matrix S) {
        return qsys_phph1(alpha, T, beta, S, 3, 100, 3);
    }

    /**
     * Converts a PH distribution to its equivalent MAP representation.
     */
    private static Pair<Matrix, Matrix> phToMap(Matrix alpha, Matrix T) {
        Matrix D0 = T.copy();

        Matrix ones = Matrix.ones(T.getNumRows(), 1);
        Matrix exitRates = T.mult(ones).scale(-1.0);

        Matrix D1 = exitRates.mult(alpha);

        return new Pair<Matrix, Matrix>(D0, D1);
    }

    /**
     * Simplified PH/PH/1 analysis using MatrixCell inputs.
     */
    public static QsysMapPhResult qsys_phph1(MatrixCell arrival, MatrixCell service) {
        Pair<Matrix, Matrix> aPair = extractPH(arrival, "arrival");
        Pair<Matrix, Matrix> sPair = extractPH(service, "service");
        return qsys_phph1(aPair.getLeft(), aPair.getRight(), sPair.getLeft(), sPair.getRight());
    }

    /**
     * Extracts PH parameters from MatrixCell.
     */
    private static Pair<Matrix, Matrix> extractPH(MatrixCell ph, String name) {
        if (ph.size() < 1) {
            throw new IllegalArgumentException(name + " PH must have at least 1 matrix");
        }

        if (ph.size() >= 2) {
            Matrix alpha = (ph.get(0).getNumRows() == 1) ? ph.get(0) : ph.get(0).transpose();
            Matrix T = ph.get(1);
            return new Pair<Matrix, Matrix>(alpha, T);
        } else {
            Matrix T = ph.get(0);
            int n = T.getNumRows();
            Matrix alpha = Matrix.ones(1, n).scale(1.0 / n);
            return new Pair<Matrix, Matrix>(alpha, T);
        }
    }

    /**
     * Analyzes a PH/PH/1 queue with exponential service (simplified E/M/1).
     */
    public static QsysMapPhResult qsys_phm1(Matrix alpha, Matrix T, double mu) {
        Matrix beta = new Matrix(1, 1);
        beta.set(0, 0, 1.0);
        Matrix S = new Matrix(1, 1);
        S.set(0, 0, -mu);

        return qsys_phph1(alpha, T, beta, S);
    }
}
