/**
 * @file MAP/PH/1 queueing system analysis
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys;

import java.util.HashMap;
import java.util.Map;

import jline.api.mc.Ctmc_solve;
import jline.lib.butools.MMAPPH1FCFS;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Qsys_mapph1 {
    private Qsys_mapph1() {}

    public static QsysMapPhResult qsys_mapph1(Matrix D0, Matrix D1, Matrix sigma, Matrix S) {
        return qsys_mapph1(D0, D1, sigma, S, 3, 100, 3);
    }

    public static QsysMapPhResult qsys_mapph1(Matrix D0, Matrix D1, Matrix sigma, Matrix S, int numQLMoms) {
        return qsys_mapph1(D0, D1, sigma, S, numQLMoms, 100, 3);
    }

    public static QsysMapPhResult qsys_mapph1(Matrix D0, Matrix D1, Matrix sigma, Matrix S,
                                              int numQLMoms, int numQLProbs) {
        return qsys_mapph1(D0, D1, sigma, S, numQLMoms, numQLProbs, 3);
    }

    /**
     * Analyzes a MAP/PH/1 queue.
     */
    public static QsysMapPhResult qsys_mapph1(Matrix D0, Matrix D1, Matrix sigma, Matrix S,
                                              int numQLMoms, int numQLProbs, int numSTMoms) {
        Matrix sigmaRow;
        if (sigma.getNumRows() > 1 && sigma.getNumCols() == 1) {
            sigmaRow = sigma.transpose();
        } else {
            sigmaRow = sigma;
        }

        MatrixCell D = new MatrixCell(2);
        D.set(0, D0);
        D.set(1, D1);

        Map<Integer, Matrix> sigmaMap = new HashMap<Integer, Matrix>();
        sigmaMap.put(0, sigmaRow);
        Map<Integer, Matrix> sMap = new HashMap<Integer, Matrix>();
        sMap.put(0, S);

        Map<String, Map<Integer, Matrix>> result = MMAPPH1FCFS.MMAPPH1FCFS(
                D, sigmaMap, sMap,
                numQLMoms, numQLProbs, numSTMoms,
                null, false, false, null, null);

        Matrix theta = Ctmc_solve.ctmc_solve(D0.add(D1));
        double lambda = theta.mult(D1).elementSum();

        Matrix negSinv = S.scale(-1.0).inv();
        Matrix ones = Matrix.ones(S.getNumRows(), 1);
        double meanService = sigmaRow.mult(negSinv).mult(ones).get(0, 0);
        double mu = 1.0 / meanService;
        double rho = lambda / mu;

        Map<Integer, Matrix> ncMomsCell = result.get("ncMoms");
        Matrix ncMoms = ncMomsCell == null ? null : ncMomsCell.get(0);
        double meanQL = ncMoms == null ? 0.0 : ncMoms.get(0, 0);

        Map<Integer, Matrix> stMomsCell = result.get("stNoms");
        Matrix stMoms = stMomsCell == null ? null : stMomsCell.get(0);
        double meanST = stMoms == null ? 0.0 : stMoms.get(0, 0);

        double meanWT = Math.max(0.0, meanST - meanService);

        Map<Integer, Matrix> ncDistrCell = result.get("ncDistr");
        Matrix ncDistr = ncDistrCell == null ? null : ncDistrCell.get(0);

        return new QsysMapPhResult(meanQL, meanWT, meanST, rho, ncDistr, ncMoms, stMoms, "BUTools:MMAPPH1FCFS");
    }

    /**
     * Simplified MAP/PH/1 analysis using MatrixCell inputs.
     */
    public static QsysMapPhResult qsys_mapph1(MatrixCell arrival, MatrixCell service) {
        Matrix D0 = arrival.get(0);
        Matrix D1 = arrival.get(1);

        Matrix sigma;
        Matrix S;
        if (service.size() >= 2) {
            sigma = service.get(0);
            S = service.get(1);
        } else {
            S = service.get(0);
            int n = S.getNumRows();
            sigma = Matrix.ones(1, n).scale(1.0 / n);
        }

        return qsys_mapph1(D0, D1, sigma, S);
    }
}
