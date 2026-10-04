/**
 * @file MAP/MAP/1 queueing system analysis
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys;

import java.util.HashMap;
import java.util.Map;

import jline.api.mc.Ctmc_solve;
import jline.lib.butools.MMAPPH1FCFS;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Qsys_mapmap1 {
    private Qsys_mapmap1() {}

    public static QsysMapPhResult qsys_mapmap1(Matrix C0, Matrix C1, Matrix D0, Matrix D1) {
        return qsys_mapmap1(C0, C1, D0, D1, 3, 100, 3);
    }

    /**
     * Analyzes a MAP/MAP/1 queue.
     */
    public static QsysMapPhResult qsys_mapmap1(Matrix C0, Matrix C1, Matrix D0, Matrix D1,
                                               int numQLMoms, int numQLProbs, int numSTMoms) {
        MatrixCell D = new MatrixCell(2);
        D.set(0, C0);
        D.set(1, C1);

        Pair<Matrix, Matrix> svc = mapToPh(D0, D1);
        Matrix sigma = svc.getLeft();
        Matrix S = svc.getRight();

        Map<Integer, Matrix> sigmaMap = new HashMap<Integer, Matrix>();
        sigmaMap.put(0, sigma);
        Map<Integer, Matrix> sMap = new HashMap<Integer, Matrix>();
        sMap.put(0, S);

        Map<String, Map<Integer, Matrix>> result = MMAPPH1FCFS.MMAPPH1FCFS(
                D, sigmaMap, sMap,
                numQLMoms, numQLProbs, numSTMoms,
                null, false, false, null, null);

        Matrix thetaArr = Ctmc_solve.ctmc_solve(C0.add(C1));
        double lambda = thetaArr.mult(C1).elementSum();

        Matrix thetaSvc = Ctmc_solve.ctmc_solve(D0.add(D1));
        double mu = thetaSvc.mult(D1).elementSum();
        double rho = lambda / mu;

        Map<Integer, Matrix> ncMomsCell = result.get("ncMoms");
        Matrix ncMoms = ncMomsCell == null ? null : ncMomsCell.get(0);
        double meanQL = ncMoms == null ? 0.0 : ncMoms.get(0, 0);

        Map<Integer, Matrix> stMomsCell = result.get("stNoms");
        Matrix stMoms = stMomsCell == null ? null : stMomsCell.get(0);
        double meanST = stMoms == null ? 0.0 : stMoms.get(0, 0);

        Matrix negSinv = S.scale(-1.0).inv();
        Matrix ones = Matrix.ones(S.getNumRows(), 1);
        double meanService = sigma.mult(negSinv).mult(ones).get(0, 0);

        double meanWT = Math.max(0.0, meanST - meanService);

        Map<Integer, Matrix> ncDistrCell = result.get("ncDistr");
        Matrix ncDistr = ncDistrCell == null ? null : ncDistrCell.get(0);

        return new QsysMapPhResult(meanQL, meanWT, meanST, rho, ncDistr, ncMoms, stMoms, "BUTools:MMAPPH1FCFS");
    }

    private static Pair<Matrix, Matrix> mapToPh(Matrix D0, Matrix D1) {
        Matrix theta = Ctmc_solve.ctmc_solve(D0.add(D1));

        Matrix sigma = theta.mult(D1);
        double total = sigma.elementSum();
        if (total > 0) {
            sigma.scaleEq(1.0 / total);
        } else {
            for (int i = 0; i < sigma.getNumCols(); i++) {
                sigma.set(0, i, 1.0 / sigma.getNumCols());
            }
        }

        return new Pair<Matrix, Matrix>(sigma, D0);
    }

    /**
     * Simplified MAP/MAP/1 analysis using MatrixCell inputs.
     */
    public static QsysMapPhResult qsys_mapmap1(MatrixCell arrival, MatrixCell service) {
        if (arrival.size() < 2) {
            throw new IllegalArgumentException("Arrival MAP must have at least 2 matrices [D0, D1]");
        }
        if (service.size() < 2) {
            throw new IllegalArgumentException("Service MAP must have at least 2 matrices [D0, D1]");
        }
        return qsys_mapmap1(arrival.get(0), arrival.get(1), service.get(0), service.get(1));
    }
}
