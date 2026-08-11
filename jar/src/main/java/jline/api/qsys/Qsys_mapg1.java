/**
 * @file MAP/G/1 queueing system analysis
 *
 * Implements analysis of MAP/G/1 queues using BUTools MMAPPH1FCFS solver.
 * The general service time is fitted to a Phase-Type distribution using moment matching.
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys;

import java.util.HashMap;
import java.util.Map;

import jline.api.mc.Ctmc_solve;
import jline.lang.processes.APH;
import jline.lib.butools.APHFrom3Moments;
import jline.lib.butools.MMAPPH1FCFS;
import jline.lib.butools.ph.PH2From3Moments;
import jline.lib.butools.ph.PH2From3Moments.PH2Representation;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Qsys_mapg1 {
    private Qsys_mapg1() {}

    public static QsysMapPhResult qsys_mapg1(Matrix D0, Matrix D1, double[] serviceMoments) {
        return qsys_mapg1(D0, D1, serviceMoments, 3, 100, 3);
    }

    public static QsysMapPhResult qsys_mapg1(Matrix D0, Matrix D1, double[] serviceMoments, int numQLMoms) {
        return qsys_mapg1(D0, D1, serviceMoments, numQLMoms, 100, 3);
    }

    public static QsysMapPhResult qsys_mapg1(Matrix D0, Matrix D1, double[] serviceMoments,
                                             int numQLMoms, int numQLProbs) {
        return qsys_mapg1(D0, D1, serviceMoments, numQLMoms, numQLProbs, 3);
    }

    /**
     * Analyzes a MAP/G/1 queue.
     */
    public static QsysMapPhResult qsys_mapg1(Matrix D0, Matrix D1, double[] serviceMoments,
                                             int numQLMoms, int numQLProbs, int numSTMoms) {
        if (serviceMoments.length < 2) {
            throw new IllegalArgumentException("At least 2 service moments are required");
        }

        // Fit service distribution to PH using moment matching
        Matrix[] phPair = fitServiceToPH(serviceMoments);
        Matrix sigma = phPair[0];
        Matrix S = phPair[1];

        // Build arrival MMAP structure for BUTools (single class)
        MatrixCell D = new MatrixCell(2);
        D.set(0, D0);
        D.set(1, D1);

        // Service parameters as maps (single class indexed by 0)
        Map<Integer, Matrix> sigmaMap = new HashMap<Integer, Matrix>();
        sigmaMap.put(0, sigma);
        Map<Integer, Matrix> sMap = new HashMap<Integer, Matrix>();
        sMap.put(0, S);

        Map<String, Map<Integer, Matrix>> result = MMAPPH1FCFS.MMAPPH1FCFS(
                D, sigmaMap, sMap,
                numQLMoms, numQLProbs, numSTMoms,
                null, false, false, null, null);

        // Compute utilization from arrival and service rates
        Matrix theta = Ctmc_solve.ctmc_solve(D0.add(D1));
        double lambda = theta.mult(D1).elementSum();

        double meanService = serviceMoments[0];
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

        return new QsysMapPhResult(meanQL, meanWT, meanST, rho, ncDistr, ncMoms, stMoms,
                "BUTools:MMAPPH1FCFS");
    }

    /**
     * Fits a general service time distribution to a Phase-Type representation.
     *
     * @return Two-element array {alpha, A}
     */
    private static Matrix[] fitServiceToPH(double[] moments) {
        double m1 = moments[0];

        if (moments.length >= 3) {
            try {
                APH aph = APHFrom3Moments.APHFrom3Moments(moments);
                Matrix alpha = aph.getInitProb();
                Matrix T = (Matrix) aph.getParam(3).getValue();
                return new Matrix[] { alpha, T };
            } catch (Exception e) {
                try {
                    PH2Representation ph2 = PH2From3Moments.ph2From3Moments(moments);
                    return new Matrix[] { ph2.alpha, ph2.A };
                } catch (Exception e2) {
                    return createExponentialPH(m1);
                }
            }
        } else if (moments.length >= 2) {
            double m2 = moments[1];
            double cv2 = m2 / (m1 * m1) - 1.0;

            if (cv2 <= 0.0) {
                return createErlangPH(m1, Math.max(1, (int) (1.0 / Math.max(cv2, 0.01))));
            } else if (cv2 < 1.0) {
                return createErlang2PH(m1, cv2);
            } else if (cv2 == 1.0) {
                return createExponentialPH(m1);
            } else {
                return createHyperexp2PH(m1, cv2);
            }
        } else {
            return createExponentialPH(m1);
        }
    }

    private static Matrix[] createExponentialPH(double mean) {
        Matrix alpha = new Matrix(1, 1);
        alpha.set(0, 0, 1.0);
        Matrix A = new Matrix(1, 1);
        A.set(0, 0, -1.0 / mean);
        return new Matrix[] { alpha, A };
    }

    private static Matrix[] createErlangPH(double mean, int k) {
        double mu = (double) k / mean;
        Matrix alpha = new Matrix(1, k);
        alpha.set(0, 0, 1.0);
        Matrix A = new Matrix(k, k);
        for (int i = 0; i < k; i++) {
            A.set(i, i, -mu);
            if (i < k - 1) {
                A.set(i, i + 1, mu);
            }
        }
        return new Matrix[] { alpha, A };
    }

    private static Matrix[] createErlang2PH(double mean, double cv2) {
        int k = Math.max(2, (int) (1.0 / cv2));
        return createErlangPH(mean, k);
    }

    private static Matrix[] createHyperexp2PH(double mean, double cv2) {
        double p = 0.5 * (1.0 + Math.sqrt((cv2 - 1.0) / (cv2 + 1.0)));
        double lambda1 = 2.0 * p / mean;
        double lambda2 = 2.0 * (1.0 - p) / mean;

        Matrix alpha = new Matrix(1, 2);
        alpha.set(0, 0, p);
        alpha.set(0, 1, 1.0 - p);

        Matrix A = new Matrix(2, 2);
        A.set(0, 0, -lambda1);
        A.set(1, 1, -lambda2);

        return new Matrix[] { alpha, A };
    }

    /**
     * Analyzes a MAP/G/1 queue with service time specified as mean and CV.
     */
    public static QsysMapPhResult qsys_mapg1(Matrix D0, Matrix D1, double meanService, double cvService) {
        double m1 = meanService;
        double m2 = m1 * m1 * (1.0 + cvService * cvService);
        return qsys_mapg1(D0, D1, new double[] { m1, m2 });
    }

    /**
     * Analyzes a MAP/G/1 queue using MatrixCell input for arrival.
     */
    public static QsysMapPhResult qsys_mapg1(MatrixCell arrival, double[] serviceMoments) {
        if (arrival.size() < 2) {
            throw new IllegalArgumentException("Arrival MAP must have at least 2 matrices [D0, D1]");
        }
        return qsys_mapg1(arrival.get(0), arrival.get(1), serviceMoments);
    }
}
