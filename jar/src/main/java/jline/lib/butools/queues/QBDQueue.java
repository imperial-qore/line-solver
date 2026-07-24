/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.queues;

import java.util.HashMap;
import java.util.Map;

import jline.lib.butools.MomsFromFactorialMoms;
import jline.lib.butools.SimilarityMatrixForVectors;
import jline.lib.butools.mam.QBDSolve;
import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class QBDQueue {
    private QBDQueue() {}

    /**
     * Returns various performance measures of a continuous time QBD queue.
     */
    public static Map<String, Object> qbdQueue(Matrix B, Matrix L, Matrix F, Matrix L0,
                                               Map<String, Object> measures, double prec) {
        Pair<Matrix, Matrix> sol = QBDSolve.qbdSolve(B, L, F, L0, prec);
        Matrix pi0 = sol.getLeft();
        Matrix R = sol.getRight();
        int N = pi0.getNumCols();
        Matrix I = Matrix.eye(N);

        boolean needST = measures.containsKey("stMoms") || measures.containsKey("stDistr")
                || measures.containsKey("stDistrME") || measures.containsKey("stDistrPH");
        Matrix eta = null;
        Matrix z = null;
        Matrix Rh = null;

        if (needST) {
            Matrix U = L.add(R.mult(B));
            Rh = U.neg().inv().mult(F);
            Matrix IminusRh = I.sub(Rh);
            eta = pi0.mult(F).mult(IminusRh.inv());
            double etaSum = eta.elementSum();
            if (etaSum > 0) {
                eta = eta.scale(1.0 / etaSum);
            }
            z = new Matrix(N * N, 1);
            for (int j = 0; j < N; j++) {
                for (int i = 0; i < N; i++) {
                    z.set(j * N + i, 0, I.get(i, j));
                }
            }
        }

        Map<String, Object> result = new HashMap<String, Object>();

        // ncMoms
        if (measures.containsKey("ncMoms")) {
            int numMoms = (Integer) measures.get("ncMoms");
            Matrix iR = I.sub(R).inv();
            Matrix factMoms = new Matrix(1, numMoms);
            for (int m = 1; m <= numMoms; m++) {
                Matrix iRpow = I;
                for (int k = 0; k <= m; k++) {
                    iRpow = iRpow.mult(iR);
                }
                Matrix Rpow = I;
                for (int k = 0; k < m; k++) {
                    Rpow = Rpow.mult(R);
                }
                double fact = 1.0;
                for (int k = 1; k <= m; k++) fact *= (double) k;
                factMoms.set(0, m - 1, fact * pi0.mult(iRpow).mult(Rpow).elementSum());
            }
            Matrix momsMat = MomsFromFactorialMoms.MomsFromFactorialMoms(factMoms);
            double[] arr = new double[numMoms];
            for (int i = 0; i < numMoms; i++) {
                arr[i] = momsMat.get(0, i);
            }
            result.put("ncMoms", arr);
        }

        // ncDistr
        if (measures.containsKey("ncDistr")) {
            int numProbs = (Integer) measures.get("ncDistr");
            double[] values = new double[numProbs];
            values[0] = pi0.elementSum();
            Matrix RPow = I;
            for (int p = 0; p < numProbs - 1; p++) {
                RPow = RPow.mult(R);
                values[p + 1] = pi0.mult(RPow).elementSum();
            }
            result.put("ncDistr", values);
        }

        // ncDistrMG
        if (measures.containsKey("ncDistrMG")) {
            Matrix iR = I.sub(R).inv();
            Matrix iRR = iR.mult(R);
            Matrix rowSums = new Matrix(N, 1);
            for (int i = 0; i < N; i++) {
                double s = 0.0;
                for (int j = 0; j < N; j++) {
                    s += iRR.get(i, j);
                }
                rowSums.set(i, 0, s);
            }
            Matrix Bsim = SimilarityMatrixForVectors.SimilarityMatrixForVectors(rowSums, Matrix.ones(N, 1));
            Matrix Bi = Bsim.inv();
            result.put("ncDistrMG_alpha", pi0.mult(Bi));
            result.put("ncDistrMG_A", Bsim.mult(R).mult(Bi));
        }

        // ncDistrDPH
        if (measures.containsKey("ncDistrDPH")) {
            Matrix iR = I.sub(R).inv();
            Matrix alpha = pi0.mult(R).mult(iR);
            double[] alphaArr = alpha.toArray1D();
            Matrix diagAlpha = Matrix.diagMatrix(alphaArr);
            Matrix diagAlphaInv = diagAlpha.inv();
            result.put("ncDistrDPH_alpha", alpha);
            result.put("ncDistrDPH_A", diagAlphaInv.mult(R.transpose()).mult(diagAlpha));
        }

        // stMoms
        if (measures.containsKey("stMoms") && needST) {
            int numMoms = (Integer) measures.get("stMoms");
            Matrix Z = L.transpose().add(F.transpose()).kron(I).add(B.transpose().kron(Rh));
            Matrix iZ = Z.neg().inv();
            Matrix keta = new Matrix(1, N * N);
            for (int j = 0; j < N; j++) {
                for (int i = 0; i < N; i++) {
                    keta.set(0, j * N + i, eta.get(0, i));
                }
            }
            double[] stMoms = new double[numMoms];
            for (int m = 1; m <= numMoms; m++) {
                double fact = 1.0;
                for (int k = 1; k <= m; k++) fact *= (double) k;
                Matrix iZpow = I.kron(I);
                Matrix iZN2 = iZ;
                for (int k = 0; k <= m; k++) {
                    iZpow = iZpow.mult(iZN2);
                }
                stMoms[m - 1] = fact * keta.mult(iZpow).mult(Z.neg()).mult(z).elementSum();
            }
            result.put("stMoms", stMoms);
        }

        // stDistr
        if (measures.containsKey("stDistr") && needST) {
            double[] points = (double[]) measures.get("stDistr");
            Matrix Z = L.transpose().add(F.transpose()).kron(I).add(B.transpose().kron(Rh));
            Matrix keta = new Matrix(1, N * N);
            for (int j = 0; j < N; j++) {
                for (int i = 0; i < N; i++) {
                    keta.set(0, j * N + i, eta.get(0, i));
                }
            }
            double[] values = new double[points.length];
            for (int p = 0; p < points.length; p++) {
                Matrix expZt = Z.scale(points[p]).expm();
                values[p] = 1.0 - keta.mult(expZt).mult(z).elementSum();
            }
            result.put("stDistr", values);
        }

        // stDistrME
        if (measures.containsKey("stDistrME") && needST) {
            Matrix Z = L.transpose().add(F.transpose()).kron(I).add(B.transpose().kron(Rh));
            Matrix keta = new Matrix(1, N * N);
            for (int j = 0; j < N; j++) {
                for (int i = 0; i < N; i++) {
                    keta.set(0, j * N + i, eta.get(0, i));
                }
            }
            Matrix Bsim = SimilarityMatrixForVectors.SimilarityMatrixForVectors(z, Matrix.ones(z.getNumRows(), 1));
            Matrix Bi = Bsim.inv();
            result.put("stDistrME_alpha", keta.mult(Bi));
            result.put("stDistrME_A", Bsim.mult(Z).mult(Bi));
        }

        return result;
    }

    public static Map<String, Object> qbdQueue(Matrix B, Matrix L, Matrix F, Matrix L0,
                                               Map<String, Object> measures) {
        return qbdQueue(B, L, F, L0, measures, 1e-14);
    }
}
