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
import jline.util.matrix.Matrix;

public final class MAPMAP1 {
    private MAPMAP1() {}

    public static Map<String, Object> mapmap1(Matrix D0, Matrix D1, Matrix S0, Matrix S1, Map<String, Object> measures) {
        return mapmap1(D0, D1, S0, S1, measures, 1e-14);
    }

    /**
     * Returns various performance measures of a continuous time MAP/MAP/1 queue.
     */
    public static Map<String, Object> mapmap1(Matrix D0, Matrix D1, Matrix S0, Matrix S1,
                                              Map<String, Object> measures, double prec) {
        int Na = D0.getNumRows();
        int Ns = S0.getNumRows();
        Matrix IA = Matrix.eye(Na);
        Matrix IS = Matrix.eye(Ns);

        Matrix B = IA.kron(S1);
        Matrix L = D0.kron(IS).add(IA.kron(S0));
        Matrix F = D1.kron(IS);
        Matrix L0 = D0.kron(IS);

        jline.util.Pair<Matrix, Matrix> sol = QBDSolve.qbdSolve(B, L, F, L0, prec);
        Matrix pi0 = sol.getFirst();
        Matrix R = sol.getSecond();
        int N = pi0.getNumCols();
        Matrix I = Matrix.eye(N);

        boolean needST = measures.containsKey("stMoms") || measures.containsKey("stDistr")
                || measures.containsKey("stDistrME") || measures.containsKey("stDistrPH");
        Matrix eta = null;
        Matrix T = null;
        Matrix Rh = null;

        if (needST) {
            Matrix U = L.add(R.mult(B));
            Rh = U.neg().inv().mult(F);
            T = IA.kron(S0).add(Rh.mult(B));
            Matrix IminusRh = I.sub(Rh);
            eta = pi0.mult(F).mult(IminusRh.inv());
            double etaSum = eta.elementSum();
            if (etaSum > 0) {
                eta = eta.scale(1.0 / etaSum);
            }
        }

        HashMap<String, Object> result = new HashMap<String, Object>();

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
            double[] outArr = new double[numMoms];
            for (int it = 0; it < numMoms; it++) outArr[it] = momsMat.get(0, it);
            result.put("ncMoms", outArr);
        }

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

        if (measures.containsKey("stMoms") && needST) {
            int numMoms = (Integer) measures.get("stMoms");
            Matrix iT = T.neg().inv();
            double[] stMoms = new double[numMoms];
            for (int m = 1; m <= numMoms; m++) {
                double fact = 1.0;
                for (int k = 1; k <= m; k++) fact *= (double) k;
                Matrix iTPow = I;
                for (int k = 0; k < m; k++) iTPow = iTPow.mult(iT);
                stMoms[m - 1] = fact * eta.mult(iTPow).elementSum();
            }
            result.put("stMoms", stMoms);
        }

        if (measures.containsKey("stDistr") && needST) {
            double[] points = (double[]) measures.get("stDistr");
            double[] values = new double[points.length];
            for (int p = 0; p < points.length; p++) {
                Matrix expTt = T.scale(points[p]).expm();
                values[p] = 1.0 - eta.mult(expTt).elementSum();
            }
            result.put("stDistr", values);
        }

        if (measures.containsKey("stDistrME") && needST) {
            result.put("stDistrME_alpha", eta);
            result.put("stDistrME_A", T);
        }

        return result;
    }
}
