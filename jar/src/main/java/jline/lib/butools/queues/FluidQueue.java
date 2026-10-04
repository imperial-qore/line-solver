/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.queues;

import java.util.HashMap;
import java.util.Map;

import jline.lib.butools.SimilarityMatrixForVectors;
import jline.lib.butools.mam.GeneralFluidSolve;
import jline.lib.butools.mam.GeneralFluidSolution;
import jline.util.matrix.Matrix;

public final class FluidQueue {
    private FluidQueue() {}

    /**
     * Returns various performance measures of a fluid queue.
     */
    public static Map<String, Object> fluidQueue(Matrix Q, Matrix Rin, Matrix Rout,
                                                 Map<String, Object> measures, Matrix Q0, double prec) {
        Matrix Rdiff = Rin.sub(Rout);
        GeneralFluidSolution sol = GeneralFluidSolve.generalFluidSolve(Q, Rdiff, Q0, prec);
        Matrix mass0 = sol.getMass0();
        Matrix ini = sol.getIni();
        Matrix K = sol.getK();
        Matrix clo = sol.getClo();

        boolean needST = measures.containsKey("stMoms") || measures.containsKey("stDistr")
                || measures.containsKey("stDistrME") || measures.containsKey("stDistrPH");

        Matrix iniKi = null;
        double lambd = 0.0;
        int Nq = Q.getNumRows();

        if (needST) {
            iniKi = K.transpose().leftMatrixDivide(ini.transpose().neg()).transpose();
            lambd = mass0.mult(Rin).elementSum() + iniKi.mult(clo).mult(Rin).elementSum();
        }

        Map<String, Object> result = new HashMap<String, Object>();

        // flMoms
        if (measures.containsKey("flMoms")) {
            int numMoms = (Integer) measures.get("flMoms");
            Matrix iK = K.neg().inv();
            double[] flMoms = new double[numMoms];
            for (int m = 1; m <= numMoms; m++) {
                double fact = 1.0;
                for (int k = 1; k <= m; k++) fact *= (double) k;
                Matrix iKpow = Matrix.eye(K.getNumRows());
                for (int k = 0; k <= m; k++) {
                    iKpow = iKpow.mult(iK);
                }
                flMoms[m - 1] = fact * ini.mult(iKpow).mult(clo).elementSum();
            }
            result.put("flMoms", flMoms);
        }

        // flDistr
        if (measures.containsKey("flDistr")) {
            double[] points = (double[]) measures.get("flDistr");
            Matrix iK = K.neg().inv();
            Matrix Ifl = Matrix.eye(K.getNumRows());
            double mass0Sum = mass0.elementSum();
            double[] values = new double[points.length];
            for (int p = 0; p < points.length; p++) {
                Matrix expKt = K.scale(points[p]).expm();
                values[p] = mass0Sum + ini.mult(Ifl.sub(expKt)).mult(iK).mult(clo).elementSum();
            }
            result.put("flDistr", values);
        }

        // flDistrPH
        if (measures.containsKey("flDistrPH")) {
            Matrix xCol = K.transpose().leftMatrixDivide(ini.transpose().neg());
            double[] deltaArr = xCol.toArray1D();
            Matrix Delta = Matrix.diagMatrix(deltaArr);
            Matrix DeltaInv = Delta.inv();
            result.put("flDistrPH_A", DeltaInv.mult(K.transpose()).mult(Delta));
            Matrix cloRowSumsPH = clo.sumRows();
            result.put("flDistrPH_alpha", cloRowSumsPH.transpose().mult(Delta));
        }

        // flDistrME
        if (measures.containsKey("flDistrME")) {
            Matrix iK = K.neg().inv();
            Matrix cloRowSums = new Matrix(clo.getNumRows(), 1);
            for (int i = 0; i < clo.getNumRows(); i++) {
                double s = 0.0;
                for (int j = 0; j < clo.getNumCols(); j++) {
                    s += clo.get(i, j);
                }
                cloRowSums.set(i, 0, s);
            }
            Matrix target = iK.mult(cloRowSums);
            Matrix Bsim = SimilarityMatrixForVectors.SimilarityMatrixForVectors(target,
                    Matrix.ones(K.getNumRows(), 1));
            Matrix Bi = Bsim.inv();
            result.put("flDistrME_alpha", ini.mult(Bi));
            result.put("flDistrME_A", Bsim.mult(K).mult(Bi));
        }

        // stMoms
        if (measures.containsKey("stMoms") && needST) {
            int numMoms = (Integer) measures.get("stMoms");
            Matrix Ifl = Matrix.eye(K.getNumRows());
            Matrix Z = Q.transpose().kron(Ifl).add(Rout.kron(K));
            Matrix iZ = Z.neg().inv();
            Matrix iniScaled = ini.scale(1.0 / lambd);
            Matrix kini = new Matrix(1, Nq * iniScaled.getNumCols());
            for (int j = 0; j < Nq; j++) {
                for (int i = 0; i < iniScaled.getNumCols(); i++) {
                    kini.set(0, j * iniScaled.getNumCols() + i, iniScaled.get(0, i));
                }
            }
            Matrix iK = K.neg().inv();
            Matrix kcloMat = iK.mult(clo).mult(Rin);
            Matrix kclo = new Matrix(Nq * ini.getNumCols(), 1);
            for (int j = 0; j < Nq; j++) {
                for (int i = 0; i < ini.getNumCols(); i++) {
                    kclo.set(j * ini.getNumCols() + i, 0, kcloMat.get(i, j));
                }
            }
            double[] stMoms = new double[numMoms];
            for (int m = 1; m <= numMoms; m++) {
                double fact = 1.0;
                for (int k = 1; k <= m; k++) fact *= (double) k;
                Matrix iZpow = Matrix.eye(Z.getNumRows());
                for (int k = 0; k <= m; k++) {
                    iZpow = iZpow.mult(iZ);
                }
                stMoms[m - 1] = fact * kini.mult(iZpow).mult(Z.neg()).mult(kclo).elementSum();
            }
            result.put("stMoms", stMoms);
        }

        // stDistr
        if (measures.containsKey("stDistr") && needST) {
            double[] points = (double[]) measures.get("stDistr");
            Matrix Ifl = Matrix.eye(K.getNumRows());
            Matrix Z = Q.transpose().kron(Ifl).add(Rout.kron(K));
            Matrix iniScaled = ini.scale(1.0 / lambd);
            Matrix kini = new Matrix(1, Nq * iniScaled.getNumCols());
            for (int j = 0; j < Nq; j++) {
                for (int i = 0; i < iniScaled.getNumCols(); i++) {
                    kini.set(0, j * iniScaled.getNumCols() + i, iniScaled.get(0, i));
                }
            }
            Matrix iK = K.neg().inv();
            Matrix kcloMat = iK.mult(clo).mult(Rin);
            Matrix kclo = new Matrix(Nq * ini.getNumCols(), 1);
            for (int j = 0; j < Nq; j++) {
                for (int i = 0; i < ini.getNumCols(); i++) {
                    kclo.set(j * ini.getNumCols() + i, 0, kcloMat.get(i, j));
                }
            }
            Matrix iZ = Z.neg().inv();
            double[] values = new double[points.length];
            for (int p = 0; p < points.length; p++) {
                Matrix expZt = Z.scale(points[p]).expm();
                values[p] = 1.0 - kini.mult(expZt).mult(iZ).mult(Z.neg()).mult(kclo).elementSum();
            }
            result.put("stDistr", values);
        }

        // stDistrME
        if (measures.containsKey("stDistrME") && needST) {
            Matrix Ifl = Matrix.eye(K.getNumRows());
            Matrix Z = Q.transpose().kron(Ifl).add(Rout.kron(K));
            Matrix iniScaled = ini.scale(1.0 / lambd);
            Matrix kini = new Matrix(1, Nq * iniScaled.getNumCols());
            for (int j = 0; j < Nq; j++) {
                for (int i = 0; i < iniScaled.getNumCols(); i++) {
                    kini.set(0, j * iniScaled.getNumCols() + i, iniScaled.get(0, i));
                }
            }
            Matrix iK = K.neg().inv();
            Matrix kcloMat = iK.mult(clo).mult(Rin);
            Matrix kclo = new Matrix(Nq * ini.getNumCols(), 1);
            for (int j = 0; j < Nq; j++) {
                for (int i = 0; i < ini.getNumCols(); i++) {
                    kclo.set(j * ini.getNumCols() + i, 0, kcloMat.get(i, j));
                }
            }
            Matrix Bsim = SimilarityMatrixForVectors.SimilarityMatrixForVectors(kclo,
                    Matrix.ones(kclo.getNumRows(), 1));
            Matrix Bi = Bsim.inv();
            result.put("stDistrME_alpha", kini.mult(Bi));
            result.put("stDistrME_A", Bsim.mult(Z).mult(Bi));
        }

        return result;
    }

    public static Map<String, Object> fluidQueue(Matrix Q, Matrix Rin, Matrix Rout,
                                                 Map<String, Object> measures, Matrix Q0) {
        return fluidQueue(Q, Rin, Rout, measures, Q0, 1e-14);
    }

    public static Map<String, Object> fluidQueue(Matrix Q, Matrix Rin, Matrix Rout,
                                                 Map<String, Object> measures) {
        return fluidQueue(Q, Rin, Rout, measures, null, 1e-14);
    }
}
