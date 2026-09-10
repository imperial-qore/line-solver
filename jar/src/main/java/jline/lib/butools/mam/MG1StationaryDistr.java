/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.mam;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import jline.lib.butools.mc.DTMCSolve;
import jline.util.matrix.Matrix;

public final class MG1StationaryDistr {
    private MG1StationaryDistr() {}

    /**
     * Returns the stationary distribution of the M/G/1 type Markov chain
     * up to a given level K.
     *
     * Uses the stable Ramaswami formula for computing the stationary distribution.
     *
     * @param A List of matrix blocks of the M/G/1 type generator in the regular part.
     * @param B List of matrix blocks at the boundary, or null if B=A is assumed.
     * @param G Matrix G of the M/G/1 type Markov chain, or null to compute from A.
     * @param K The stationary distribution is returned up to this level.
     * @param prec Numerical precision.
     * @return Stationary probability vector up to level K, shape (1, (K+1)*N)
     */
    public static Matrix mg1StationaryDistr(List<Matrix> A, List<Matrix> B, Matrix G, int K, double prec) {
        Matrix Gmat = (G != null) ? G : MG1FundamentalMatrix.mg1FundamentalMatrix(A, prec);

        Matrix g = DTMCSolve.dtmcSolve(Gmat);

        int m = A.get(0).getNumRows();
        Matrix I = Matrix.eye(m);
        int dega = A.size() - 1;

        Matrix[] Awork = new Matrix[A.size()];
        for (int i = 0; i < A.size(); i++) {
            Awork[i] = A.get(i).copy();
        }

        Matrix sumA = Awork[dega].copy();
        Matrix beta = sumA.sumRows();

        for (int i = dega - 1; i >= 1; i--) {
            sumA = sumA.add(1.0, Awork[i]);
            Awork[i] = Awork[i].add(1.0, Awork[i + 1].mult(Gmat));
            beta = beta.add(1.0, sumA.sumRows());
        }
        sumA = sumA.add(1.0, Awork[0]);

        Matrix theta = DTMCSolve.dtmcSolve(sumA);
        double drift = theta.mult(beta).get(0, 0);

        if (drift >= 1.0) {
            throw new IllegalStateException("MG1StationaryDistr: The Markov chain is not positive recurrent (drift=" + drift + ")");
        }

        boolean useBoundary = B != null;

        if (useBoundary) {
            int mb = B.size();
            int degb = mb - 1;
            Matrix[] Bwork = new Matrix[mb];
            for (int i = 0; i < mb; i++) {
                Bwork[i] = B.get(i).copy();
            }

            Matrix sumBB0 = Bwork[degb].copy();
            Matrix Bbeta = Matrix.zeros(m, 1);
            for (int i = degb - 1; i >= 1; i--) {
                Matrix rowSums = Matrix.zeros(m, 1);
                for (int r = 0; r < m; r++) {
                    double sum = 0.0;
                    for (int c = 0; c < m; c++) {
                        sum += sumBB0.get(r, c);
                    }
                    rowSums.set(r, 0, sum);
                }
                for (int r = 0; r < m; r++) {
                    Bbeta.set(r, 0, Bbeta.get(r, 0) + rowSums.get(r, 0));
                }
                sumBB0 = sumBB0.add(1.0, Bwork[i]);
                Bwork[i] = Bwork[i].add(1.0, Bwork[i + 1].mult(Gmat));
            }

            Matrix Km = Bwork[0].add(1.0, Bwork[1].mult(Gmat));
            Matrix kappa = DTMCSolve.dtmcSolve(Km);

            Matrix onesM = Matrix.ones(m, 1);
            Matrix betaOnes = beta.copy();
            Matrix gColTimesOnesRow = new Matrix(m, m);
            for (int r = 0; r < m; r++) {
                for (int c = 0; c < m; c++) {
                    gColTimesOnesRow.set(r, c, onesM.get(r, 0) * g.get(0, c));
                }
            }
            Matrix betaMinusOnes = new Matrix(m, 1);
            for (int r = 0; r < m; r++) {
                betaMinusOnes.set(r, 0, onesM.get(r, 0) - betaOnes.get(r, 0));
            }
            Matrix outerProd = new Matrix(m, m);
            for (int r = 0; r < m; r++) {
                for (int c = 0; c < m; c++) {
                    outerProd.set(r, c, betaMinusOnes.get(r, 0) * g.get(0, c));
                }
            }
            Matrix tempMat = I.sub(sumA).sub(1.0, outerProd);
            Matrix tempInv = tempMat.inv();
            Matrix tempSum = tempInv.sumRows();

            Matrix IminusA0A1 = I.sub(Awork[0]).sub(1.0, Awork[1]);
            Matrix psi1 = IminusA0A1.mult(tempSum);
            Matrix A0sum = Awork[0].sumRows();
            double scaleFactor = 1.0 / (1.0 - drift);
            for (int r = 0; r < m; r++) {
                psi1.set(r, 0, psi1.get(r, 0) + scaleFactor * A0sum.get(r, 0));
            }

            Matrix sumBB0minusB1 = sumBB0.sub(Bwork[1]);
            Matrix psi2 = sumBB0minusB1.mult(tempSum);
            for (int r = 0; r < m; r++) {
                psi2.set(r, 0, psi2.get(r, 0) + 1.0 + scaleFactor * Bbeta.get(r, 0));
            }

            Matrix IminusA1inv = I.sub(Awork[1]).inv();
            Matrix tildekappa1 = psi2.add(1.0, Bwork[1].mult(IminusA1inv).mult(psi1));

            double normVal = kappa.mult(tildekappa1).get(0, 0);
            Matrix pi0 = kappa.scale(1.0 / normVal);

            Matrix invbarA1 = I.sub(Awork[1]).inv();
            List<Matrix> piList = new ArrayList<Matrix>();
            piList.add(pi0);
            double sumpi = pi0.elementSum();
            int numit = 1;

            while (sumpi < 1.0 - 1e-10 && numit <= K) {
                Matrix pix;
                if (numit <= degb) {
                    pix = pi0.mult(Bwork[numit]);
                } else {
                    pix = Matrix.zeros(1, m);
                }
                int jMax = (numit < dega) ? numit : dega;
                for (int j = 1; j < jMax; j++) {
                    pix = pix.add(1.0, piList.get(numit - j).mult(Awork[j + 1]));
                }
                pix = pix.mult(invbarA1);
                sumpi += pix.elementSum();
                piList.add(pix);
                numit++;
            }

            Matrix result = Matrix.zeros(1, (K + 1) * m);
            int kCap = Math.min(piList.size(), K + 1);
            for (int k = 0; k < kCap; k++) {
                for (int j = 0; j < m; j++) {
                    result.set(0, k * m + j, piList.get(k).get(0, j));
                }
            }
            return result;
        } else {
            Matrix pi0 = g.scale(1.0 - drift);

            Matrix invbarA1 = I.sub(Awork[1]).inv();
            List<Matrix> piList = new ArrayList<Matrix>();
            piList.add(pi0);
            double sumpi = pi0.elementSum();
            int numit = 1;

            while (sumpi < 1.0 - 1e-10 && numit <= K) {
                Matrix pix;
                if (numit <= dega) {
                    pix = pi0.mult(Awork[numit]);
                } else {
                    pix = Matrix.zeros(1, m);
                }
                int jMax = (numit < dega) ? numit : dega;
                for (int j = 1; j < jMax; j++) {
                    pix = pix.add(1.0, piList.get(numit - j).mult(Awork[j + 1]));
                }
                pix = pix.mult(invbarA1);
                sumpi += pix.elementSum();
                piList.add(pix);
                numit++;
            }

            Matrix result = Matrix.zeros(1, (K + 1) * m);
            int kCap = Math.min(piList.size(), K + 1);
            for (int k = 0; k < kCap; k++) {
                for (int j = 0; j < m; j++) {
                    result.set(0, k * m + j, piList.get(k).get(0, j));
                }
            }
            return result;
        }
    }

    public static Matrix mg1StationaryDistr(List<Matrix> A) {
        return mg1StationaryDistr(A, null, null, 500, 1e-14);
    }

    public static Matrix mg1StationaryDistr(List<Matrix> A, List<Matrix> B) {
        return mg1StationaryDistr(A, B, null, 500, 1e-14);
    }

    public static Matrix mg1StationaryDistr(List<Matrix> A, List<Matrix> B, Matrix G) {
        return mg1StationaryDistr(A, B, G, 500, 1e-14);
    }

    public static Matrix mg1StationaryDistr(List<Matrix> A, List<Matrix> B, Matrix G, int K) {
        return mg1StationaryDistr(A, B, G, K, 1e-14);
    }

    /**
     * Overload accepting Matrix[] for A and B.
     */
    public static Matrix mg1StationaryDistr(Matrix[] A, Matrix[] B, Matrix G, int K, double prec) {
        return mg1StationaryDistr(
                Arrays.asList(A),
                B != null ? Arrays.asList(B) : null,
                G, K, prec);
    }

    public static Matrix mg1StationaryDistr(Matrix[] A) {
        return mg1StationaryDistr(A, null, null, 500, 1e-14);
    }

    public static Matrix mg1StationaryDistr(Matrix[] A, Matrix[] B) {
        return mg1StationaryDistr(A, B, null, 500, 1e-14);
    }

    public static Matrix mg1StationaryDistr(Matrix[] A, Matrix[] B, Matrix G) {
        return mg1StationaryDistr(A, B, G, 500, 1e-14);
    }

    public static Matrix mg1StationaryDistr(Matrix[] A, Matrix[] B, Matrix G, int K) {
        return mg1StationaryDistr(A, B, G, K, 1e-14);
    }
}
