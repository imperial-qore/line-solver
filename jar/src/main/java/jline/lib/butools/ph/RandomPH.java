/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import jline.lib.butools.mc.CTMCSolve;
import jline.util.matrix.Matrix;

public final class RandomPH {
    private RandomPH() {}

    public static PHRepresentation randomPH(int order) {
        return randomPH(order, 1.0, 0, 1000, 1e-7, new Random());
    }

    public static PHRepresentation randomPH(int order, double mean) {
        return randomPH(order, mean, 0, 1000, 1e-7, new Random());
    }

    public static PHRepresentation randomPH(int order, double mean, int zeroEntries) {
        return randomPH(order, mean, zeroEntries, 1000, 1e-7, new Random());
    }

    public static PHRepresentation randomPH(int order, double mean, int zeroEntries,
                                            int maxTrials, double prec, Random random) {
        if (zeroEntries > (order + 1) * (order - 1)) {
            throw new IllegalArgumentException("RandomPH: Too many zero entries requested!");
        }

        List<int[]> zeroDistr = allZeroDistr(order, zeroEntries);

        int trials = 0;
        while (trials < maxTrials) {
            List<Integer> indices = new ArrayList<Integer>();
            for (int i = 0; i < zeroDistr.size(); i++) indices.add(i);
            Collections.shuffle(indices, random);

            for (int zdix : indices) {
                int[] zDistr = zeroDistr.get(zdix);

                Matrix B = Matrix.zeros(order, order + 2);
                for (int i = 0; i < order; i++) {
                    List<Integer> rp = new ArrayList<Integer>();
                    for (int p = 0; p < order + 1; p++) rp.add(p);
                    Collections.shuffle(rp, random);
                    double[] a = new double[order + 1];
                    for (int j = 0; j < order + 1 - zDistr[i]; j++) {
                        a[rp.get(j)] = random.nextDouble();
                    }
                    for (int j = 0; j < i; j++) {
                        B.set(i, j, a[j]);
                    }
                    for (int j = i; j < order + 1; j++) {
                        B.set(i, j + 1, a[j]);
                    }
                }

                Matrix A = Matrix.zeros(order, order);
                for (int i = 0; i < order; i++) {
                    for (int j = 0; j < order; j++) {
                        A.set(i, j, B.get(i, j));
                    }
                }
                double[] aVec = new double[order];
                double[] alphaVec = new double[order];
                for (int i = 0; i < order; i++) {
                    aVec[i] = B.get(i, order + 1);
                    alphaVec[i] = B.get(i, order);
                }

                for (int i = 0; i < order; i++) {
                    double rowSum = 0.0;
                    for (int j = 0; j < order; j++) {
                        rowSum += A.get(i, j);
                    }
                    A.set(i, i, A.get(i, i) - rowSum - aVec[i]);
                }

                boolean allAZero = true;
                boolean allAlphaZero = true;
                boolean allAVecZero = true;
                for (int i = 0; i < order; i++) {
                    for (int j = 0; j < order; j++) {
                        if (A.get(i, j) != 0.0) allAZero = false;
                    }
                    if (alphaVec[i] != 0.0) allAlphaZero = false;
                    if (aVec[i] != 0.0) allAVecZero = false;
                }

                if (allAZero || allAlphaZero || allAVecZero) {
                    continue;
                }

                double alphaSum = 0.0;
                for (int i = 0; i < order; i++) {
                    alphaSum += alphaVec[i];
                }
                Matrix alpha = new Matrix(1, order);
                for (int i = 0; i < order; i++) {
                    alpha.set(0, i, alphaVec[i] / alphaSum);
                }

                Matrix D = Matrix.zeros(order, order);
                for (int i = 0; i < order; i++) {
                    for (int j = 0; j < order; j++) {
                        D.set(i, j, A.get(i, j) + aVec[i] * alpha.get(0, j));
                    }
                }

                if (D.rank() == order - 1) {
                    Matrix pi = CTMCSolve.ctmcSolve(D);
                    double minPi = Double.MAX_VALUE;
                    for (int i = 0; i < pi.length(); i++) {
                        if (Math.abs(pi.get(i)) < minPi) minPi = Math.abs(pi.get(i));
                    }

                    if (minPi > prec) {
                        double[] moms = MomentsFromME.momentsFromME(alpha, A, 1);
                        double scaleFactor = moms[0] / mean;
                        Matrix scaledA = A.scale(scaleFactor);
                        return new PHRepresentation(alpha, scaledA);
                    }
                }
                trials++;
            }
        }

        throw new IllegalArgumentException("No feasible random PH found with such many zero entries!");
    }

    private static List<int[]> allZeroDistr(int states, int zeros) {
        List<int[]> result = new ArrayList<int[]>();
        if (states == 1) {
            result.add(new int[]{zeros});
            return result;
        }

        for (int iz = 0; iz <= zeros; iz++) {
            List<int[]> subResults = allZeroDistr(states - 1, zeros - iz);
            for (int[] subResult : subResults) {
                int[] xt = new int[subResult.length + 1];
                System.arraycopy(subResult, 0, xt, 0, subResult.length);
                xt[subResult.length] = iz;
                java.util.Arrays.sort(xt);
                boolean found = false;
                for (int[] existing : result) {
                    if (java.util.Arrays.equals(existing, xt)) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    result.add(xt);
                }
            }
        }
        return result;
    }
}
