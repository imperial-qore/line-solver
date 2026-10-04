/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import jline.util.matrix.Matrix;

public final class RandomDPH {
    private RandomDPH() {}

    /**
     * Returns a random discrete phase-type distribution with a given mean value.
     */
    public static MGFromMoments.MGRepresentation randomDPH(int order, double mean, int zeroEntries,
                                                              int maxTrials, double prec, Random random) {
        if (zeroEntries > (order + 1) * (order - 1)) {
            throw new IllegalArgumentException("RandomDPH: Too many zero entries requested!");
        }

        // Generate all possible distributions of zero entries among rows
        List<int[]> zeroDistr = allZeroDistr(order, zeroEntries);

        int trials = 0;
        while (trials < maxTrials) {
            // Randomly select a configuration
            List<Integer> zdixList = new ArrayList<Integer>();
            for (int idx = 0; idx < zeroDistr.size(); idx++) zdixList.add(idx);
            Collections.shuffle(zdixList, random);

            for (int zdix : zdixList) {
                int[] zDistr = zeroDistr.get(zdix);

                // Create B matrix (order x (order+2))
                Matrix B = Matrix.zeros(order, order + 2);

                for (int i = 0; i < order; i++) {
                    List<Integer> rp = new ArrayList<Integer>();
                    for (int j = 0; j < order + 1; j++) rp.add(j);
                    Collections.shuffle(rp, random);
                    double[] a = new double[order + 1];

                    for (int j = 0; j < order + 1 - zDistr[i]; j++) {
                        a[rp.get(j)] = random.nextDouble();
                    }

                    // B(i, 1:i-1) = a(1:i-1)
                    for (int j = 0; j < i; j++) {
                        B.set(i, j, a[j]);
                    }
                    // B(i, i+1:end) = a(i:end)
                    for (int j = i; j < order + 1; j++) {
                        B.set(i, j + 1, a[j]);
                    }
                }

                // Construct DPH parameters
                Matrix A = Matrix.zeros(order, order);
                for (int i = 0; i < order; i++) {
                    for (int j = 0; j < order; j++) {
                        A.set(i, j, B.get(i, j));
                    }
                }

                Matrix aVec = new Matrix(order, 1);
                for (int i = 0; i < order; i++) {
                    aVec.set(i, 0, B.get(i, order + 1));
                }

                // Scale rows
                double[] sc = new double[order];
                boolean hasZeroRow = false;
                for (int i = 0; i < order; i++) {
                    double rowSum = 0.0;
                    for (int j = 0; j < order; j++) {
                        rowSum += A.get(i, j);
                    }
                    rowSum += aVec.get(i, 0);
                    sc[i] = rowSum;
                    if (rowSum == 0.0) {
                        hasZeroRow = true;
                        break;
                    }
                }

                if (hasZeroRow) continue;

                for (int i = 0; i < order; i++) {
                    for (int j = 0; j < order; j++) {
                        A.set(i, j, A.get(i, j) / sc[i]);
                    }
                    aVec.set(i, 0, aVec.get(i, 0) / sc[i]);
                }

                // Extract alpha
                Matrix alpha = new Matrix(1, order);
                for (int i = 0; i < order; i++) {
                    alpha.set(0, i, B.get(i, order));
                }

                // Check for all-zero matrices
                if (A.elementMax() == 0.0 || alpha.elementMax() == 0.0 || aVec.elementMax() == 0.0) {
                    continue;
                }

                // Normalize alpha
                double alphaSum = alpha.elementSum();
                if (alphaSum == 0.0) continue;
                alpha = alpha.scale(1.0 / alphaSum);

                // Check irreducibility: rank(I - A) should be order
                Matrix I = Matrix.eye(order);
                Matrix IminusA = I.sub(A);
                int rank = IminusA.rank();

                if (rank == order) {
                    // Check if alpha * inv(I - A) has all positive elements
                    Matrix invIminusA = IminusA.inv();
                    Matrix alphaInv = alpha.mult(invIminusA);
                    boolean allPositive = true;
                    for (int i = 0; i < order; i++) {
                        if (Math.abs(alphaInv.get(0, i)) <= prec) {
                            allPositive = false;
                            break;
                        }
                    }

                    if (allPositive) {
                        // Scale diagonals to achieve target mean
                        double[] d = new double[order];
                        for (int i = 0; i < order; i++) d[i] = random.nextDouble();

                        // Compute current mean with diagonal modification
                        Matrix scaledA = Matrix.zeros(order, order);
                        for (int i = 0; i < order; i++) {
                            for (int j = 0; j < order; j++) {
                                if (i == j) {
                                    scaledA.set(i, j, d[i]);
                                } else {
                                    scaledA.set(i, j, (1 - d[i]) * A.get(i, j));
                                }
                            }
                        }

                        double[] testMoms = MomentsFromDPH.momentsFromDPH(alpha, scaledA, 1);
                        double currentMean = testMoms[0];

                        // Scale to target mean
                        for (int i = 0; i < order; i++) {
                            d[i] = 1 - (1 - d[i]) * currentMean / mean;
                        }

                        // Rebuild A with scaled diagonals
                        for (int i = 0; i < order; i++) {
                            for (int j = 0; j < order; j++) {
                                if (i == j) {
                                    A.set(i, j, d[i]);
                                } else {
                                    A.set(i, j, (1 - d[i]) * A.get(i, j));
                                }
                            }
                        }

                        if (CheckDPHRepresentation.checkDPHRepresentation(alpha, A, prec)) {
                            return new MGFromMoments.MGRepresentation(alpha, A);
                        }
                    }
                }

                trials++;
            }
        }

        throw new IllegalArgumentException("RandomDPH: No feasible random DPH found! Try increasing mean or maxTrials.");
    }

    public static MGFromMoments.MGRepresentation randomDPH(int order, double mean, int zeroEntries,
                                                              int maxTrials, double prec) {
        return randomDPH(order, mean, zeroEntries, maxTrials, prec, new Random());
    }

    public static MGFromMoments.MGRepresentation randomDPH(int order, double mean, int zeroEntries, int maxTrials) {
        return randomDPH(order, mean, zeroEntries, maxTrials, 1e-7, new Random());
    }

    public static MGFromMoments.MGRepresentation randomDPH(int order, double mean, int zeroEntries) {
        return randomDPH(order, mean, zeroEntries, 1000, 1e-7, new Random());
    }

    public static MGFromMoments.MGRepresentation randomDPH(int order, double mean) {
        return randomDPH(order, mean, 0, 1000, 1e-7, new Random());
    }

    public static MGFromMoments.MGRepresentation randomDPH(int order) {
        return randomDPH(order, 10.0, 0, 1000, 1e-7, new Random());
    }

    /**
     * Generate all possible distributions of zero entries among rows.
     */
    private static List<int[]> allZeroDistr(int states, int zeros) {
        if (states == 1) {
            return Collections.singletonList(new int[]{zeros});
        }

        List<int[]> result = new ArrayList<int[]>();
        for (int iz = 0; iz <= zeros; iz++) {
            List<int[]> subDistrs = allZeroDistr(states - 1, zeros - iz);
            for (int[] subDistr : subDistrs) {
                int[] combined = new int[subDistr.length + 1];
                System.arraycopy(subDistr, 0, combined, 0, subDistr.length);
                combined[subDistr.length] = iz;
                Arrays.sort(combined);
                boolean found = false;
                for (int[] existing : result) {
                    if (Arrays.equals(existing, combined)) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    result.add(combined);
                }
            }
        }
        return result;
    }
}
