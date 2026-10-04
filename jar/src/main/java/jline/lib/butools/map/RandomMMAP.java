/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import jline.lib.butools.mc.CTMCSolve;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class RandomMMAP {
    private RandomMMAP() {}

    public static MatrixCell randomMMAP(int order, int types) {
        return randomMMAP(order, types, 1.0, 0, 1000, 1e-7, new Random());
    }

    public static MatrixCell randomMMAP(int order, int types, double mean) {
        return randomMMAP(order, types, mean, 0, 1000, 1e-7, new Random());
    }

    public static MatrixCell randomMMAP(int order, int types, double mean, int zeroEntries) {
        return randomMMAP(order, types, mean, zeroEntries, 1000, 1e-7, new Random());
    }

    /**
     * Returns a random continuous marked Markovian arrival process.
     */
    public static MatrixCell randomMMAP(int order, int types, double mean, int zeroEntries,
                                        int maxTrials, double prec, Random random) {
        if (types < 1) {
            throw new IllegalArgumentException("RandomMMAP: 'types' must be positive integer!");
        }

        if (zeroEntries > (order + 1) * (order - 1) + types * (order * order - 1)) {
            throw new IllegalArgumentException("RandomMMAP: Too many zero entries requested!");
        }

        List<int[]> zeroDistr = allZeroDistributionsMMAP(order, zeroEntries);

        int trials = 0;
        while (trials < maxTrials) {
            List<Integer> indices = new ArrayList<Integer>();
            for (int i = 0; i < zeroDistr.size(); i++) indices.add(i);
            Collections.shuffle(indices, random);

            for (int zdix : indices) {
                int[] zDistr = zeroDistr.get(zdix);

                boolean bad = false;
                for (int z : zDistr) {
                    if (z >= (types + 1) * order - 1) {
                        bad = true;
                        break;
                    }
                }
                if (bad) {
                    trials++;
                    continue;
                }

                Matrix B = new Matrix(order, (types + 1) * order);
                for (int i = 0; i < order; i++) {
                    List<Integer> rp = new ArrayList<Integer>();
                    for (int p = 0; p < (types + 1) * order - 1; p++) rp.add(p);
                    Collections.shuffle(rp, random);
                    double[] a = new double[(types + 1) * order - 1];
                    for (int j = 0; j < (types + 1) * order - 1 - zDistr[i]; j++) {
                        a[rp.get(j)] = random.nextDouble();
                    }
                    for (int j = 0; j < i; j++) {
                        B.set(i, j, a[j]);
                    }
                    for (int j = i + 1; j < (types + 1) * order; j++) {
                        B.set(i, j, a[j - 1]);
                    }
                }

                MatrixCell D = new MatrixCell(types + 1);
                double[] sc = new double[order];
                for (int k = 0; k <= types; k++) {
                    D.set(k, new Matrix(order, order));
                    for (int i = 0; i < order; i++) {
                        for (int j = 0; j < order; j++) {
                            D.get(k).set(i, j, B.get(i, k * order + j));
                            sc[i] += D.get(k).get(i, j);
                        }
                    }
                }

                boolean anyZero = false;
                for (double s : sc) if (s == 0.0) { anyZero = true; break; }
                if (anyZero) continue;

                for (int i = 0; i < order; i++) {
                    D.get(0).set(i, i, -sc[i]);
                }

                Matrix sumD = new Matrix(order, order);
                for (int k = 0; k <= types; k++) {
                    sumD = sumD.add(D.get(k));
                }

                if (sumD.rank() == order - 1) {
                    Matrix alpha = CTMCSolve.ctmcSolve(sumD);
                    double minAlpha = Double.POSITIVE_INFINITY;
                    for (int j = 0; j < order; j++) {
                        minAlpha = Math.min(minAlpha, Math.abs(alpha.get(0, j)));
                    }

                    if (minAlpha > prec) {
                        boolean fullZero = false;
                        for (int k = 1; k <= types; k++) {
                            boolean allZero = true;
                            outer:
                            for (int i = 0; i < order; i++) {
                                for (int j = 0; j < order; j++) {
                                    if (D.get(k).get(i, j) != 0.0) {
                                        allZero = false;
                                        break outer;
                                    }
                                }
                            }
                            if (allZero) {
                                fullZero = true;
                                break;
                            }
                        }

                        if (!fullZero) {
                            double[] m = MarginalMomentsFromMMAP.marginalMomentsFromMMAP(D, 1, prec);
                            double scaleFactor = m[0] / mean;
                            for (int k = 0; k <= types; k++) {
                                for (int i = 0; i < order; i++) {
                                    for (int j = 0; j < order; j++) {
                                        D.get(k).set(i, j, D.get(k).get(i, j) * scaleFactor);
                                    }
                                }
                            }

                            if (CheckMMAPRepresentation.checkMMAPRepresentation(D, prec)) {
                                return D;
                            }
                        }
                    }
                }
                trials++;
            }
        }

        throw new IllegalStateException("No feasible random MMAP found!");
    }

    private static List<int[]> allZeroDistributionsMMAP(int states, int zeros) {
        if (states == 1) {
            List<int[]> result = new ArrayList<int[]>();
            result.add(new int[]{zeros});
            return result;
        }
        List<int[]> result = new ArrayList<int[]>();
        for (int iz = 0; iz <= zeros; iz++) {
            List<int[]> subDistributions = allZeroDistributionsMMAP(states - 1, zeros - iz);
            for (int[] sub : subDistributions) {
                int[] combined = new int[sub.length + 1];
                System.arraycopy(sub, 0, combined, 0, sub.length);
                combined[sub.length] = iz;
                java.util.Arrays.sort(combined);
                boolean exists = false;
                for (int[] r : result) {
                    if (java.util.Arrays.equals(r, combined)) {
                        exists = true;
                        break;
                    }
                }
                if (!exists) {
                    result.add(combined);
                }
            }
        }
        return result;
    }
}
