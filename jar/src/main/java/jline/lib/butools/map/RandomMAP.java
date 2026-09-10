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

import jline.util.Pair;

import jline.lib.butools.mc.CTMCSolve;
import jline.util.matrix.Matrix;

public final class RandomMAP {
    private RandomMAP() {}

    /**
     * Returns a random Markovian arrival process with given mean value.
     */
    public static MAPRepresentation randomMAP(int order, double mean, int zeroEntries, int maxTrials,
                                               double prec, Random random) {
        if (zeroEntries > 2 * order * (order - 1)) {
            throw new IllegalArgumentException("RandomMAP: Too many zero entries requested!");
        }

        List<int[]> zeroDistr = allZeroDistr(order, zeroEntries);

        int trials = 0;
        while (trials < maxTrials) {
            List<int[]> zdixList = new ArrayList<int[]>(zeroDistr);
            Collections.shuffle(zdixList, random);

            for (int z = 0; z < zdixList.size(); z++) {
                Matrix D0 = Matrix.zeros(order, order);
                Matrix D1 = Matrix.zeros(order, order);

                int totalEntries = 2 * order * order;
                int nonZeroEntries = totalEntries - zeroEntries;

                List<Pair<Integer, Integer>> allPositions = new ArrayList<Pair<Integer, Integer>>();
                for (int i = 0; i < order; i++) {
                    for (int j = 0; j < order; j++) {
                        allPositions.add(new Pair<Integer, Integer>(Integer.valueOf(i), Integer.valueOf(j)));
                        allPositions.add(new Pair<Integer, Integer>(Integer.valueOf(i), Integer.valueOf(j + order)));
                    }
                }

                Collections.shuffle(allPositions, random);
                int limit = Math.min(nonZeroEntries, allPositions.size());
                for (int idx = 0; idx < limit; idx++) {
                    Pair<Integer, Integer> pos = allPositions.get(idx);
                    int i = pos.getFirst().intValue();
                    int j = pos.getSecond().intValue();
                    double value = random.nextDouble();
                    if (j < order) {
                        if (i != j) {
                            D0.set(i, j, value);
                        }
                    } else {
                        D1.set(i, j - order, value);
                    }
                }

                for (int i = 0; i < order; i++) {
                    double rowSum = 0.0;
                    for (int j = 0; j < order; j++) {
                        if (i != j) rowSum += D0.get(i, j);
                        rowSum += D1.get(i, j);
                    }
                    D0.set(i, i, -rowSum);
                }

                boolean allD0Zero = true;
                boolean allD1Zero = true;
                for (int i = 0; i < order; i++) {
                    for (int j = 0; j < order; j++) {
                        if (D0.get(i, j) != 0.0) allD0Zero = false;
                        if (D1.get(i, j) != 0.0) allD1Zero = false;
                    }
                }

                if (allD0Zero || allD1Zero) {
                    continue;
                }

                Matrix D = D0.add(D1);
                if (D.rank() == order - 1) {
                    Matrix pi = CTMCSolve.ctmcSolve(D);
                    double minPi = Double.MAX_VALUE;
                    for (int i = 0; i < pi.length(); i++) {
                        if (Math.abs(pi.get(i)) < minPi) minPi = Math.abs(pi.get(i));
                    }

                    if (minPi > prec) {
                        double[] moms = MarginalMomentsFromMAP.marginalMomentsFromMAP(D0, D1, 1);
                        double scaleFactor = moms[0] / mean;
                        return new MAPRepresentation(D0.scale(scaleFactor), D1.scale(scaleFactor));
                    }
                }
                trials++;
            }
        }

        throw new IllegalArgumentException("No feasible random MAP found!");
    }

    public static MAPRepresentation randomMAP(int order, double mean, int zeroEntries, int maxTrials, double prec) {
        return randomMAP(order, mean, zeroEntries, maxTrials, prec, new Random());
    }

    public static MAPRepresentation randomMAP(int order, double mean, int zeroEntries, int maxTrials) {
        return randomMAP(order, mean, zeroEntries, maxTrials, 1e-7, new Random());
    }

    public static MAPRepresentation randomMAP(int order, double mean, int zeroEntries) {
        return randomMAP(order, mean, zeroEntries, 1000, 1e-7, new Random());
    }

    public static MAPRepresentation randomMAP(int order, double mean) {
        return randomMAP(order, mean, 0, 1000, 1e-7, new Random());
    }

    public static MAPRepresentation randomMAP(int order) {
        return randomMAP(order, 1.0, 0, 1000, 1e-7, new Random());
    }

    /**
     * Helper function to generate all zero distributions.
     */
    private static List<int[]> allZeroDistr(int states, int zeros) {
        if (states == 1) {
            return Collections.singletonList(new int[]{zeros});
        }

        List<int[]> result = new ArrayList<int[]>();
        for (int iz = 0; iz <= zeros; iz++) {
            List<int[]> subResults = allZeroDistr(states - 1, zeros - iz);
            for (int s = 0; s < subResults.size(); s++) {
                int[] subResult = subResults.get(s);
                int[] xt = new int[subResult.length + 1];
                System.arraycopy(subResult, 0, xt, 0, subResult.length);
                xt[subResult.length] = iz;
                java.util.Arrays.sort(xt);
                boolean found = false;
                for (int e = 0; e < result.size(); e++) {
                    int[] existing = result.get(e);
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
