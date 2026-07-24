/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import jline.lib.butools.mc.DTMCSolve;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class RandomDMMAP {
    private RandomDMMAP() {}

    /**
     * Returns a random discrete marked Markovian arrival process.
     */
    public static MatrixCell randomDMMAP(int order, int types, double mean, int zeroEntries,
                                         int maxTrials, double prec, Random random) {
        if (types < 1) {
            throw new IllegalArgumentException("RandomDMMAP: 'types' must be positive integer!");
        }

        if (zeroEntries > (order + 1) * (order - 1) + types * (order * order - 1)) {
            throw new IllegalArgumentException("RandomDMMAP: Too many zero entries requested!");
        }

        // Generate all possible zero distributions among rows
        List<int[]> zeroDistr = allZeroDistributions(order, zeroEntries);

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
                    for (int j = 0; j < (types + 1) * order - 1; j++) rp.add(j);
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
                for (int kk = 0; kk <= types; kk++) {
                    D.set(kk, new Matrix(order, order));
                    for (int i = 0; i < order; i++) {
                        for (int j = 0; j < order; j++) {
                            D.get(kk).set(i, j, B.get(i, kk * order + j));
                            sc[i] += D.get(kk).get(i, j);
                        }
                    }
                }

                boolean anyZero = false;
                for (double v : sc) {
                    if (v == 0.0) {
                        anyZero = true;
                        break;
                    }
                }
                if (anyZero) continue;

                for (int kk = 0; kk <= types; kk++) {
                    for (int i = 0; i < order; i++) {
                        for (int j = 0; j < order; j++) {
                            D.get(kk).set(i, j, D.get(kk).get(i, j) / sc[i]);
                        }
                    }
                }

                Matrix sumD = new Matrix(order, order);
                for (int kk = 0; kk <= types; kk++) {
                    sumD = sumD.add(D.get(kk));
                }

                Matrix I = Matrix.eye(order);
                if (D.get(0).rank() == order && I.sub(sumD).rank() == order - 1) {
                    Matrix alpha = DTMCSolve.dtmcSolve(sumD);
                    double minAlpha = Double.POSITIVE_INFINITY;
                    for (int j = 0; j < order; j++) {
                        minAlpha = Math.min(minAlpha, Math.abs(alpha.get(0, j)));
                    }

                    if (minAlpha > prec) {
                        boolean fullZero = false;
                        for (int kk = 0; kk <= types; kk++) {
                            boolean allZero = true;
                            outer:
                            for (int i = 0; i < order; i++) {
                                for (int j = 0; j < order; j++) {
                                    if (D.get(kk).get(i, j) != 0.0) {
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
                            double[] d = new double[order];
                            for (int idx = 0; idx < order; idx++) d[idx] = random.nextDouble();

                            MatrixCell Dv = new MatrixCell(types + 1);
                            for (int kk = 0; kk <= types; kk++) {
                                Dv.set(kk, new Matrix(order, order));
                                for (int i = 0; i < order; i++) {
                                    for (int j = 0; j < order; j++) {
                                        Dv.get(kk).set(i, j, (1.0 - d[i]) * D.get(kk).get(i, j));
                                    }
                                }
                            }
                            for (int i = 0; i < order; i++) {
                                Dv.get(0).set(i, i, Dv.get(0).get(i, i) + d[i]);
                            }

                            double[] m = MarginalMomentsFromDMMAP.marginalMomentsFromDMMAP(Dv, 1, prec);

                            for (int i = 0; i < order; i++) {
                                d[i] = 1.0 - (1.0 - d[i]) * m[0] / mean;
                            }

                            for (int kk = 0; kk <= types; kk++) {
                                for (int i = 0; i < order; i++) {
                                    for (int j = 0; j < order; j++) {
                                        D.get(kk).set(i, j, (1.0 - d[i]) * D.get(kk).get(i, j));
                                    }
                                }
                            }
                            for (int i = 0; i < order; i++) {
                                D.get(0).set(i, i, D.get(0).get(i, i) + d[i]);
                            }

                            if (CheckDMMAPRepresentation.checkDMMAPRepresentation(D, prec)) {
                                return D;
                            }
                        }
                    }
                }
                trials++;
            }
        }

        throw new IllegalStateException("No feasible random DMMAP found!");
    }

    public static MatrixCell randomDMMAP(int order, int types, double mean, int zeroEntries,
                                         int maxTrials, double prec) {
        return randomDMMAP(order, types, mean, zeroEntries, maxTrials, prec, new Random());
    }

    public static MatrixCell randomDMMAP(int order, int types, double mean, int zeroEntries, int maxTrials) {
        return randomDMMAP(order, types, mean, zeroEntries, maxTrials, 1e-7, new Random());
    }

    public static MatrixCell randomDMMAP(int order, int types, double mean, int zeroEntries) {
        return randomDMMAP(order, types, mean, zeroEntries, 1000, 1e-7, new Random());
    }

    public static MatrixCell randomDMMAP(int order, int types, double mean) {
        return randomDMMAP(order, types, mean, 0, 1000, 1e-7, new Random());
    }

    public static MatrixCell randomDMMAP(int order, int types) {
        return randomDMMAP(order, types, 10.0, 0, 1000, 1e-7, new Random());
    }

    private static List<int[]> allZeroDistributions(int states, int zeros) {
        if (states == 1) {
            return Collections.singletonList(new int[]{zeros});
        }
        List<int[]> result = new ArrayList<int[]>();
        for (int iz = 0; iz <= zeros; iz++) {
            List<int[]> subDistributions = allZeroDistributions(states - 1, zeros - iz);
            for (int[] sub : subDistributions) {
                int[] combined = new int[sub.length + 1];
                System.arraycopy(sub, 0, combined, 0, sub.length);
                combined[sub.length] = iz;
                java.util.Arrays.sort(combined);
                boolean exists = false;
                for (int[] e : result) {
                    if (java.util.Arrays.equals(e, combined)) {
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
