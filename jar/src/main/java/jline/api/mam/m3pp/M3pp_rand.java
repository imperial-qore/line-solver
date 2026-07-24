/**
 * @file M3PP random process generator
 *
 * @since LINE 3.0
 */
package jline.api.mam.m3pp;

import java.util.Random;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class M3pp_rand {
    private M3pp_rand() {}

    public static MatrixCell m3pp_rand(int order, int classes) {
        return m3pp_rand(order, classes, null);
    }

    /**
     * Generates a random M3PP with specified order and number of classes.
     */
    public static MatrixCell m3pp_rand(int order, int classes, Long seed) {
        Random random = (seed != null) ? new Random(seed) : new Random();

        MatrixCell mmpp = mmpp_rand(order, random);

        MatrixCell mmap = new MatrixCell(classes + 2);
        mmap.set(0, mmpp.get(0));
        mmap.set(1, mmpp.get(1));

        for (int i = 0; i < order; i++) {
            double[] p = new double[classes];
            double pSum = 0.0;
            for (int j = 0; j < classes; j++) {
                p[j] = random.nextDouble();
                pSum += p[j];
            }

            for (int j = 0; j < p.length; j++) {
                p[j] /= pSum;
            }

            for (int c = 0; c < classes; c++) {
                if (mmap.get(2 + c) == null) {
                    mmap.set(2 + c, new Matrix(order, order));
                }

                Matrix classMatrix = mmap.get(2 + c);
                Matrix baseMatrix = mmpp.get(1);

                for (int row = 0; row < order; row++) {
                    for (int col = 0; col < order; col++) {
                        classMatrix.set(row, col, baseMatrix.get(row, col) * p[c]);
                    }
                }
            }
        }

        return mmap;
    }

    private static MatrixCell mmpp_rand(int order, Random random) {
        Matrix D0 = new Matrix(order, order);
        Matrix D1 = new Matrix(order, order);

        for (int i = 0; i < order; i++) {
            for (int j = 0; j < order; j++) {
                if (i != j) {
                    D0.set(i, j, random.nextDouble() * 2.0);
                }
            }
        }

        for (int i = 0; i < order; i++) {
            D1.set(i, i, random.nextDouble() * 5.0 + 0.1);
        }

        // Normalize so that D0's diagonal makes (D0 + D1) a proper generator;
        // matches MATLAB mmpp_rand, which ends with map_normalize. Without it
        // the returned process is not a valid MAP (positive row sums).
        return jline.api.mam.Map_normalize.map_normalize(D0, D1);
    }

    public static MatrixCell m3pp_rand_targeted(int order, int classes) {
        return m3pp_rand_targeted(order, classes, 1.0, 2.0, null);
    }

    public static MatrixCell m3pp_rand_targeted(int order, int classes, double targetRate,
                                                double targetSCV, Long seed) {
        Random random = (seed != null) ? new Random(seed) : new Random();
        MatrixCell bestMmap = null;
        double bestError = Double.MAX_VALUE;

        for (int iter = 0; iter < 20; iter++) {
            MatrixCell candidate = m3pp_rand(order, classes, random.nextLong());

            double candidateRate = computeArrivalRate(candidate);
            double candidateSCV = computeSCV(candidate);

            double rateError = Math.abs(candidateRate - targetRate) / targetRate;
            double scvError = Math.abs(candidateSCV - targetSCV) / targetSCV;
            double totalError = rateError + scvError;

            if (totalError < bestError) {
                bestError = totalError;
                bestMmap = candidate;
            }
        }

        return bestMmap != null ? bestMmap : m3pp_rand(order, classes, seed);
    }

    private static double computeArrivalRate(MatrixCell mmap) {
        Matrix D1 = mmap.get(1);
        return D1.elementSum();
    }

    private static double computeSCV(MatrixCell mmap) {
        return 1.5;
    }
}
