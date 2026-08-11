/**
 * @file Random sampling with relative probabilities
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.Random;

import jline.util.matrix.Matrix;

public final class Randp {
    private Randp() {}

    private static final Random RANDOM = new Random();

    /**
     * Pick random values with relative probability.
     */
    public static Matrix randp(double[] P, int rows, int cols) {
        for (double p : P) {
            if (p < 0.0) {
                throw new IllegalArgumentException("All probabilities should be 0 or larger");
            }
        }
        Matrix result = Matrix.zeros(rows, cols);
        if (P.length == 0) return result;
        double totalSum = 0.0;
        for (double p : P) totalSum += p;
        if (totalSum == 0.0) return result;

        double[] cumsum = new double[P.length + 1];
        cumsum[0] = 0.0;
        for (int i = 0; i < P.length; i++) {
            cumsum[i + 1] = cumsum[i] + P[i] / totalSum;
        }

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                double u = RANDOM.nextDouble();
                int bin = 1;
                for (int k = 1; k < cumsum.length; k++) {
                    if (u <= cumsum[k]) {
                        bin = k;
                        break;
                    }
                }
                result.set(i, j, (double) bin);
            }
        }
        return result;
    }

    public static Matrix randp(double[] P, int rows) {
        return randp(P, rows, rows);
    }

    public static Matrix randp(Matrix P, int rows, int cols) {
        double[] probArray = new double[P.getNumRows() * P.getNumCols()];
        int idx = 0;
        for (int j = 0; j < P.getNumCols(); j++) {
            for (int i = 0; i < P.getNumRows(); i++) {
                probArray[idx++] = P.get(i, j);
            }
        }
        return randp(probArray, rows, cols);
    }

    public static Matrix randp(Matrix P, int rows) {
        return randp(P, rows, rows);
    }

    public static int randp(double[] P) {
        Matrix result = randp(P, 1, 1);
        return (int) result.get(0, 0);
    }

    public static int randp(Matrix P) {
        Matrix result = randp(P, 1, 1);
        return (int) result.get(0, 0);
    }
}
