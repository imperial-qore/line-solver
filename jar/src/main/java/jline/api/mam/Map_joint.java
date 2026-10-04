/**
 * @file Markovian Arrival Process joint moment analysis
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_joint {
    private Map_joint() {}

    /**
     * Computes the joint moments of a Markovian Arrival Process (MAP).
     */
    public static double map_joint(MatrixCell MAP, int[] a, int[] i) {
        return map_joint(MAP.get(0), MAP.get(1), a, i);
    }

    /**
     * Computes the joint moments of a Markovian Arrival Process (MAP).
     */
    public static double map_joint(Matrix D0, Matrix D1, int[] a, int[] i) {
        if (a.length != i.length) {
            throw new IllegalArgumentException("Vectors a and i must have the same length");
        }

        int[] aCumSum = new int[a.length];
        aCumSum[0] = a[0];
        for (int k = 1; k < a.length; k++) {
            aCumSum[k] = aCumSum[k - 1] + a[k];
        }

        Matrix minusD0 = D0.scale(-1.0);
        Matrix invD0 = minusD0.inv();
        Matrix P = invD0.mult(D1);

        Matrix JM = Matrix.eye(D0.getNumRows());
        int K = a.length;

        for (int k = 0; k < (K - 1); k++) {
            double factorial_i_k = (double) factorial(i[k]);
            Matrix invD0_power = matrixPower(invD0, i[k]);
            Matrix P_power = matrixPower(P, (aCumSum[k + 1] - aCumSum[k]));

            JM = JM.mult(invD0_power.scale(factorial_i_k)).mult(P_power);
        }

        Matrix pie = Map_pie.map_pie(D0, D1);
        double factorial_i_K = (double) factorial(i[K - 1]);
        Matrix invD0_power_final = matrixPower(invD0, i[K - 1]);
        Matrix ones = Matrix.ones(P.getNumRows(), 1);

        Matrix finalResult = pie.mult(JM)
                .mult(invD0_power_final.scale(factorial_i_K))
                .mult(ones);

        return finalResult.get(0, 0);
    }

    private static long factorial(int n) {
        if (n < 0) throw new IllegalArgumentException("Factorial is not defined for negative numbers");
        if (n == 0 || n == 1) return 1;
        long result = 1L;
        for (int v = 2; v <= n; v++) {
            result *= v;
        }
        return result;
    }

    private static Matrix matrixPower(Matrix matrix, int power) {
        if (power < 0) throw new IllegalArgumentException("Negative powers not supported");
        if (power == 0) return Matrix.eye(matrix.getNumRows());
        if (power == 1) return matrix.copy();

        Matrix result = Matrix.eye(matrix.getNumRows());
        Matrix base = matrix.copy();
        int exp = power;

        while (exp > 0) {
            if (exp % 2 == 1) {
                result = result.mult(base);
            }
            base = base.mult(base);
            exp /= 2;
        }
        return result;
    }
}
