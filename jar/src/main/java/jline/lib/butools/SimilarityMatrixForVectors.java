package jline.lib.butools;

import jline.util.matrix.Matrix;

import java.util.Comparator;
import java.util.stream.IntStream;

public final class SimilarityMatrixForVectors {
    private SimilarityMatrixForVectors() {}

    public static Matrix SimilarityMatrixForVectors(Matrix vecA, Matrix vecB) {
        int m = vecA.length();
        Matrix neg_vecA = vecA.copy();
        neg_vecA.scaleEq(-1.0);
        final double[] neg_vecA_array = neg_vecA.toArray1D();
        int[] ix = IntStream.range(0, neg_vecA_array.length).boxed()
                .sorted(Comparator.comparingDouble((Integer i) -> neg_vecA_array[i]))
                .mapToInt(Integer::intValue).toArray();
        Matrix P = new Matrix(m, m, m);
        for (int i = 0; i < m; i++) {
            P.set(i, ix[i], 1.0);
        }
        Matrix cp = P.mult(vecA);

        Matrix B = new Matrix(m, m, m * m);
        for (int i = 0; i < m; i++) {
            double cp_sum = 0.0;
            for (int j = 0; j <= i; j++) {
                cp_sum += cp.get(j, 0);
            }
            for (int j = 0; j <= i; j++) {
                B.set(i, j, vecB.get(i) / cp_sum);
            }
        }
        return B.mult(P);
    }
}
