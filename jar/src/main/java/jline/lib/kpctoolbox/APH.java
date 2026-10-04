package jline.lib.kpctoolbox;

import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Facade class for Acyclic Phase-Type (APH) functions.
 */
public final class APH {

    private APH() {}

    public static Pair<Matrix, Matrix> fromMoments(Matrix moments) {
        double[] m = moments.toArray1D();
        Pair<MatrixCell, ?> result = jline.lib.kpctoolbox.aph.APH.aph_fit(m[0], m[1], m[2]);
        MatrixCell MAP = result.getFirst();
        Matrix D0 = MAP.get(0);
        Matrix D1 = MAP.get(1);
        int n = D0.getNumRows();

        Matrix alpha = new Matrix(1, n);
        double sumD1 = 0.0;
        for (int j = 0; j < n; j++) {
            double colSum = 0.0;
            for (int i = 0; i < n; i++) {
                colSum += D1.get(i, j);
            }
            alpha.set(0, j, colSum);
            sumD1 += colSum;
        }
        if (sumD1 > 0) {
            for (int j = 0; j < n; j++) {
                alpha.set(0, j, alpha.get(0, j) / sumD1);
            }
        }

        return new Pair<Matrix, Matrix>(alpha, D0);
    }

    public static double bounds2ndOrder(double m1, double m2, double m3) {
        return 1.5 * m2 * m2 / m1;
    }

    public static double bounds3rdOrder(double m1, double m2, double m3) {
        double n2 = m2 / (m1 * m1);
        double n3_lb = (4.0 / 3.0) * n2 * n2 - (1.0 / 3.0) * n2;
        return n3_lb * m1 * m2;
    }
}
