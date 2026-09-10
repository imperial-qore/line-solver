/**
 * @file MAP construction via Bernstein polynomial approximation
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.function.DoubleUnaryOperator;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_bernstein {
    private Map_bernstein() {}

    /**
     * Converts a distribution to a MAP via Bernstein polynomial approximation.
     */
    public static MatrixCell map_bernstein(DoubleUnaryOperator f, int n) {
        // Bernstein approximation normalizing constant
        double c = 0.0;
        for (int i = 1; i <= n; i++) {
            double xi = -Math.log((double) i / n);
            double fi = f.applyAsDouble(xi);
            if (Double.isFinite(fi) && fi > 0) {
                c += fi / i;
            }
        }

        if (c <= 0 || !Double.isFinite(c)) {
            return Map_erlang.map_erlang(1.0, n);
        }

        Matrix T = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            T.set(i, i, -(double) (i + 1));
            if (i < n - 1) {
                T.set(i, i + 1, (double) (i + 1));
            }
        }

        Matrix alpha = new Matrix(1, n);
        double alphaSum = 0.0;
        for (int i = 1; i <= n; i++) {
            double xi = -Math.log((double) i / n);
            double fi = f.applyAsDouble(xi);
            if (Double.isFinite(fi) && fi > 0) {
                double v = fi / (i * c);
                alpha.set(0, i - 1, v);
                alphaSum += v;
            }
        }

        if (alphaSum > 0) {
            for (int j = 0; j < n; j++) {
                alpha.set(0, j, alpha.get(0, j) / alphaSum);
            }
        } else {
            alpha.set(0, 0, 1.0);
        }

        Matrix P = Matrix.ones(n, 1).mult(alpha);
        Matrix negT = T.copy();
        negT.scaleEq(-1.0);
        Matrix D1 = negT.mult(P);

        MatrixCell MAP = new MatrixCell();
        MAP.set(0, T);
        MAP.set(1, D1);
        return MAP;
    }

    public static MatrixCell map_bernstein(DoubleUnaryOperator f) {
        return map_bernstein(f, 20);
    }
}
