package jline.lib.butools;

import jline.util.matrix.Matrix;

public final class Poly {
    private Poly() {}

    public static Matrix poly(Matrix x) {
        int n = x.length();
        Matrix c = Matrix.concatColumns(Matrix.singleton(1.0), new Matrix(1, n, n), null);
        for (int j = 0; j < n; j++) {
            for (int i = j + 1; i >= 1; i--) {
                c.set(i, c.get(i) - x.get(j) * c.get(i - 1));
            }
        }
        return c;
    }
}
