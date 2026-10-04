package jline.api.moment;

import jline.util.matrix.Matrix;

/**
 * Inverse binomial transform of a sequence.
 *
 * <p>Computes the inverse binomial transform of the sequence y_0,y_1,...,y_n,
 *
 * <pre>
 *   x_n = sum_{k=0}^{n} nchoosek(n,k) * y_k
 * </pre>
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003, eq. (9).
 *
 * @since LINE 3.0
 */
public final class Moment_binotransinv {
    private Moment_binotransinv() {}

    /**
     * Inverse binomial transform of a sequence.
     *
     * @param y column vector of length n+1 holding y_0,...,y_n, i.e. element i
     *          is the element of order i
     * @return column vector of length n+1 holding x_0,...,x_n
     */
    public static Matrix moment_binotransinv(Matrix y) {
        int n = y.length() - 1;
        Matrix x = new Matrix(n + 1, 1);
        for (int i = 0; i <= n; i++) {
            double xi = 0.0;
            for (int k = 0; k <= i; k++) {
                xi += Moment_binotrans.nchoosek(i, k) * y.get(k);
            }
            x.set(i, 0, xi);
        }
        return x;
    }
}
