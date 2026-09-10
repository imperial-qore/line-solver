package jline.api.moment;

import jline.util.matrix.Matrix;

/**
 * Binomial transform of a sequence.
 *
 * <p>Computes the binomial transform of the sequence x_0,x_1,...,x_n into
 * y_0,y_1,...,y_n,
 *
 * <pre>
 *   y_n = sum_{k=0}^{n} (-1)^(n-k) * nchoosek(n,k) * x_k
 * </pre>
 *
 * <p>Applied to a moment sequence m_i = E[X^i] it returns the moments of the
 * unit downshift, y_i = E[(X-1)^i]. It is not an involution: its inverse is
 * {@link Moment_binotransinv}, the unsigned transform.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003, eq. (8).
 *
 * @since LINE 3.0
 */
public final class Moment_binotrans {
    private Moment_binotrans() {}

    /**
     * Binomial transform of a sequence.
     *
     * @param x column vector of length n+1 holding x_0,...,x_n, i.e. element i
     *          is the element of order i
     * @return column vector of length n+1 holding y_0,...,y_n
     */
    public static Matrix moment_binotrans(Matrix x) {
        int n = x.length() - 1;
        Matrix y = new Matrix(n + 1, 1);
        for (int i = 0; i <= n; i++) {
            double yi = 0.0;
            for (int k = 0; k <= i; k++) {
                double sign = ((i - k) % 2 == 0) ? 1.0 : -1.0;
                yi += sign * nchoosek(i, k) * x.get(k);
            }
            y.set(i, 0, yi);
        }
        return y;
    }

    /**
     * Binomial coefficient computed with a stable multiplicative loop that
     * avoids the overflow of the explicit factorial form.
     *
     * @param n total items
     * @param k items chosen
     * @return the binomial coefficient nchoosek(n,k)
     */
    static double nchoosek(int n, int k) {
        if (k < 0 || k > n) {
            return 0.0;
        }
        int kk = Math.min(k, n - k);
        double c = 1.0;
        for (int i = 0; i < kk; i++) {
            c = c * (n - i) / (i + 1);
        }
        return Math.rint(c);
    }
}
