package jline.api.moment;

import jline.util.matrix.Matrix;

/**
 * Factorial moments of a total count from the joint factorial moments of its
 * parts.
 *
 * <p>For N = N_1+...+N_d the Vandermonde convolution of falling factorials
 * gives
 *
 * <pre>
 *   f_n = sum_(|a|=n) (n! / prod_j a_j!) * F_a
 * </pre>
 *
 * <p>which holds for ANY joint law of the parts, marked or not, and is the
 * inverse direction of {@link Moment_joint_marking} whenever the marking is
 * multinomial. The order reached is limited by the smallest per-class order in
 * F, since the term a = n*e_j must be available for every j.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003.
 *
 * @since LINE 3.0
 */
public final class Moment_joint_aggregate {
    private Moment_joint_aggregate() {}

    /**
     * Factorial moments of the total count.
     *
     * @param F flattened joint factorial moments of the parts, row-major order
     * @param dims extents of the array, dims[j] = n_j+1
     * @return column vector of length min_j(n_j)+1 holding the factorial
     *         moments of the total
     */
    public static Matrix moment_joint_aggregate(double[] F, int[] dims) {
        int d = dims.length;
        int nel = 1;
        int nmax = Integer.MAX_VALUE;
        for (int j = 0; j < d; j++) {
            nel *= dims[j];
            nmax = Math.min(nmax, dims[j] - 1);
        }
        if (F.length != nel) {
            throw new IllegalArgumentException("moment_joint_aggregate: The array does not match "
                    + "dims.");
        }
        Matrix f = new Matrix(nmax + 1, 1);
        int[] a = new int[d];
        for (int ia = 0; ia < nel; ia++) {
            int n = 0;
            for (int j = 0; j < d; j++) {
                n += a[j];
            }
            if (n <= nmax) {
                double coef = 1.0;
                for (int t = 2; t <= n; t++) {
                    coef *= t;
                }
                for (int j = 0; j < d; j++) {
                    for (int t = 2; t <= a[j]; t++) {
                        coef /= t;
                    }
                }
                f.set(n, 0, f.get(n) + coef * F[ia]);
            }
            for (int l = d - 1; l >= 0; l--) {
                a[l]++;
                if (a[l] < dims[l]) {
                    break;
                }
                a[l] = 0;
            }
        }
        return f;
    }
}
