package jline.api.moment;

/**
 * Joint power (raw) moments from joint cumulants.
 *
 * <p>Runs the multivariate exponential-formula recursion of
 * {@link Moment_joint_cumulant_from_raw} forward,
 *
 * <pre>
 *   m_a = sum_(0&lt;b&lt;=a) prod_l nchoosek(a_l-[l=j], b_l-[l=j]) kappa_b m_(a-b)
 * </pre>
 *
 * <p>with m_0 = 1 and j the first dimension in which a is nonzero.
 *
 * <p>Reference:
 * V. P. Leonov and A. N. Shiryaev. On a method of calculation of
 * semi-invariants. Theory of Probability and its Applications,
 * 4(3):319-329, 1959.
 *
 * @since LINE 3.0
 */
public final class Moment_joint_raw_from_cumulant {
    private Moment_joint_raw_from_cumulant() {}

    /**
     * Converts joint cumulants into joint power moments.
     *
     * @param kappa flattened joint cumulants in row-major order; element 0 is
     *              ignored
     * @param dims extents of the array, dims[j] = n_j+1
     * @return flattened joint power moments, same layout, element 0 being 1
     */
    public static double[] moment_joint_raw_from_cumulant(double[] kappa, int[] dims) {
        return recur(kappa, dims, false);
    }

    /**
     * Shared sweep of the multivariate exponential formula. In the inverse
     * direction the input holds the moments and the output the cumulants; in
     * the forward direction the roles are exchanged. The two differ only in
     * which array is being filled in, so a single sweep serves both.
     *
     * @param in flattened input array in row-major order
     * @param dims extents of the array, dims[j] = n_j+1
     * @param inverse true to compute cumulants from moments, false for the
     *                forward direction
     * @return flattened output array, same layout
     */
    static double[] recur(double[] in, int[] dims, boolean inverse) {
        int d = dims.length;
        int nel = 1;
        for (int i = 0; i < d; i++) {
            nel *= dims[i];
        }
        if (in.length != nel) {
            throw new IllegalArgumentException("moment_joint_raw_from_cumulant: The array does "
                    + "not match dims.");
        }
        int[] stride = new int[d];
        stride[d - 1] = 1;
        for (int i = d - 2; i >= 0; i--) {
            stride[i] = stride[i + 1] * dims[i + 1];
        }
        double[] kappa = new double[nel];
        double[] mom = new double[nel];
        if (inverse) {
            System.arraycopy(in, 0, mom, 0, nel);
        } else {
            System.arraycopy(in, 0, kappa, 0, nel);
            mom[0] = 1.0;
        }
        int[] a = new int[d];
        int[] b = new int[d];
        for (int ia = 0; ia < nel; ia++) {
            int nz = -1;
            for (int l = 0; l < d; l++) {
                if (a[l] > 0) {
                    nz = l;
                    break;
                }
            }
            if (nz >= 0) {
                double acc = 0.0;
                int nb = 1;
                for (int l = 0; l < d; l++) {
                    nb *= a[l] + 1;
                    b[l] = 0;
                }
                for (int ib = 0; ib < nb; ib++) {
                    boolean zero = true;
                    boolean full = true;
                    for (int l = 0; l < d; l++) {
                        if (b[l] > 0) {
                            zero = false;
                        }
                        if (b[l] != a[l]) {
                            full = false;
                        }
                    }
                    if (!zero && b[nz] > 0 && !(inverse && full)) {
                        double c = 1.0;
                        for (int l = 0; l < d; l++) {
                            int shift = (l == nz) ? 1 : 0;
                            c *= Moment_binotrans.nchoosek(a[l] - shift, b[l] - shift);
                        }
                        if (c != 0.0) {
                            int ib2 = 0;
                            int irest = 0;
                            for (int l = 0; l < d; l++) {
                                ib2 += b[l] * stride[l];
                                irest += (a[l] - b[l]) * stride[l];
                            }
                            acc += c * kappa[ib2] * mom[irest];
                        }
                    }
                    for (int l = d - 1; l >= 0; l--) {
                        b[l]++;
                        if (b[l] <= a[l]) {
                            break;
                        }
                        b[l] = 0;
                    }
                }
                if (inverse) {
                    kappa[ia] = mom[ia] - acc;
                } else {
                    mom[ia] = acc;
                }
            }
            for (int l = d - 1; l >= 0; l--) {
                a[l]++;
                if (a[l] < dims[l]) {
                    break;
                }
                a[l] = 0;
            }
        }
        return inverse ? kappa : mom;
    }
}
