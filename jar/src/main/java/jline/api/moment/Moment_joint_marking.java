package jline.api.moment;

import jline.util.matrix.Matrix;

/**
 * Joint factorial moments of the per-class counts under multinomial marking.
 *
 * <p>If a count N is marked independently, every event receiving class j with
 * probability p_j, then the per-class counts (N_1,...,N_d) have joint factorial
 * moments
 *
 * <pre>
 *   E[prod_j (N_j)_(a_j)] = (prod_j p_j^(a_j)) * f_(|a|)
 * </pre>
 *
 * <p>where f is the factorial moment sequence of the aggregate count N and
 * |a| = a_1+...+a_d. This is the counting-process counterpart of the marking
 * (class-splitting) formulas of the M3A fitters, and it is exact for the
 * per-class counts of a MAP marked in this i.i.d. way, in particular for an
 * MMAP whose marking probabilities do not depend on the phase.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003.
 *
 * @since LINE 3.0
 */
public final class Moment_joint_marking {
    private Moment_joint_marking() {}

    /**
     * Joint factorial moments of the per-class counts.
     *
     * @param f column vector of length n+1 holding the factorial moments of the
     *          aggregate count
     * @param p marking probabilities, one per class
     * @param dims maximum order per class; their sum must not exceed n
     * @return flattened joint factorial moments in row-major order, of extent
     *         (dims[0]+1)*...*(dims[d-1]+1)
     */
    public static double[] moment_joint_marking(Matrix f, double[] p, int[] dims) {
        int d = p.length;
        if (dims.length != d) {
            throw new IllegalArgumentException("moment_joint_marking: The arrays p and dims must "
                    + "have the same length.");
        }
        int tot = 0;
        for (int j = 0; j < d; j++) {
            if (dims[j] < 0) {
                throw new IllegalArgumentException("moment_joint_marking: The maximum orders must "
                        + "be nonnegative.");
            }
            tot += dims[j];
        }
        if (tot > f.length() - 1) {
            throw new IllegalArgumentException("moment_joint_marking: The aggregate factorial "
                    + "moments must reach order sum(dims) = " + tot + ".");
        }
        int nel = 1;
        for (int j = 0; j < d; j++) {
            nel *= dims[j] + 1;
        }
        double[] out = new double[nel];
        int[] a = new int[d];
        for (int ia = 0; ia < nel; ia++) {
            double coef = 1.0;
            int ord = 0;
            for (int j = 0; j < d; j++) {
                coef *= Math.pow(p[j], a[j]);
                ord += a[j];
            }
            out[ia] = coef * f.get(ord);
            for (int l = d - 1; l >= 0; l--) {
                a[l]++;
                if (a[l] <= dims[l]) {
                    break;
                }
                a[l] = 0;
            }
        }
        return out;
    }
}
