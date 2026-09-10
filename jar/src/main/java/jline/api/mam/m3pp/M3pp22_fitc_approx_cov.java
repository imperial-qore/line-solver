/**
 * @file M3PP(2,2) fitting using covariance approximation
 *
 * Implements parameter fitting for second-order Marked Markov Modulated Poisson Process
 * with two arrival classes, utilizing covariance approximation methods for accurate
 * multi-class arrival modeling.
 *
 * @since LINE 3.0
 */
package jline.api.mam.m3pp;

import jline.api.mam.Mmpp2_fitc_approx;
import jline.util.matrix.Matrix;

public final class M3pp22_fitc_approx_cov {
    private M3pp22_fitc_approx_cov() {}

    /**
     * Fits a second-order Marked MMPP for two classes using covariance approximation.
     *
     * @param a arrival rate
     * @param bt1 IDC at scale t1
     * @param bt2 IDC at scale t2
     * @param binf IDC for t-&gt;inf
     * @param m3t2 third central moment at scale t2
     * @param t1 first time scale
     * @param t2 second time scale
     * @param ai rates of the two classes
     * @param st3 count covariance between the two classes at scale t3
     * @param t3 third time scale
     * @return Array of matrices {D0, D1, D1_class1, D1_class2} representing the M3PP(2,2)
     */
    public static Matrix[] m3pp22_fitc_approx_cov(
            double a,
            double bt1,
            double bt2,
            double binf,
            double m3t2,
            double t1,
            double t2,
            double[] ai,
            double st3,
            double t3) {

        // Check consistency of per-class arrival rates
        double sum = 0.0;
        for (int k = 0; k < ai.length; k++) {
            sum += ai[k];
        }
        if (Math.abs(a - sum) > 1e-8) {
            throw new IllegalArgumentException("Inconsistent per-class arrival rates.");
        }

        // Fit underlying MMPP(2)
        Matrix[] mmppFit = Mmpp2_fitc_approx.mmpp2_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2);

        // Fit M3PP(2,2)
        return M3pp22_fitc_approx_cov_multiclass.m3pp22_fitc_approx_cov_multiclass(
                mmppFit, ai, st3, t3);
    }
}
