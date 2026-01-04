/**
 * @file M3PP(2,2) fitting using covariance approximation
 * 
 * Implements parameter fitting for second-order Marked Markov Modulated Poisson Process 
 * with two arrival classes, utilizing covariance approximation methods for accurate 
 * multi-class arrival modeling.
 * 
 * @since LINE 3.0
 */
package jline.api.mam.m3pp

import jline.util.matrix.Matrix
import jline.api.mam.Mmpp2_fitc_approx

/**
 * Fits a second-order Marked MMPP for two classes using covariance approximation.
 * 
 * @param a arrival rate
 * @param bt1 IDC at scale t1
 * @param bt2 IDC at scale t2
 * @param binf IDC for t->inf
 * @param m3t2 third central moment at scale t2
 * @param t1 first time scale
 * @param t2 second time scale
 * @param ai rates of the two classes
 * @param st3 count covariance between the two classes at scale t3
 * @param t3 third time scale
 * @return Array of matrices {D0, D1, D1_class1, D1_class2} representing the M3PP(2,2)
 */
class M3pp22_fitc_approx_cov {
    companion object {
        @JvmStatic
        fun m3pp22_fitc_approx_cov(
            a: Double,
            bt1: Double,
            bt2: Double,
            binf: Double,
            m3t2: Double,
            t1: Double,
            t2: Double,
            ai: DoubleArray,
            st3: Double,
            t3: Double
        ): Array<Matrix> {
            
            // Check consistency of per-class arrival rates
            if (kotlin.math.abs(a - ai.sum()) > 1e-8) {
                throw IllegalArgumentException("Inconsistent per-class arrival rates.")
            }
            
            // Fit underlying MMPP(2)
            val mmppFit = Mmpp2_fitc_approx.mmpp2_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2)
            
            // Fit M3PP(2,2)
            return M3pp22_fitc_approx_cov_multiclass.m3pp22_fitc_approx_cov_multiclass(
                mmppFit, ai, st3, t3
            )
        }
    }
}
/**
 * M3Pp22 Fit Count Approx Cov algorithms
 */
@Suppress("unused")
class M3pp22FitCountApproxCovAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}