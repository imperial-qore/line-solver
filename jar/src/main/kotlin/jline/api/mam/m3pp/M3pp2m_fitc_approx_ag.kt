/**
 * @file M3PP(2,m) auto-gamma approximation fitting
 * 
 * Implements auto-gamma approximation method for M3PP(2,m) parameter fitting.
 * Uses automatic gamma distribution matching for improved robustness in 
 * multi-class arrival process characterization with variance modeling.
 * 
 * @since LINE 3.0
 */
package jline.api.mam.m3pp

import jline.util.matrix.Matrix
import jline.api.mam.Mmpp2_fitc_approx

/**
 * Fits a second-order Marked MMPP using auto-gamma approach.
 * 
 * @param a arrival rate
 * @param bt1 IDC at scale t1
 * @param bt2 IDC at scale t2
 * @param binf IDC for t->inf
 * @param m3t2 third central moment at scale t2
 * @param t1 first time scale
 * @param t2 second time scale
 * @param ai array where i-th element is the rate of class i
 * @param gt3 array where i-th element is the sum of the variance of class i and its
 *            covariance with all the other classes combined, at scale t3
 * @param t3 third time scale
 * @return Array of matrices {D0, D1, D1_class1, D1_class2, ..., D1_classm} representing the M3PP(2,m)
 */
class M3pp2m_fitc_approx_ag {
    companion object {
        @JvmStatic
        fun m3pp2m_fitc_approx_ag(
            a: Double,
            bt1: Double,
            bt2: Double,
            binf: Double,
            m3t2: Double,
            t1: Double,
            t2: Double,
            ai: DoubleArray,
            gt3: DoubleArray,
            t3: Double
        ): Array<Matrix> {
            
            // Check consistency of per-class arrival rates
            if (kotlin.math.abs(a - ai.sum()) > 1e-8) {
                throw IllegalArgumentException("Inconsistent per-class arrival rates.")
            }
            
            // Fit underlying MMPP(2)
            val mmppFit = Mmpp2_fitc_approx.mmpp2_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2)
            
            // Fit M3PP(2,m)
            return M3pp2m_fitc_approx_ag_multiclass.m3pp2m_fitc_approx_ag_multiclass(
                mmppFit, ai, gt3, t3
            )
        }
    }
}
/**
 * M3Pp2M Fit Count Approx Ag algorithms
 */
@Suppress("unused")
class M3pp2mFitCountApproxAgAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}