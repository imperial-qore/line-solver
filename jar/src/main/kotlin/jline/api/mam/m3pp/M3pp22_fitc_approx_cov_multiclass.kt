/**
 * @file M3PP(2,2) multiclass fitting with covariance constraint optimization
 * 
 * Implements constrained optimization for fitting M3PP(2,2) parameters given an underlying 
 * MMPP(2), using covariance constraints between arrival classes to ensure feasibility 
 * and accurate multiclass statistical characterization.
 * 
 * @since LINE 3.0
 */
package jline.api.mam.m3pp

import jline.GlobalConstants.Inf
import jline.GlobalConstants.NegInf

import jline.util.matrix.Matrix
import kotlin.math.*

/**
 * Fits a M3PP(2,2) given the underlying MMPP(2).
 * 
 * @param mmpp underlying MMPP(2) as array {D0, D1}
 * @param ai rates of the two classes
 * @param st3 count covariance between the two classes at scale t3
 * @param t3 third time scale
 * @return Array of matrices {D0, D1, D1_class1, D1_class2} representing the M3PP(2,2)
 */
class M3pp22_fitc_approx_cov_multiclass {
    companion object {
        @JvmStatic
        fun m3pp22_fitc_approx_cov_multiclass(
            mmpp: Array<Matrix>,
            ai: DoubleArray,
            st3: Double,
            t3: Double
        ): Array<Matrix> {
            
            val m = ai.size
            
            if (m > 2) {
                throw IllegalArgumentException("No more than two classes supported")
            }
            
            var fit = mmpp
            
            // Degenerate case: Poisson process
            if (fit[0].getNumRows() == 1) {
                val D0 = fit[0]
                val D1 = fit[1]
                val result = Array(2 + m) { Matrix(1, 1) }
                result[0] = D0
                result[1] = D1
                val a = ai.sum()
                val pi = ai.map { it / a }
                for (i in 0 until m) {
                    result[2 + i] = D1.scale(pi[i])
                }
                return result
            }
            
            // Just a single class
            if (m == 1) {
                return arrayOf(fit[0], fit[1], fit[1])
            }
            
            // Extract parameters
            val l1 = fit[1].get(0, 0)
            val l2 = fit[1].get(1, 1)
            val r1 = fit[0].get(0, 1)
            val r2 = fit[0].get(1, 0)
            val t = t3
            val a1 = ai[0]
            
            // Calculate coefficients
            val w0 = (2 * r1 * (1 - exp(-(r1 + r2) * t) - (r1 + r2) * t) * 
                     (a1.pow(2) * (r1 + r2) - a1 * r2 * (l1 - l2))) / 
                     (r2 * (r1 + r2).pow(3))
            
            val w1 = -(2 * r1 * (1 - exp(-(r1 + r2) * t) - (r1 + r2) * t) * 
                      (2 * a1 * l2 * (r1 + r2) - l2 * r2 * (l1 - l2))) / 
                      (r2 * (r1 + r2).pow(3))
            
            val w2 = ((2 * l2.pow(2) * r2 * t) * (r1 + r2) + 
                     2 * l2.pow(2) * r1 * (1 - exp(-(r1 + r2) * t))) / 
                     (r2 * (r1 + r2).pow(2)) - (2 * l2.pow(2) * t) / r2
            
            val w3 = (r1 + r2) / (l2 * r1)
            val w4 = (l1 * r2) / (l2 * r1)
            
            // Bounds for the first and second root
            var L1 = NegInf
            var L2 = NegInf
            var U1 = Inf
            var U2 = Inf
            
            // If set to true, the first (second) root is never feasible
            var infeasible1 = false
            var infeasible2 = false
            
            // Defined for convenience
            val z = w0 - w1.pow(2) / (4 * w2)
            
            // Impose square root argument is >= 0
            if (w2 > 0) {
                L1 = max(L1, z)
                L2 = max(L2, z)
            } else if (w2 < 0) {
                U1 = min(U1, z)
                U2 = min(U2, z)
            }
            
            // Impose q2 >= 0
            // First root
            if (w1 >= 0) {
                L1 = max(L1, w0)
            } else if (w2 < 0) {
                infeasible1 = true
            }
            // Second root
            if (w1 <= 0) {
                U2 = min(U2, w0)
            } else if (w2 > 0) {
                infeasible2 = true
            }
            
            // Impose q1 >= 0
            var tmp = 2 * a1 * w3 * w2 + w1
            // First root
            if (tmp >= 0) {
                U1 = min(U1, z + tmp.pow(2) / (4 * w2))
            } else if (w2 > 0) {
                infeasible1 = true
            }
            // Second root
            if (tmp <= 0) {
                L2 = max(L2, z + tmp.pow(2) / (4 * w2))
            } else if (w2 < 0) {
                infeasible2 = true
            }
            
            // Impose q2 <= 1
            tmp = 2 * w2 + w1
            // First root
            if (tmp >= 0) {
                U1 = min(U1, z + tmp.pow(2) / (4 * w2))
            } else if (w2 > 0) {
                infeasible1 = true
            }
            // Second root
            if (tmp <= 0) {
                L2 = max(L2, z + tmp.pow(2) / (4 * w2))
            } else if (w2 < 0) {
                infeasible2 = true
            }
            
            // Impose q1 <= 1
            tmp = 2 * a1 * w2 * w3 - 2 * w2 * w4 + w1
            // First root
            if (tmp >= 0) {
                L1 = max(L1, z + tmp.pow(2) / (4 * w2))
            } else if (w2 < 0) {
                infeasible1 = true
            }
            // Second root
            if (tmp <= 0) {
                U2 = min(U2, z + tmp.pow(2) / (4 * w2))
            } else if (w2 > 0) {
                infeasible2 = true
            }
            
            if (infeasible1 && infeasible2) {
                throw IllegalStateException("Empty feasibility region. This should not happen.")
            }
            
            // Compute feasible covariance
            val sigma: Double
            val root: Int
            
            if (infeasible2) {
                sigma = max(min(st3, U1), L1)
                root = 1
            } else if (infeasible1) {
                sigma = max(min(st3, U2), L2)
                root = 2
            } else {
                val sigma1 = max(min(st3, U1), L1)
                val sigma2 = max(min(st3, U2), L2)
                if (abs(sigma1 - st3) < abs(sigma2 - st3)) {
                    sigma = sigma1
                    root = 1
                } else {
                    sigma = sigma2
                    root = 2
                }
            }
            
            // Compute parameters
            val q2 = if (root == 1) {
                (-w1 + sqrt(w1.pow(2) - 4 * w2 * (w0 - sigma))) / (2 * w2)
            } else {
                (-w1 - sqrt(w1.pow(2) - 4 * w2 * (w0 - sigma))) / (2 * w2)
            }
            
            val q1 = (a1 * (r1 + r2) - l2 * q2 * r1) / (l1 * r2)
            
            // Check feasibility
            val tol = 1e-8
            val q1Final: Double
            val q2Final: Double
            
            if (q1 >= -tol && q1 <= 1 + tol && q2 >= -tol && q2 <= 1 + tol) {
                q1Final = min(max(q1, 0.0), 1.0)
                q2Final = min(max(q2, 0.0), 1.0)
            } else {
                throw IllegalStateException("Parameters are infeasible. This should not happen.")
            }
            
            // Assemble M3PP[2]
            val D0 = fit[0]
            val D1 = fit[1]
            
            val D1_class1 = Matrix(2, 2)
            D1_class1.set(0, 0, D1.get(0, 0) * q1Final)
            D1_class1.set(1, 1, D1.get(1, 1) * q2Final)
            
            val D1_class2 = Matrix(2, 2)
            D1_class2.set(0, 0, D1.get(0, 0) * (1 - q1Final))
            D1_class2.set(1, 1, D1.get(1, 1) * (1 - q2Final))
            
            return arrayOf(D0, D1, D1_class1, D1_class2)
        }
    }
}
/**
 * M3Pp22 Fit Count Approx Cov Multiclass algorithms
 */
@Suppress("unused")
class M3pp22FitCountApproxCovMulticlassAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}