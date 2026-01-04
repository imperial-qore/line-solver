/**
 * @file Markov Modulated Poisson Process counting-based fitting
 * 
 * Fits MMPP(2) models using count statistics and index of dispersion criteria.
 * Implements Heffes and Lucantoni methodology for MMPP parameter estimation.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import org.apache.commons.math3.analysis.UnivariateFunction
import org.apache.commons.math3.analysis.solvers.BisectionSolver
import kotlin.math.*

/**
 * Fits a MMPP(2) according to [Heffes and Lucantoni, 1986].
 * 
 * @param mu arrival rate
 * @param bt1 IDC at scale t1
 * @param bt2 IDC at scale t2
 * @param binf IDC for t->inf
 * @param m3t2 third central moment at scale t2
 * @param t1 first time scale
 * @param t2 second time scale
 * @return Array of matrices {D0, D1} representing the MMPP(2)
 */
class Mmpp2_fitc {
    companion object {
        @JvmStatic
        fun mmpp2_fitc(
            mu: Double,
            bt1: Double,
            bt2: Double,
            binf: Double,
            m3t2: Double,
            t1: Double,
            t2: Double
        ): Array<Matrix> {
            
            // Degenerate case: r2 = 0, fit a marked Poisson process
            if (abs(binf - 1) < 1e-8 && abs(binf - bt1) < 1e-8) {
                return arrayOf(
                    Matrix(doubleArrayOf(-mu)),
                    Matrix(doubleArrayOf(mu))
                )
            }
            
            // Check feasibility
            if (!(binf > bt1 && bt1 > 1)) {
                // Return Poisson process for infeasible case
                return arrayOf(
                    Matrix(doubleArrayOf(-mu)),
                    Matrix(doubleArrayOf(mu))
                )
            }
            
            // Calculate c = (binf-1)/(binf-bt1)
            val c = (binf - 1) / (binf - bt1)
            
            // Solve for w using Lambert W function approximation
            // We need to solve: z = w * exp(w) where z = -c * exp(-c)
            val z = -c * exp(-c)
            
            // Use numerical solver for Lambert W function
            val solver = BisectionSolver(1e-12)
            val function = object : UnivariateFunction {
                override fun value(w: Double): Double {
                    return z - w * exp(w)
                }
            }
            
            // Initial guess and bounds for w
            val w = try {
                solver.solve(10000, function, -1.0, 10.0, 1.0)
            } catch (e: Exception) {
                // Fallback to approximation if solver fails
                if (z > -1/E) {
                    // Principal branch approximation
                    -1 + sqrt(2 * (1 + E * z))
                } else {
                    1.0
                }
            }
            
            val x = (w + c) / t1
            
            // Calculate coefficients
            val k1 = mu.pow(3) * t2.pow(3)
            val k2 = 3 * mu.pow(2) * (binf - 1) * t2.pow(2)
            val k3 = 3 * mu * (binf - 1) / x * t2
            val k4 = 3 * mu / x.pow(2) * (binf - 1) * t2 * exp(-x * t2)
            val k5 = 6 * mu / x.pow(3) * (binf - 1) * (1 - exp(-x * t2))
            val g1t2 = m3t2 + 3 * mu * t2 * (mu * t2 - 1) * bt2 + mu * t2 * (mu * t2 - 1) * (mu * t2 - 2)
            val h = (g1t2 - k1 - k2 - k3 * (-mu) - k4 * mu * x) / ((k3 / x) + k4 - k5)
            
            val r1: Double
            val r2: Double
            val l1: Double
            val l2: Double
            
            if (abs(h) < 1e-4) {
                // h = 0 case
                r1 = x / 2
                r2 = x / 2
                l2 = mu - 0.5 * sqrt(2 * (binf - 1) * mu * x)
                l1 = mu + 0.5 * sqrt(2 * (binf - 1) * mu * x)
            } else {
                val y = (binf - 1) * mu * x.pow(3) / (2 * h.pow(2))
                var tempR1 = x / 2 * (1 + 1 / sqrt(4 * y + 1))
                var tempR2 = x - tempR1
                
                // Ensure r1 >= r2
                if (tempR1 < tempR2) {
                    val tmp = tempR1
                    tempR1 = tempR2
                    tempR2 = tmp
                }
                r1 = tempR1
                r2 = tempR2
                
                val w_val = h / (r1 - r2)
                val w_min = -mu / r1 * (r1 + r2)
                val w_max = mu / r2 * (r1 + r2)
                
                if (w_val < w_min || w_val > w_max) {
                    // Ignore third moment to achieve feasibility
                    val z_val = (binf - 1) * x.pow(3) * mu
                    val u = x * z_val / (2 * mu.pow(2) * x.pow(2) + z_val)
                    val newR1 = u + (x - u) / 2
                    val newR2 = x - newR1
                    val delta = sqrt(z_val / (2 * newR1 * newR2))
                    l2 = mu - newR2 / x * delta
                    l1 = l2 + delta
                } else {
                    l2 = mu - h / (r1 - r2) * (r2 / (r1 + r2))
                    l1 = h / (r1 - r2) + l2
                }
            }
            
            // Construct MMPP matrices
            val D0 = Matrix(2, 2)
            D0.set(0, 0, -(r1 + l1))
            D0.set(0, 1, r1)
            D0.set(1, 0, r2)
            D0.set(1, 1, -(r2 + l2))
            
            val D1 = Matrix(2, 2)
            D1.set(0, 0, l1)
            D1.set(0, 1, 0.0)
            D1.set(1, 0, 0.0)
            D1.set(1, 1, l2)
            
            return arrayOf(D0, D1)
        }
    }
}
/**
 * Mmpp2 Fitc algorithms
 */
@Suppress("unused")
class Mmpp2FitcAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}