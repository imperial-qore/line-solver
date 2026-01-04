package jline.api.mam

import jline.util.matrix.Matrix
import jline.api.mam.map_scale
import org.apache.commons.math3.optim.InitialGuess
import org.apache.commons.math3.optim.MaxEval
import org.apache.commons.math3.optim.PointValuePair
import org.apache.commons.math3.optim.SimpleBounds
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer
import kotlin.math.*

/**
 * Fits a second-order Marked MMPP using optimization.
 * 
 * @param a arrival rate
 * @param bt1 IDC at scale t1
 * @param bt2 IDC at scale t2
 * @param binf IDC for t->inf
 * @param m3t2 third central moment at scale t2
 * @param t1 first time scale
 * @param t2 second time scale
 * @return Array of matrices {D0, D1} representing the MMPP(2)
 */
class Mmpp2_fit_count_approx {
    companion object {
        @JvmStatic
        fun mmpp2_fit_count_approx(
            a: Double,
            bt1: Double,
            bt2: Double,
            binf: Double,
            m3t2: Double,
            t1: Double,
            t2: Double
        ): Array<Matrix> {
            
            // Define the objective function
            val objectiveFunction = ObjectiveFunction { params ->
                val l1 = params[0]
                val l2 = params[1]
                val r1 = params[2]
                val r2 = params[3]
                
                // Compute characteristics
                val xa = (l1 * r2 + l2 * r1) / (r1 + r2)
                val factor = a / xa
                
                // Compute xbt1
                val xbt1 = (r1 * (2 * l1.pow(2) * r2.pow(2) * t1 * factor - 2 * l2.pow(2) * r2 - 2 * l1.pow(2) * r2 + 
                           2 * l2.pow(2) * r2.pow(2) * t1 * factor + 4 * l1 * l2 * r2 + 
                           2 * l1.pow(2) * r2 * exp(- r1 * t1 * factor - r2 * t1 * factor) + 
                           2 * l2.pow(2) * r2 * exp(- r1 * t1 * factor - r2 * t1 * factor) - 
                           4 * l1 * l2 * r2.pow(2) * t1 * factor - 
                           4 * l1 * l2 * r2 * exp(- r1 * t1 * factor - r2 * t1 * factor)) + 
                           r1.pow(2) * (2 * r2 * t1 * factor * l1.pow(2) - 4 * r2 * t1 * factor * l1 * l2 + 
                           2 * r2 * t1 * factor * l2.pow(2))) / 
                           (t1 * factor * (r1 + r2).pow(3) * (l1 * r2 + l2 * r1)) + 1
                
                // Compute xbt2 if t1 != t2
                val xbt2 = if (abs(t1 - t2) > 1e-10) {
                    (r1 * (2 * l1.pow(2) * r2.pow(2) * t2 * factor - 2 * l2.pow(2) * r2 - 2 * l1.pow(2) * r2 + 
                           2 * l2.pow(2) * r2.pow(2) * t2 * factor + 4 * l1 * l2 * r2 + 
                           2 * l1.pow(2) * r2 * exp(- r1 * t2 * factor - r2 * t2 * factor) + 
                           2 * l2.pow(2) * r2 * exp(- r1 * t2 * factor - r2 * t2 * factor) - 
                           4 * l1 * l2 * r2.pow(2) * t2 * factor - 
                           4 * l1 * l2 * r2 * exp(- r1 * t2 * factor - r2 * t2 * factor)) + 
                           r1.pow(2) * (2 * r2 * t2 * factor * l1.pow(2) - 4 * r2 * t2 * factor * l1 * l2 + 
                           2 * r2 * t2 * factor * l2.pow(2))) / 
                           (t2 * factor * (r1 + r2).pow(3) * (l1 * r2 + l2 * r1)) + 1
                } else {
                    xbt1
                }
                
                // Compute xbinf
                val xbinf = ((2 * r2 * l1.pow(2) - 4 * r2 * l1 * l2 + 2 * r2 * l2.pow(2)) * r1.pow(2) + 
                            (2 * l1.pow(2) * r2.pow(2) - 4 * l1 * l2 * r2.pow(2) + 2 * l2.pow(2) * r2.pow(2)) * r1) / 
                            ((r1 + r2).pow(3) * (l1 * r2 + l2 * r1)) + 1
                
                // Compute third moment
                val t = t2 * factor
                val d = r1 + r2
                val p = (l1 - l2) * (r1 - r2)
                val xg3t = xa.pow(3) * t.pow(3) + 
                          3 * xa.pow(2) * (xbinf - 1) * t.pow(2) + 
                          3 * xa * (xbinf - 1) / d * (p / d - xa) * t + 
                          3 * xa / d.pow(2) * (xbinf - 1) * (p + xa * d) * t * exp(-t * d) - 
                          6 * xa / d.pow(3) * (xbinf - 1) * p * (1 - exp(-t * d))
                val xm3t2 = xg3t - 3 * xa * t * (xa * t - 1) * xbt2 - xa * t * (xa * t - 1) * (xa * t - 2)
                
                // Compute objective
                var obj = 0.0
                obj += (xa / a - 1).pow(2)
                obj += (xbt1 / bt1 - 1).pow(2)
                if (abs(t1 - t2) > 1e-10) {
                    obj += (xbt2 / bt2 - 1).pow(2)
                }
                obj += (xbinf / binf - 1).pow(2)
                obj += (xm3t2 / m3t2 - 1).pow(2)
                
                obj
            }
            
            // Initial guess
            val initialGuess = doubleArrayOf(
                a * 0.75,    // l1
                a * 1.5,     // l2
                1.0 / 3.0,   // r1
                2.0 / 3.0    // r2
            )
            
            // Bounds
            val lowerBounds = doubleArrayOf(1e-6, 1e-6, 0.0, 0.0)
            val upperBounds = doubleArrayOf(Double.MAX_VALUE, Double.MAX_VALUE, Double.MAX_VALUE, Double.MAX_VALUE)
            
            // Create optimizer
            val optimizer: MultivariateOptimizer = BOBYQAOptimizer(9) // 2*n+1 interpolation points
            
            // Optimize
            val optimum: PointValuePair = optimizer.optimize(
                objectiveFunction,
                GoalType.MINIMIZE,
                InitialGuess(initialGuess),
                SimpleBounds(lowerBounds, upperBounds),
                MaxEval(10000)
            )
            
            val xopt = optimum.point
            val l1 = xopt[0]
            val l2 = xopt[1]
            val r1 = xopt[2]
            val r2 = xopt[3]
            
            // Assemble MMPP
            val D0 = Matrix(2, 2)
            D0.set(0, 0, -(l1 + r1))
            D0.set(0, 1, r1)
            D0.set(1, 0, r2)
            D0.set(1, 1, -(l2 + r2))
            
            val D1 = Matrix(2, 2)
            D1.set(0, 0, l1)
            D1.set(0, 1, 0.0)
            D1.set(1, 0, 0.0)
            D1.set(1, 1, l2)
            
            // Scale to match exact rate
            val fit = map_scale(D0, D1, 1.0 / a)
            return arrayOf(fit[0], fit[1])
        }
    }
}