/**
 * @file M3PP(2,m) auto-gamma multiclass fitting with variance constraints
 * 
 * Implements multiclass auto-gamma fitting for M3PP(2,m) with variance and covariance 
 * constraints. Provides advanced statistical matching for complex multi-class arrival 
 * processes using feasibility checks and MMAP statistical functions.
 * 
 * @since LINE 3.0
 */
package jline.api.mam.m3pp

import jline.GlobalConstants
import jline.VerboseLevel
import jline.io.line_warning
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import jline.api.mam.map_count_mean
import jline.api.mam.mmap_isfeasible
import jline.api.mam.mmap_count_mean
import jline.api.mam.mmap_count_var
import jline.api.mam.mmap_count_mcov
import kotlin.math.*
import com.quantego.josqp.Model
import com.quantego.josqp.OSQP

/**
 * Fits a M3PP(2,m) given the underlying MMPP(2).
 * 
 * @param mmpp underlying MMPP(2) as array {D0, D1}
 * @param ai array where i-th element is the rate of class i
 * @param gt3 array where i-th element is the sum of the variance of class i and its
 *            covariance with all the other classes combined
 * @param t3 third time scale
 * @return Array of matrices {D0, D1, D1_class1, D1_class2, ..., D1_classm} representing the M3PP(2,m)
 */
class M3pp2m_fitc_approx_ag_multiclass {
    companion object {
        @JvmStatic
        fun m3pp2m_fitc_approx_ag_multiclass(
            mmpp: Array<Matrix>,
            ai: DoubleArray,
            gt3: DoubleArray,
            t3: Double
        ): Array<Matrix> {
            
            val m = ai.size
            var fit = mmpp
            
            // Total rate
            val a = map_count_mean(MatrixCell(mmpp[0], mmpp[1]), 1.0)
            
            // Check per-class rates consistency
            if (kotlin.math.abs(a - ai.sum()) > 1e-8) {
                throw IllegalArgumentException("Inconsistent per-class arrival rates.")
            }
            
            // Degenerate case: Poisson process
            if (fit[0].getNumRows() == 1) {
                val D0 = fit[0]
                val D1 = fit[1]
                val result = Array(2 + m) { Matrix(1, 1) }
                result[0] = D0
                result[1] = D1
                val pi = DoubleArray(ai.size) { i -> ai[i] / a }
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
            
            // Calculate f1 and f2
            val f1num = l1 * r2 * (2 * l2 * r1 - 2 * l1 * r1 + r1.pow(3) * t + r2.pow(3) * t + 
                2 * l1 * r1.pow(2) * t - 2 * l2 * r1.pow(2) * t + 3 * r1 * r2.pow(2) * t + 
                3 * r1.pow(2) * r2 * t + 2 * l1 * r1 * exp(-r1 * t - r2 * t) - 
                2 * l2 * r1 * exp(-r1 * t - r2 * t) + 2 * l1 * r1 * r2 * t - 2 * l2 * r1 * r2 * t)
            val f1 = f1num / (r1 + r2).pow(4)
            
            val f2num = l2 * r1 * (2 * l1 * r2 - 2 * l2 * r2 + r1.pow(3) * t + r2.pow(3) * t - 
                2 * l1 * r2.pow(2) * t + 2 * l2 * r2.pow(2) * t + 3 * r1 * r2.pow(2) * t + 
                3 * r1.pow(2) * r2 * t - 2 * l1 * r2 * exp(-r1 * t - r2 * t) + 
                2 * l2 * r2 * exp(-r1 * t - r2 * t) - 2 * l1 * r1 * r2 * t + 2 * l2 * r1 * r2 * t)
            val f2 = f2num / (r1 + r2).pow(4)
            
            val tmp = f1 * l2 * r1 - f2 * l1 * r2
            
            // Calculate coefficients
            val q1i_ai = -(f2 * (r1 + r2)) / tmp
            val q1i_gi = (l2 * r1) / tmp
            val q1i_const = 0.0
            
            val q2i_ai = (f1 * (r1 + r2)) / tmp
            val q2i_gi = -(l1 * r2) / tmp
            val q2i_const = 0.0
            
            // Set up quadratic programming problem
            // min sum_i (1/gt3[i])^2 * (x[i] - gt3[i])^2
            // s.t. constraints on q values
            
            val n = m  // Number of variables
            
            // Handle infinite or NaN values in gt3
            val cleanGt3 = DoubleArray(n) { i -> 
                if (!gt3[i].isFinite() || gt3[i] <= 0) 1.0 else gt3[i]
            }
            
            // Construct Hessian matrix H
            val H = Array(n) { i -> 
                DoubleArray(n) { j -> 
                    if (i == j) 2.0 / (cleanGt3[i] * cleanGt3[i]) else 0.0
                }
            }
            
            // Linear term h
            val h = DoubleArray(n) { i -> -2.0 / cleanGt3[i] }
            
            // Remove inf/nan values
            for (i in 0 until n) {
                if (!h[i].isFinite()) h[i] = 0.0
                for (j in 0 until n) {
                    if (!H[i][j].isFinite()) H[i][j] = if (i == j) 1.0 else 0.0
                }
            }
            
            // Inequality constraints: ensure 0 <= q[j][i] <= 1 for all i,j
            val constraintsList = mutableListOf<DoubleArray>()
            val boundsList = mutableListOf<Double>()
            
            for (i in 0 until m) {
                // q[0][i] >= 0: -q1i_ai * ai[i] - q1i_gi * x[i] - q1i_const <= 0
                val c1 = DoubleArray(n) { 0.0 }
                c1[i] = -q1i_gi
                constraintsList.add(c1)
                boundsList.add(q1i_ai * ai[i] + q1i_const)
                
                // q[0][i] <= 1: q1i_ai * ai[i] + q1i_gi * x[i] + q1i_const <= 1
                val c2 = DoubleArray(n) { 0.0 }
                c2[i] = q1i_gi
                constraintsList.add(c2)
                boundsList.add(1.0 - q1i_ai * ai[i] - q1i_const)
                
                // q[1][i] >= 0: -q2i_ai * ai[i] - q2i_gi * x[i] - q2i_const <= 0
                val c3 = DoubleArray(n) { 0.0 }
                c3[i] = -q2i_gi
                constraintsList.add(c3)
                boundsList.add(q2i_ai * ai[i] + q2i_const)
                
                // q[1][i] <= 1: q2i_ai * ai[i] + q2i_gi * x[i] + q2i_const <= 1
                val c4 = DoubleArray(n) { 0.0 }
                c4[i] = q2i_gi
                constraintsList.add(c4)
                boundsList.add(1.0 - q2i_ai * ai[i] - q2i_const)
            }
            
            // Equality constraints: sum of q values = 1 for each state
            val eqConstraintsList = mutableListOf<DoubleArray>()
            val eqBoundsList = mutableListOf<Double>()
            
            // sum_i q[0][i] = 1
            val eq1 = DoubleArray(n) { i -> q1i_gi }
            eqConstraintsList.add(eq1)
            eqBoundsList.add(1.0 - m * q1i_const - q1i_ai * a)
            
            // sum_i q[1][i] = 1
            val eq2 = DoubleArray(n) { i -> q2i_gi }
            eqConstraintsList.add(eq2)
            eqBoundsList.add(1.0 - m * q2i_const - q2i_ai * a)
            
            // Clean up constraint values
            for (i in boundsList.indices) {
                if (!boundsList[i].isFinite()) boundsList[i] = 0.0
            }
            for (i in eqBoundsList.indices) {
                if (!eqBoundsList[i].isFinite()) eqBoundsList[i] = 1.0
            }
            
            val A = constraintsList.toTypedArray()
            val b = boundsList.toDoubleArray()
            val Aeq = eqConstraintsList.toTypedArray()
            val beq = eqBoundsList.toDoubleArray()
            
            // Clean up constraint matrices
            for (i in A.indices) {
                for (j in A[i].indices) {
                    if (!A[i][j].isFinite()) A[i][j] = 0.0
                }
            }
            for (i in Aeq.indices) {
                for (j in Aeq[i].indices) {
                    if (!Aeq[i][j].isFinite()) Aeq[i][j] = 0.0
                }
            }
            
            // Solve QP problem
            val optimalG = solveQP(H, h, A, b, Aeq, beq, cleanGt3)
            
            // Calculate q values
            val q = Array(2) { DoubleArray(m) }
            for (i in 0 until m) {
                val Ai = ai[i]
                val Gi = optimalG[i]
                q[0][i] = q1i_ai * Ai + q1i_gi * Gi + q1i_const
                q[1][i] = q2i_ai * Ai + q2i_gi * Gi + q2i_const
            }
            
            // Ensure q values are valid probabilities
            for (j in 0..1) {
                // Ensure non-negative
                for (i in 0 until m) {
                    q[j][i] = max(0.0, q[j][i])
                }
                // Ensure sum <= 1
                val sum = q[j].sum()
                if (sum > 1.0) {
                    for (i in 0 until m) {
                        q[j][i] /= sum
                    }
                }
            }
            
            // Create result
            val D0 = fit[0]
            val D1 = fit[1]
            val result = Array(2 + m) { Matrix(2, 2) }
            result[0] = D0
            result[1] = D1
            
            for (i in 0 until m) {
                val diag = Matrix(2, 2)
                diag.set(0, 0, q[0][i])
                diag.set(1, 1, q[1][i])
                result[2 + i] = D1.elementMult(diag)
            }
            
            // Check feasibility
            if (!mmap_isfeasible(MatrixCell(result))) {
                // Warning: Infeasible fitted M3PP, but continue
                line_warning("m3pp2m_fitc_approx_ag_multiclass", "Infeasible fitted M3PP")
            }
            
            return result
        }
        
        /**
         * Solve quadratic programming problem using JOSQP
         */
        private fun solveQP(
            H: Array<DoubleArray>, h: DoubleArray,
            A: Array<DoubleArray>, b: DoubleArray,
            Aeq: Array<DoubleArray>, beq: DoubleArray,
            initialGuess: DoubleArray
        ): DoubleArray {
            
            val n = h.size
            
            // Create model builder
            val builder = Model.getBuilder()
            
            // Add variables with lower bounds of 0
            val vars = Array(n) { builder.addVariable().lb(0.0) }
            
            // Set quadratic objective: 0.5 * x' * H * x + h' * x
            val obj = builder.setObjective()
            
            // Add linear terms
            for (i in 0 until n) {
                obj.add(h[i], vars[i])
            }
            
            // Add quadratic terms
            for (i in 0 until n) {
                for (j in i until n) {
                    val coeff = if (i == j) H[i][j] else H[i][j] + H[j][i]
                    if (kotlin.math.abs(coeff) > 1e-12) {
                        obj.add(if (i == j) coeff / 2.0 else coeff / 2.0, vars[i], vars[j])
                    }
                }
            }
            
            obj.minimize()
            
            // Add inequality constraints A * x <= b
            for (i in A.indices) {
                val ctr = builder.addConstraint()
                for (j in 0 until n) {
                    if (kotlin.math.abs(A[i][j]) > 1e-12) {
                        ctr.add(A[i][j], vars[j])
                    }
                }
                ctr.leq(b[i])
            }
            
            // Add equality constraints Aeq * x = beq
            for (i in Aeq.indices) {
                val ctr = builder.addConstraint()
                for (j in 0 until n) {
                    if (kotlin.math.abs(Aeq[i][j]) > 1e-12) {
                        ctr.add(Aeq[i][j], vars[j])
                    }
                }
                ctr.eq(beq[i])
            }
            
            // Build and solve model
            val model = builder.build()
            model.param.verbose = GlobalConstants.getVerbose() != VerboseLevel.SILENT
            val status = model.solve()

            return if (status == OSQP.Status.SOLVED) {
                DoubleArray(n) { i -> model.getSolution(vars[i].index) }
            } else {
                // Return initial guess if OSQP fails
                initialGuess
            }
        }
    }
}
/**
 * M3Pp2M Fit Count Approx Ag Multiclass algorithms
 */
@Suppress("unused")
class M3pp2mFitCountApproxAgMulticlassAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}