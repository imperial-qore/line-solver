/**
 * @file M3PP(2,m) approximate count-based fitting with optimization
 * 
 * Implements approximation-based fitting for M3PP(2,m) using optimization methods. 
 * Employs quadratic programming and feasibility checks to fit multi-class arrival 
 * processes with improved computational efficiency over exact methods.
 * 
 * @since LINE 3.0
 */
package jline.api.mam.m3pp

import jline.GlobalConstants
import jline.VerboseLevel
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import jline.api.mam.Mmpp2_fitc_approx
import jline.api.mam.mmap_isfeasible
import kotlin.math.*
import com.quantego.josqp.Model
import com.quantego.josqp.OSQP

/**
 * Fits a second-order Marked MMPP using approximation with optimization.
 * 
 * @param a arrival rate
 * @param bt1 IDC at scale t1
 * @param bt2 IDC at scale t2
 * @param binf IDC for t->inf
 * @param m3t2 third central moment at scale t2
 * @param t1 first time scale
 * @param t2 second time scale
 * @param ai array where i-th element is the rate of class i
 * @param dvt3 array where i-th element is the delta of variance of class i and the variance
 *             of all other classes combined, at resolution t3
 * @param t3 third time scale
 * @return Array of matrices {D0, D1, D1_class1, D1_class2, ..., D1_classm} representing the M3PP(2,m)
 */
class M3pp2m_fitc_approx {
    companion object {
        @JvmStatic
        fun m3pp2m_fitc_approx(
            a: Double,
            bt1: Double,
            bt2: Double,
            binf: Double,
            m3t2: Double,
            t1: Double,
            t2: Double,
            ai: DoubleArray,
            dvt3: DoubleArray,
            t3: Double
        ): Array<Matrix> {
            
            // Check consistency of per-class arrival rates
            if (abs(a - ai.sum()) > 1e-8) {
                throw IllegalArgumentException("Inconsistent per-class arrival rates.")
            }
            
            val m = ai.size
            
            // Fit underlying MMPP(2)
            val mmppFit = Mmpp2_fitc_approx.mmpp2_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2)
            
            // Degenerate case: Poisson process
            if (mmppFit[0].getNumRows() == 1) {
                val D0 = mmppFit[0]
                val D1 = mmppFit[1]
                val fit = Array(2 + m) { Matrix(1, 1) }
                fit[0] = D0
                fit[1] = D1
                val pi = ai.map { it / a }
                for (i in 0 until m) {
                    fit[2 + i] = D1.scale(pi[i])
                }
                return fit
            }
            
            // Just a single class
            if (m == 1) {
                return arrayOf(mmppFit[0], mmppFit[1], mmppFit[1])
            }
            
            // Extract parameters
            val l1 = mmppFit[1].get(0, 0)
            val l2 = mmppFit[1].get(1, 1)
            val r1 = mmppFit[0].get(0, 1)
            val r2 = mmppFit[0].get(1, 0)
            val t = t3
            
            // Calculate coefficients
            val sinhTerm = sinh((r1 * t) / 2 + (r2 * t) / 2)
            val expTerm = exp(-(r1 * t) / 2 - (r2 * t) / 2)
            val sinhExp = sinhTerm * expTerm
            
            // Calculate q1i coefficients
            val q1i_ai = ((r1.pow(4) * t) / 2 + (r2.pow(4) * t) / 2 - l1 * r2.pow(3) * t + 
                l2 * r2.pow(3) * t + 2 * r1 * r2.pow(3) * t + 2 * r1.pow(3) * r2 * t + 
                3 * r1.pow(2) * r2.pow(2) * t + 2 * l1 * r2.pow(2) * sinhExp - 
                2 * l2 * r2.pow(2) * sinhExp - 2 * l1 * r1 * r2.pow(2) * t - 
                l1 * r1.pow(2) * r2 * t + 2 * l2 * r1 * r2.pow(2) * t + 
                l2 * r1.pow(2) * r2 * t + 2 * l1 * r1 * r2 * sinhExp - 
                2 * l2 * r1 * r2 * sinhExp) / 
                (l1 * r2 * (r1 + r2) * (2 * l1 * sinhExp - 2 * l2 * sinhExp - 
                l1 * r1 * t - l1 * r2 * t + l2 * r1 * t + l2 * r2 * t))
            
            val q1i_dvi = -(r1.pow(4) + 4 * r1.pow(3) * r2 + 6 * r1.pow(2) * r2.pow(2) + 
                4 * r1 * r2.pow(3) + r2.pow(4)) / 
                (4 * l1 * r2 * (r1 + r2) * (2 * l1 * sinhExp - 2 * l2 * sinhExp - 
                l1 * r1 * t - l1 * r2 * t + l2 * r1 * t + l2 * r2 * t))
            
            val q1i_const = -(l1 * r2.pow(4) * t + l2 * r1.pow(4) * t + 3 * l1 * r1 * r2.pow(3) * t + 
                l1 * r1.pow(3) * r2 * t + l2 * r1 * r2.pow(3) * t + 3 * l2 * r1.pow(3) * r2 * t + 
                3 * l1 * r1.pow(2) * r2.pow(2) * t + 2 * l1.pow(2) * r1 * r2.pow(2) * t + 
                2 * l1.pow(2) * r1.pow(2) * r2 * t + 3 * l2 * r1.pow(2) * r2.pow(2) * t + 
                2 * l2.pow(2) * r1 * r2.pow(2) * t + 2 * l2.pow(2) * r1.pow(2) * r2 * t - 
                4 * l1.pow(2) * r1 * r2 * sinhExp - 4 * l2.pow(2) * r1 * r2 * sinhExp - 
                4 * l1 * l2 * r1 * r2.pow(2) * t - 4 * l1 * l2 * r1.pow(2) * r2 * t + 
                8 * l1 * l2 * r1 * r2 * sinhExp) / 
                (4 * l1 * r2 * (r1 + r2) * (2 * l1 * sinhExp - 2 * l2 * sinhExp - 
                l1 * r1 * t - l1 * r2 * t + l2 * r1 * t + l2 * r2 * t))
            
            // Calculate q2i coefficients
            val denom2 = (r1 + r2) * (l2.pow(2) * r1.pow(2) * t - 2 * l2.pow(2) * r1 * sinhExp - 
                l1 * l2 * r1.pow(2) * t + l2.pow(2) * r1 * r2 * t + 
                2 * l1 * l2 * r1 * sinhExp - l1 * l2 * r1 * r2 * t)
            
            val q2i_ai = -((r1.pow(4) * t) / 2 + (r2.pow(4) * t) / 2 + l1 * r1.pow(3) * t - 
                l2 * r1.pow(3) * t + 2 * r1 * r2.pow(3) * t + 2 * r1.pow(3) * r2 * t + 
                3 * r1.pow(2) * r2.pow(2) * t - 2 * l1 * r1.pow(2) * sinhExp + 
                2 * l2 * r1.pow(2) * sinhExp + l1 * r1 * r2.pow(2) * t + 
                2 * l1 * r1.pow(2) * r2 * t - l2 * r1 * r2.pow(2) * t - 
                2 * l2 * r1.pow(2) * r2 * t - 2 * l1 * r1 * r2 * sinhExp + 
                2 * l2 * r1 * r2 * sinhExp) / denom2
            
            val q2i_dvi = (r1.pow(4) + 4 * r1.pow(3) * r2 + 6 * r1.pow(2) * r2.pow(2) + 
                4 * r1 * r2.pow(3) + r2.pow(4)) / (4 * denom2)
            
            val q2i_const = (l1 * r2.pow(4) * t + l2 * r1.pow(4) * t + 3 * l1 * r1 * r2.pow(3) * t + 
                l1 * r1.pow(3) * r2 * t + l2 * r1 * r2.pow(3) * t + 3 * l2 * r1.pow(3) * r2 * t + 
                3 * l1 * r1.pow(2) * r2.pow(2) * t + 2 * l1.pow(2) * r1 * r2.pow(2) * t + 
                2 * l1.pow(2) * r1.pow(2) * r2 * t + 3 * l2 * r1.pow(2) * r2.pow(2) * t + 
                2 * l2.pow(2) * r1 * r2.pow(2) * t + 2 * l2.pow(2) * r1.pow(2) * r2 * t - 
                4 * l1.pow(2) * r1 * r2 * sinhExp - 4 * l2.pow(2) * r1 * r2 * sinhExp - 
                4 * l1 * l2 * r1 * r2.pow(2) * t - 4 * l1 * l2 * r1.pow(2) * r2 * t + 
                8 * l1 * l2 * r1 * r2 * sinhExp) / (4 * denom2)
            
            // Set up quadratic programming problem
            // min sum_i (1/dvt3[i])^2 * (x[i] - dvt3[i])^2
            // s.t. constraints on q values
            
            val n = m  // Number of variables (delta variances)
            
            // Construct Hessian matrix H
            val H = Array(n) { i -> 
                DoubleArray(n) { j -> 
                    if (i == j) 2.0 / (dvt3[i] * dvt3[i]) else 0.0
                }
            }
            
            // Linear term h
            val h = DoubleArray(n) { i -> -2.0 / dvt3[i] }
            
            // Inequality constraints: ensure 0 <= q[j][i] <= 1 for all i,j
            // This gives us 2*m*2 = 4*m inequality constraints
            val constraintsList = mutableListOf<DoubleArray>()
            val boundsList = mutableListOf<Double>()
            
            for (i in 0 until m) {
                // q[0][i] >= 0: -q1i_ai * ai[i] - q1i_dvi * x[i] - q1i_const <= 0
                val c1 = DoubleArray(n) { 0.0 }
                c1[i] = -q1i_dvi
                constraintsList.add(c1)
                boundsList.add(q1i_ai * ai[i] + q1i_const)
                
                // q[0][i] <= 1: q1i_ai * ai[i] + q1i_dvi * x[i] + q1i_const <= 1
                val c2 = DoubleArray(n) { 0.0 }
                c2[i] = q1i_dvi
                constraintsList.add(c2)
                boundsList.add(1.0 - q1i_ai * ai[i] - q1i_const)
                
                // q[1][i] >= 0: -q2i_ai * ai[i] - q2i_dvi * x[i] - q2i_const <= 0
                val c3 = DoubleArray(n) { 0.0 }
                c3[i] = -q2i_dvi
                constraintsList.add(c3)
                boundsList.add(q2i_ai * ai[i] + q2i_const)
                
                // q[1][i] <= 1: q2i_ai * ai[i] + q2i_dvi * x[i] + q2i_const <= 1
                val c4 = DoubleArray(n) { 0.0 }
                c4[i] = q2i_dvi
                constraintsList.add(c4)
                boundsList.add(1.0 - q2i_ai * ai[i] - q2i_const)
            }
            
            // Equality constraints: sum of q values = 1 for each state
            val eqConstraintsList = mutableListOf<DoubleArray>()
            val eqBoundsList = mutableListOf<Double>()
            
            // sum_i q[0][i] = 1
            val eq1 = DoubleArray(n) { i -> q1i_dvi }
            eqConstraintsList.add(eq1)
            eqBoundsList.add(1.0 - m * q1i_const - q1i_ai * a)
            
            // sum_i q[1][i] = 1
            val eq2 = DoubleArray(n) { i -> q2i_dvi }
            eqConstraintsList.add(eq2)
            eqBoundsList.add(1.0 - m * q2i_const - q2i_ai * a)
            
            val A = constraintsList.toTypedArray()
            val b = boundsList.toDoubleArray()
            val Aeq = eqConstraintsList.toTypedArray()
            val beq = eqBoundsList.toDoubleArray()
            
            // Solve QP problem
            val optimalDv = solveQP(H, h, A, b, Aeq, beq, dvt3)
            
            // Calculate q values
            val q = Array(2) { DoubleArray(m) }
            for (i in 0 until m) {
                val Ai = ai[i]
                val Dvi = optimalDv[i]
                q[0][i] = q1i_ai * Ai + q1i_dvi * Dvi + q1i_const
                q[1][i] = q2i_ai * Ai + q2i_dvi * Dvi + q2i_const
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
            val D0 = mmppFit[0]
            val D1 = mmppFit[1]
            val fit = Array(2 + m) { Matrix(2, 2) }
            fit[0] = D0
            fit[1] = D1
            
            for (i in 0 until m) {
                val diag = Matrix(2, 2)
                diag.set(0, 0, q[0][i])
                diag.set(1, 1, q[1][i])
                fit[2 + i] = D1.elementMult(diag)
            }
            
            // Check feasibility
            if (!mmap_isfeasible(MatrixCell(fit))) {
                throw IllegalStateException("Infeasible fitted M3PP")
            }
            
            return fit
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
            
            // Add variables with lower bounds of 1e-6
            val vars = Array(n) { builder.addVariable().lb(1e-6) }
            
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
 * M3Pp2M Fit Count Approx algorithms
 */
@Suppress("unused")
class M3pp2mFitCountApproxAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}