package jline.lib.mom.solver

import jline.lib.mom.util.MomUtils
import jline.util.Maths
import jline.util.matrix.Matrix
import org.apache.commons.math3.linear.*
import kotlin.math.pow

/**
 * Linear solver implementation for MOM using double precision arithmetic
 * This is a port of the MATLAB linear.m solver
 */
class LinearSolver {
    
    /**
     * Solve queueing network using linear system approach with double precision
     * 
     * @param L Service rate matrix (M×R)
     * @param N Population vector (1×R)
     * @param Z Think time vector (1×R)
     * @return MomSolverResult containing performance measures
     */
    fun solve(L: RealMatrix, N: IntArray, Z: DoubleArray): MomSolverResult {
        val M = L.rowDimension
        val R = L.columnDimension
        
        // Initialize population tracker
        val n = IntArray(R) { 0 }
        
        // Process each class
        var g: RealVector? = null
        var gr: RealVector? = null
        
        for (r in 0 until R) {
            println("Linear solver: Setting up linear system for R*=${r + 1}")
            
            // Setup linear system
            val ls = SetupLinearSystem.setup(L, N, Z, r)
            
            if (r == 0) {
                // First class - g is all ones
                g = MatrixUtils.createRealVector(DoubleArray(M + 1) { 1.0 })
            } else {
                // Expand g for subsequent classes
                val G = expandG(g!!, gr!!, M, r)
                
                // Solve for new components
                val Gk = solveDense(ls.C, ls.Cg.operate(G).mapMultiply(-1.0))
                
                // Combine solutions
                g = combineVectors(Gk, G)
            }
            
            // Precompute matrix products for efficiency
            val CgDr = ls.Cg.multiply(ls.Dr)
            
            // Iterate through population levels
            for (nr in 0 until N[r]) {
                n[r]++
                
                // Compute G = Dr * g / n(r)
                val G = ls.Dr.operate(g!!).mapDivide(n[r].toDouble())
                
                // Compute right-hand side
                val rhs = ls.D.subtract(CgDr.scalarMultiply(1.0 / n[r])).operate(g)
                
                // Solve linear system
                val Gk = solveDense(ls.C, rhs)
                
                // Update g
                g = combineVectors(Gk, G)
            }
            
            // Save gr for next iteration
            gr = g
            n[r]++
            
            // Final G computation
            val G = ls.Dr.operate(g!!).mapDivide(n[r].toDouble())
            val Gk = solveDense(ls.C, ls.D.operate(g))
            g = combineVectors(Gk, G)
        }
        
        // Compute performance measures
        return computePerformanceMeasures(L, N, Z, g!!, n)
    }
    
    /**
     * Expand g vector when moving to next class
     */
    private fun expandG(g: RealVector, gr: RealVector, M: Int, r: Int): RealVector {
        val numNetworks = Maths.binomialCoeff(M + r - 2, r - 1).toInt()
        val newSize = numNetworks * r
        val G = MatrixUtils.createRealVector(DoubleArray(newSize))
        
        for (i in 0 until numNetworks) {
            // Copy previous r-1 components
            for (s in 0 until r - 1) {
                G.setEntry(i * r + s, g.getEntry(i * (r - 1) + s))
            }
            // Add component from gr
            G.setEntry(i * r + r - 1, gr.getEntry(i))
        }
        
        return G
    }
    
    /**
     * Solve dense linear system C * x = b
     */
    private fun solveDense(C: RealMatrix, b: RealVector): RealVector {
        return try {
            val solver = LUDecomposition(C).solver
            solver.solve(b)
        } catch (e: SingularMatrixException) {
            // Handle singular matrix - use pseudoinverse
            val svd = SingularValueDecomposition(C)
            svd.solver.solve(b)
        }
    }
    
    /**
     * Combine two vectors into one
     */
    private fun combineVectors(v1: RealVector, v2: RealVector): RealVector {
        val combined = MatrixUtils.createRealVector(DoubleArray(v1.dimension + v2.dimension))
        combined.setSubVector(0, v1)
        combined.setSubVector(v1.dimension, v2)
        return combined
    }
    
    /**
     * Compute performance measures from final g vector
     */
    private fun computePerformanceMeasures(
        L: RealMatrix,
        N: IntArray,
        Z: DoubleArray,
        g: RealVector,
        n: IntArray
    ): MomSolverResult {
        val M = L.rowDimension
        val R = L.columnDimension
        
        // Compute normalizing constant
        val G0 = g.getEntry(g.dimension - 1)
        
        // Initialize result matrices
        val X = MatrixUtils.createRealMatrix(M, R)
        val Q = MatrixUtils.createRealMatrix(M, R)
        
        // Compute throughput and queue lengths
        // This is simplified - full implementation needs proper marginal computation
        val totalPop = n.sum().toDouble()
        
        for (i in 0 until M) {
            for (r in 0 until R) {
                // Throughput calculation
                val lambda = N[r].toDouble() / (Z[r] + totalPop / L.getEntry(i, r))
                X.setEntry(i, r, lambda)
                
                // Queue length via Little's Law
                Q.setEntry(i, r, lambda / L.getEntry(i, r))
            }
        }
        
        // Normalize
        val norm = X.getFrobeniusNorm()
        if (norm > 0) {
            X.scalarMultiply(G0 / norm)
            Q.scalarMultiply(G0 / norm)
        }
        
        return MomSolverResult(X, Q, g.toArray())
    }
}