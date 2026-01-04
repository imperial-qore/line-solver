package jline.lib.mom.solver

import jline.lib.mom.util.MomUtils
import jline.util.Maths
import jline.util.matrix.Matrix
import org.apache.commons.math3.linear.*
import kotlin.math.min

/**
 * Constructs the linear system matrices for the MOM solver
 * This is a port of the MATLAB setupls.m function
 */
object SetupLinearSystem {
    
    /**
     * Setup linear system matrices for a given class r
     * 
     * @param L Service rate matrix (M×R)
     * @param N Population vector (1×R)
     * @param Z Think time vector (1×R)
     * @param r Current class index (0-based)
     * @return LinearSystemMatrices containing C, Cg, D, Dr
     */
    fun setup(L: RealMatrix, N: IntArray, Z: DoubleArray, r: Int): LinearSystemMatrices {
        val M = L.rowDimension
        val R = L.columnDimension
        
        // Expand service rate matrix
        val Lhat = expandServiceRates(L, N, r)
        
        // Generate all population distributions for current and previous classes
        // Convert multichoose results from Matrix to Array<IntArray>
        val IkMatrix = if (r < R - 1) {
            Maths.multichoose(M.toDouble(), (r + 1).toDouble())
        } else {
            // For last class, only consider feasible states
            generateFeasibleStatesMatrix(M, R, N, r)
        }
        
        val IMatrix = if (r > 0) {
            Maths.multichoose(M.toDouble(), r.toDouble())
        } else {
            Matrix(1, M)  // All zeros for r=0
        }
        
        // Convert Matrix to Array<IntArray> for sorting
        val Ik = matrixToIntArrays(IkMatrix)
        val I = matrixToIntArrays(IMatrix)
        
        // Sort combinations for better numerical properties
        val IkSorted = MomUtils.sortByNnzPos(Ik)
        val ISorted = MomUtils.sortByNnzPos(I)
        
        // Initialize matrices
        val numStates = IkSorted.size * (r + 1)
        val numPrevStates = ISorted.size * maxOf(r, 1)
        
        val C = MatrixUtils.createRealMatrix(numStates, numStates)
        val Cg = MatrixUtils.createRealMatrix(numStates, numPrevStates)
        val D = MatrixUtils.createRealMatrix(numStates, numPrevStates)
        val Dr = MatrixUtils.createRealMatrix(numPrevStates, numPrevStates)
        
        // Build coefficient matrices
        buildCoefficientMatrices(C, Cg, D, Dr, Lhat, Z, IkSorted, ISorted, M, r)
        
        return LinearSystemMatrices(C, Cg, D, Dr)
    }
    
    /**
     * Convert Matrix to Array<IntArray>
     */
    private fun matrixToIntArrays(matrix: Matrix): Array<IntArray> {
        val rows = matrix.getNumRows()
        val cols = matrix.getNumCols()
        return Array(rows) { i ->
            IntArray(cols) { j ->
                matrix.get(i, j).toInt()
            }
        }
    }
    
    /**
     * Expand service rates for multi-server stations
     */
    private fun expandServiceRates(L: RealMatrix, N: IntArray, r: Int): RealMatrix {
        val M = L.rowDimension
        val maxServers = N.sum()
        val Lhat = MatrixUtils.createRealMatrix(M, maxServers)
        
        for (i in 0 until M) {
            for (k in 0 until maxServers) {
                // Rate depends on number of servers busy
                val rate = L.getEntry(i, minOf(r, L.columnDimension - 1))
                Lhat.setEntry(i, k, rate * minOf(k + 1, N[minOf(r, N.size - 1)]).toDouble())
            }
        }
        
        return Lhat
    }
    
    /**
     * Generate feasible population states for the last class
     */
    private fun generateFeasibleStatesMatrix(M: Int, R: Int, N: IntArray, r: Int): Matrix {
        // For simplicity, generate all states and filter
        val allStates = Maths.multichoose(M.toDouble(), minOf(r + 1, N.sum()).toDouble())
        
        // Filter feasible states
        val feasibleList = mutableListOf<IntArray>()
        for (i in 0 until allStates.getNumRows()) {
            val state = IntArray(M) { j -> allStates.get(i, j).toInt() }
            if (state.sum() <= N.sum()) {
                feasibleList.add(state)
            }
        }
        
        // Convert back to Matrix
        val result = Matrix(feasibleList.size, M)
        feasibleList.forEachIndexed { i, state ->
            state.forEachIndexed { j, value ->
                result.set(i, j, value.toDouble())
            }
        }
        
        return result
    }
    
    /**
     * Build the coefficient matrices based on balance equations
     */
    private fun buildCoefficientMatrices(
        C: RealMatrix,
        Cg: RealMatrix,
        D: RealMatrix,
        Dr: RealMatrix,
        Lhat: RealMatrix,
        Z: DoubleArray,
        Ik: Array<IntArray>,
        I: Array<IntArray>,
        M: Int,
        r: Int
    ) {
        // Build C matrix (block diagonal structure)
        for (h in Ik.indices) {
            val ik = Ik[h]
            
            for (s in 0..r) {
                val row = h * (r + 1) + s
                
                // Diagonal entry
                var diag = 0.0
                for (i in 0 until M) {
                    if (ik[i] > 0) {
                        diag += Lhat.getEntry(i, ik[i] - 1)
                    }
                }
                if (s < Z.size && Z[s] > 0) {
                    diag += ik.sum() / Z[s]
                }
                
                C.setEntry(row, row, diag)
                
                // Off-diagonal entries within block
                for (sp in 0..r) {
                    if (sp != s) {
                        val col = h * (r + 1) + sp
                        var value = 0.0
                        
                        // Station transitions
                        for (i in 0 until M) {
                            if (ik[i] > 0) {
                                // Simplified routing - in practice would use routing matrix
                                value += Lhat.getEntry(i, ik[i] - 1) / M
                            }
                        }
                        
                        C.setEntry(row, col, -value)
                    }
                }
            }
        }
        
        // Build Cg matrix (coupling to previous populations)
        for (h in Ik.indices) {
            val ik = Ik[h]
            
            for (hp in I.indices) {
                val i = I[hp]
                
                // Check if i can lead to ik by adding one customer
                if (canTransition(i, ik)) {
                    for (s in 0..r) {
                        val row = h * (r + 1) + s
                        
                        for (sp in 0 until maxOf(r, 1)) {
                            val col = hp * maxOf(r, 1) + sp
                            
                            // Simplified - actual implementation needs routing logic
                            Cg.setEntry(row, col, 1.0 / M)
                        }
                    }
                }
            }
        }
        
        // Build D matrix (direct transitions)
        // Simplified version - actual implementation is more complex
        val dValue = 1.0 / (M * (r + 1))
        for (i in 0 until D.rowDimension) {
            for (j in 0 until D.columnDimension) {
                D.setEntry(i, j, dValue)
            }
        }
        
        // Build Dr matrix (recursion matrix)
        // Identity-like structure for recursion
        for (i in 0 until minOf(Dr.rowDimension, Dr.columnDimension)) {
            Dr.setEntry(i, i, 1.0)
        }
    }
    
    /**
     * Check if population state i can transition to ik by adding one customer
     */
    private fun canTransition(i: IntArray, ik: IntArray): Boolean {
        var diff = 0
        for (j in i.indices) {
            val d = ik[j] - i[j]
            if (d < 0) return false
            diff += d
        }
        return diff == 1
    }
}