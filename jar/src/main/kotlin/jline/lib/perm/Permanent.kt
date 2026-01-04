package jline.lib.perm

import jline.util.matrix.Matrix
import kotlin.math.pow

/**
 * Implementation of the MATLAB perm.m permanent computation algorithm.
 * 
 * This implementation mirrors the MATLAB perm.m function which computes the permanent
 * of a matrix by applying computational savings when some rows or columns are repeated.
 * It uses the inclusion-exclusion principle with multiplicities to efficiently handle
 * matrices with duplicate rows or columns.
 * 
 * The algorithm:
 * 1. Detects and counts repeated columns/rows
 * 2. Uses inclusion-exclusion with multinomial coefficients
 * 3. Applies pprod iterator for generating all combinations
 * 
 * @param matrix The matrix for which to compute the permanent
 * @param solve Whether to automatically run solve() after construction (default: false)
 */
class Permanent(
    matrix: Matrix,
    solve: Boolean = false
) : PermSolver(matrix) {
    
    init {
        if (solve) {
            solve()
        }
    }
    
    override fun compute() {
        value = computeWithMultiplicities()
    }
    
    /**
     * Computes the permanent using the MATLAB perm.m algorithm with multiplicities.
     * 
     * @return The permanent value
     */
    private fun computeWithMultiplicities(): Double {
        // Find unique columns and their multiplicities
        val (uniqueMatrix, multiplicities) = findUniqueColumnsWithMultiplicities()
        
        val R = multiplicities.size
        val n = multiplicities.sum()
        var value = 0.0
        
        // Initialize pprod iterator
        var f = IntArray(R) { 0 }
        
        do {
            var term = (-1.0).pow(f.sum())
            
            // Multinomial coefficients
            for (j in 0 until R) {
                term *= binomialCoefficient(multiplicities[j], f[j])
            }
            
            // Product term
            for (i in 0 until n) {
                var sumTerm = 0.0
                for (k in 0 until R) {
                    sumTerm += f[k] * uniqueMatrix.get(i, k)
                }
                term *= sumTerm
            }
            
            value += term
            f = pprodNext(f, multiplicities)
        } while (f.isNotEmpty() && !f.all { it == -1 })
        
        return (-1.0).pow(n) * value
    }
    
    /**
     * Finds unique columns in the matrix and their multiplicities.
     * If no repeated columns, checks for repeated rows instead.
     * 
     * @return Pair of (unique matrix, multiplicities array)
     */
    private fun findUniqueColumnsWithMultiplicities(): Pair<Matrix, IntArray> {
        // Find unique columns
        val columnMap = mutableMapOf<List<Double>, Int>()
        val uniqueColumns = mutableListOf<List<Double>>()
        val multiplicities = mutableListOf<Int>()
        
        for (j in 0 until matrix.getNumCols()) {
            val column = (0 until matrix.getNumRows()).map { i -> matrix.get(i, j) }
            val index = columnMap[column]
            if (index == null) {
                columnMap[column] = uniqueColumns.size
                uniqueColumns.add(column)
                multiplicities.add(1)
            } else {
                multiplicities[index]++
            }
        }
        
        // If no repeated columns (all multiplicities are 1), check for repeated rows
        if (multiplicities.all { it == 1 }) {
            val rowMap = mutableMapOf<List<Double>, Int>()
            val uniqueRows = mutableListOf<List<Double>>()
            val rowMultiplicities = mutableListOf<Int>()
            
            for (i in 0 until matrix.getNumRows()) {
                val row = (0 until matrix.getNumCols()).map { j -> matrix.get(i, j) }
                val index = rowMap[row]
                if (index == null) {
                    rowMap[row] = uniqueRows.size
                    uniqueRows.add(row)
                    rowMultiplicities.add(1)
                } else {
                    rowMultiplicities[index]++
                }
            }
            
            // If there are repeated rows, transpose and use row multiplicities
            if (rowMultiplicities.any { it > 1 }) {
                val transposedMatrix = Matrix(uniqueRows[0].size, uniqueRows.size)
                for (i in uniqueRows.indices) {
                    for (j in uniqueRows[i].indices) {
                        transposedMatrix.set(j, i, uniqueRows[i][j])
                    }
                }
                return Pair(transposedMatrix, rowMultiplicities.toIntArray())
            } else {
                // No repeated columns or rows, use original unique columns
                val uniqueMatrix = Matrix(matrix.getNumRows(), uniqueColumns.size)
                for (j in uniqueColumns.indices) {
                    for (i in uniqueColumns[j].indices) {
                        uniqueMatrix.set(i, j, uniqueColumns[j][i])
                    }
                }
                return Pair(uniqueMatrix, multiplicities.toIntArray())
            }
        } else {
            // There are repeated columns, use unique columns
            val uniqueMatrix = Matrix(matrix.getNumRows(), uniqueColumns.size)
            for (j in uniqueColumns.indices) {
                for (i in uniqueColumns[j].indices) {
                    uniqueMatrix.set(i, j, uniqueColumns[j][i])
                }
            }
            return Pair(uniqueMatrix, multiplicities.toIntArray())
        }
    }
    
    /**
     * Computes binomial coefficient C(n, k).
     * 
     * @param n The upper parameter
     * @param k The lower parameter
     * @return The binomial coefficient
     */
    private fun binomialCoefficient(n: Int, k: Int): Double {
        if (k > n || k < 0) return 0.0
        if (k == 0 || k == n) return 1.0
        
        var result = 1.0
        for (i in 1..minOf(k, n - k)) {
            result = result * (n - i + 1) / i
        }
        return result
    }
    
    /**
     * MATLAB pprod iterator - generates the next state in the sequence.
     * Returns a sequence of non-negative vectors less than a given vector.
     * 
     * @param current Current state vector
     * @param bounds Upper bounds vector
     * @return Next state vector, or array of -1s if sequence is complete
     */
    private fun pprodNext(current: IntArray, bounds: IntArray): IntArray {
        val n = current.copyOf()
        val N = bounds
        val R = N.size
        
        // Check if we've reached the maximum state
        if (n.zip(N).all { (a, b) -> a == b }) {
            return intArrayOf(-1)
        }
        
        var s = R - 1
        while (s >= 0 && n[s] == N[s]) {
            n[s] = 0
            s--
        }
        
        if (s < 0) {
            return intArrayOf(-1)
        }
        
        n[s]++
        return n
    }
}