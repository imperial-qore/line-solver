package jline.lib.perm

import jline.util.matrix.Matrix
import kotlin.math.pow

/**
 * Implementation of Ryzer's algorithm to calculate the permanent.
 * 
 * Ryzer's algorithm computes the permanent using the inclusion-exclusion principle.
 * This implementation provides two modes:
 * - Gray code version: Uses Gray code ordering for efficient bit manipulation (default)
 * - Naive version: Direct implementation using combinations
 * 
 * The Gray code version is generally more efficient as it incrementally updates
 * row sums rather than recomputing them for each subset.
 * 
 * @param matrix The matrix for which to compute the permanent
 * @param mode The algorithm mode: "graycode" (default) or "naive"
 * @param solve Whether to automatically run solve() after construction (default: false)
 */
class RyzerPermanent(
    matrix: Matrix,
    private val mode: String = "graycode",
    solve: Boolean = false
) : PermSolver(matrix) {
    
    init {
        if (solve) {
            solve()
        }
    }
    
    override fun compute() {
        value = when (mode) {
            "graycode" -> ryzerAlgorithmGraycode()
            else -> ryzerAlgorithmNaive()
        }
    }
    
    /**
     * Gray code version of Ryzer's algorithm.
     * 
     * This version uses Gray code ordering to efficiently iterate through
     * all subsets of columns, incrementally updating row sums.
     * 
     * @return The permanent value
     */
    private fun ryzerAlgorithmGraycode(): Double {
        var permanent = 0.0
        
        // Start with empty row sum
        val rowSum = DoubleArray(n) { 0.0 }
        
        // Generate the bit sequence to modify in Gray code order
        val bitToModifyList = bitToModify(n)
        
        // The bit word that corresponds to matrix columns
        val currentBit = BooleanArray(n) { false }
        
        for (bitIndex in bitToModifyList) {
            // Change bit
            currentBit[bitIndex] = !currentBit[bitIndex]
            
            // Modify row sum
            val multiplier = if (currentBit[bitIndex]) 1.0 else -1.0
            for (i in 0 until n) {
                rowSum[i] += multiplier * matrix.get(i, bitIndex)
            }
            
            // Obtain the sign in the sum
            val nbCol = currentBit.count { it }
            
            // Calculate product of row sums
            var product = 1.0
            for (i in 0 until n) {
                product *= rowSum[i]
            }
            
            // Add to the permanent
            permanent += (-1.0).pow(nbCol) * product
        }
        
        permanent *= (-1.0).pow(n)
        return permanent
    }
    
    /**
     * Naive version of Ryzer's algorithm.
     * 
     * This version directly implements the inclusion-exclusion principle
     * using combinations of columns.
     * 
     * @return The permanent value
     */
    private fun ryzerAlgorithmNaive(): Double {
        var permanent = 0.0
        
        // Iterate over subset sizes
        for (i in 0..n) {
            val combinations = generateCombinations(n, n - i)
            
            for (combination in combinations) {
                // Calculate row sums for this combination of columns
                val rowSums = DoubleArray(n) { row ->
                    combination.sumOf { col -> matrix.get(row, col) }
                }
                
                // Calculate product of row sums
                var product = 1.0
                for (sum in rowSums) {
                    product *= sum
                }
                
                permanent += (-1.0).pow(n - i) * product
            }
        }
        
        permanent *= (-1.0).pow(n)
        return permanent
    }
    
    /**
     * Create the order of bits to modify to satisfy Gray code ordering.
     * 
     * @param m Gray code size
     * @return List of bit indices to modify
     */
    private fun bitToModify(m: Int): List<Int> {
        return if (m == 1) {
            listOf(0)
        } else {
            bitToModify(m - 1) + listOf(m - 1) + bitToModify(m - 1)
        }
    }
    
    /**
     * Generate all combinations of k elements from n elements.
     * 
     * @param n Total number of elements
     * @param k Number of elements to choose
     * @return List of all combinations
     */
    private fun generateCombinations(n: Int, k: Int): List<List<Int>> {
        if (k == 0) return listOf(emptyList())
        if (k > n) return emptyList()
        
        val result = mutableListOf<List<Int>>()
        val current = mutableListOf<Int>()
        
        fun backtrack(start: Int) {
            if (current.size == k) {
                result.add(current.toList())
                return
            }
            
            for (i in start until n) {
                current.add(i)
                backtrack(i + 1)
                current.removeAt(current.size - 1)
            }
        }
        
        backtrack(0)
        return result
    }
}