package jline.lib.perm

import jline.util.matrix.Matrix

/**
 * Implementation of the naive exact permanent computation.
 * 
 * This solver computes the exact permanent of a matrix by iterating over all
 * possible permutations and summing the products of matrix elements according
 * to each permutation. While this gives the exact result, it has O(n!) complexity
 * and is only practical for small matrices.
 * 
 * The permanent of an n×n matrix A is defined as:
 * perm(A) = Σ_π ∏_{i=1}^n a_{i,π(i)}
 * where the sum is over all permutations π of {1,2,...,n}.
 * 
 * @param matrix The matrix for which to compute the exact permanent
 * @param solve Whether to automatically run solve() after construction (default: false)
 */
class NaivePermanent(
    matrix: Matrix,
    solve: Boolean = false
) : PermSolver(matrix) {
    
    init {
        if (solve) {
            solve()
        }
    }
    
    override fun compute() {
        value = naive()
    }
    
    /**
     * Compute the exact permanent using the naive permutation-based method.
     * 
     * @return The exact permanent value
     */
    private fun naive(): Double {
        var permanent = 0.0
        
        // Generate all permutations of column indices
        val indices = IntArray(n) { it }
        val permutations = generatePermutations(indices)
        
        // Iterate over all permutations
        for (permutation in permutations) {
            var product = 1.0
            
            // Compute product for this permutation
            for (i in 0 until n) {
                product *= matrix.get(i, permutation[i])
            }
            
            permanent += product
        }
        
        return permanent
    }
    
    /**
     * Generate all permutations of the given array using Heap's algorithm.
     * 
     * @param array The array to permute
     * @return List of all permutations
     */
    private fun generatePermutations(array: IntArray): List<IntArray> {
        val result = mutableListOf<IntArray>()
        heapPermutation(array.copyOf(), array.size, result)
        return result
    }
    
    /**
     * Heap's algorithm for generating permutations.
     * 
     * @param array Current permutation being built
     * @param size Current size being permuted
     * @param result List to store all generated permutations
     */
    private fun heapPermutation(array: IntArray, size: Int, result: MutableList<IntArray>) {
        if (size == 1) {
            result.add(array.copyOf())
            return
        }
        
        for (i in 0 until size) {
            heapPermutation(array, size - 1, result)
            
            // If size is odd, swap first and last element
            if (size % 2 == 1) {
                val temp = array[0]
                array[0] = array[size - 1]
                array[size - 1] = temp
            } else {
                // If size is even, swap ith and last element
                val temp = array[i]
                array[i] = array[size - 1]
                array[size - 1] = temp
            }
        }
    }
}