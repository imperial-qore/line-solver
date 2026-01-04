package jline.lib.mom.util

/**
 * MOM-specific utility functions for combinatorial operations
 * Complements the existing jline.util.Maths class
 */
object MomUtils {
    
    /**
     * Sort combinations by number of non-zeros and their positions
     * This is used to optimize the order of processing in the solver
     * 
     * @param combinations Array of combinations to sort
     * @return Sorted array of combinations
     */
    fun sortByNnzPos(combinations: Array<IntArray>): Array<IntArray> {
        return combinations.sortedWith(compareBy(
            { row -> row.count { it != 0 } }, // First by number of non-zeros
            { row -> 
                // Then by position of non-zeros (lexicographic)
                val nonZeroPositions = row.indices.filter { row[it] != 0 }
                nonZeroPositions.joinToString("")
            }
        )).toTypedArray()
    }
    
    /**
     * Count non-zero elements in an array
     * 
     * @param array The array to count non-zeros in
     * @return Number of non-zero elements
     */
    fun countNonZeros(array: IntArray): Int {
        return array.count { it != 0 }
    }
    
    /**
     * Hash a population vector to a unique index
     * This is used for efficient lookup of population states
     * 
     * @param population The population vector
     * @param maxPop Maximum population per class (for hash calculation)
     * @return Hash index
     */
    fun hashPop(population: IntArray, maxPop: IntArray): Int {
        var hash = 0
        var multiplier = 1
        
        for (i in population.indices) {
            hash += population[i] * multiplier
            multiplier *= (maxPop[i] + 1)
        }
        
        return hash
    }
}