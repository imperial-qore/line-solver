package jline.api.mam

import jline.util.matrix.MatrixCell
import jline.util.matrix.Matrix

/**
 * Creates a second-order MMAP mixture from a collection of MMAPs.
 * This function reduces higher-order MMAPs to order 2 and creates a mixture.
 *
 * @param mmaps List of MMAPs to mix
 * @param weights Mixing weights (must sum to 1)
 * @return Second-order MMAP mixture
 */
fun mmap_mixture_order2(mmaps: List<MatrixCell>, weights: DoubleArray): MatrixCell {
    require(mmaps.isNotEmpty()) {
        "Must provide at least one MMAP"
    }
    require(weights.size == mmaps.size) {
        "Number of weights must match number of MMAPs"
    }
    require(Math.abs(weights.sum() - 1.0) < 1e-10) {
        "Weights must sum to 1"
    }
    
    // First, reduce each MMAP to order 2 if necessary
    val order2Mmaps = mmaps.map { mmap ->
        if (mmap[0].getNumRows() <= 2) {
            mmap // Already order 2 or less
        } else {
            reduceToOrder2(mmap)
        }
    }
    
    // Create mixture using the regular mixture function
    val mapsMap = mutableMapOf<Int, MatrixCell>()
    order2Mmaps.forEachIndexed { index, mmap ->
        mapsMap[index] = mmap
    }
    val alpha = Matrix(1, weights.size)
    weights.forEachIndexed { index, weight ->
        alpha.set(0, index, weight)
    }
    return mmap_mixture(alpha, mapsMap)
}

/**
 * Reduces an MMAP to order 2 using moment matching.
 *
 * @param mmap Original MMAP of any order
 * @return MMAP of order 2 with matched characteristics
 */
private fun reduceToOrder2(mmap: MatrixCell): MatrixCell {
    val originalOrder = mmap[0].getNumRows()
    val numClasses = mmap.size() - 2
    
    if (originalOrder <= 2) {
        return mmap
    }
    
    // Extract key characteristics from original MMAP
    val M1 = map_moment(mmap, 1)
    val M2 = map_moment(mmap, 2)
    val M3 = map_moment(mmap, 3)
    val classProbs = mmap_pc(mmap)
    
    // Create order-2 approximation using AMAP2 fitting
    val baseAmap2 = amap2_fit_gamma(M1, M2, M3, 0.0) // Start with zero correlation
    val baseOrder2 = baseAmap2.first ?: throw RuntimeException("Failed to create order-2 approximation")
    
    // Convert to multi-class MMAP by distributing classes
    val order2Mmap = MatrixCell(numClasses + 2)
    order2Mmap[0] = baseOrder2[0] // D0
    order2Mmap[1] = Matrix(2, 2) // Total arrivals matrix
    
    // Distribute class arrivals
    for (c in 0 until numClasses) {
        val classMatrix = Matrix(2, 2)
        val classWeight = if (c < classProbs.length()) classProbs.get(0, c) else 1.0 / numClasses
        
        // Scale base arrival matrix by class probability
        val D1 = baseOrder2[1]
        for (i in 0 until 2) {
            for (j in 0 until 2) {
                val value = D1[i, j] * classWeight
                classMatrix[i, j] = value
                order2Mmap[1][i, j] += value // Add to total arrivals
            }
        }
        
        order2Mmap[c + 2] = classMatrix
    }
    
    return order2Mmap
}

/**
 * Creates a second-order MMAP mixture with automatic weight selection.
 * Weights are chosen to minimize the approximation error.
 *
 * @param mmaps List of MMAPs to mix
 * @param targetCharacteristics Target moments to match (M1, M2, M3)
 * @return Pair of (optimized mixture, optimal weights)
 */
fun mmap_mixture_order2_optimal(
    mmaps: List<MatrixCell>,
    targetCharacteristics: DoubleArray
): Pair<MatrixCell, DoubleArray> {
    
    require(mmaps.isNotEmpty()) {
        "Must provide at least one MMAP"
    }
    require(targetCharacteristics.size >= 3) {
        "Must provide at least 3 target characteristics (M1, M2, M3)"
    }
    
    val n = mmaps.size
    
    // Compute characteristics of each component MMAP
    val componentCharacteristics = mmaps.map { mmap ->
        doubleArrayOf(
            map_moment(mmap, 1),
            map_moment(mmap, 2), 
            map_moment(mmap, 3)
        )
    }
    
    // Find optimal weights using simple optimization
    var bestWeights = DoubleArray(n) { 1.0 / n } // Start with uniform weights
    var bestError = Double.MAX_VALUE
    
    // Grid search for optimal weights (simplified approach)
    val gridSize = 20
    generateWeightCombinations(n, gridSize) { weights ->
        // Compute mixture characteristics
        val mixtureChar = DoubleArray(3)
        for (i in 0 until n) {
            for (j in 0 until 3) {
                mixtureChar[j] += weights[i] * componentCharacteristics[i][j]
            }
        }
        
        // Compute error
        var error = 0.0
        for (j in 0 until 3) {
            val relativeError = Math.abs(mixtureChar[j] - targetCharacteristics[j]) / 
                               Math.max(targetCharacteristics[j], 1e-6)
            error += relativeError * relativeError
        }
        
        if (error < bestError) {
            bestError = error
            bestWeights = weights.clone()
        }
    }
    
    // Create mixture with optimal weights
    val mixture = mmap_mixture_order2(mmaps, bestWeights)
    
    return Pair(mixture, bestWeights)
}

/**
 * Generate weight combinations for optimization
 */
private fun generateWeightCombinations(n: Int, gridSize: Int, callback: (DoubleArray) -> Unit) {
    val weights = DoubleArray(n)
    
    fun generateRecursive(index: Int, remainingWeight: Double) {
        if (index == n - 1) {
            weights[index] = remainingWeight
            callback(weights)
            return
        }
        
        val step = remainingWeight / (gridSize - index)
        for (i in 0..gridSize - index) {
            weights[index] = i * step
            generateRecursive(index + 1, remainingWeight - weights[index])
        }
    }
    
    generateRecursive(0, 1.0)
}
/**
 * MMAP mixture order2 algorithms
 */
@Suppress("unused")
class MmapMixtureOrder2Algo {
    companion object {
        // Class documentation marker for Dokka
    }
}