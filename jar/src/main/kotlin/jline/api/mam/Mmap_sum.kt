/**
 * @file Marked Markovian Arrival Process superposition operations
 * 
 * Computes superposition of MMAP processes creating independent multiclass arrival streams.
 * Essential for modeling combined traffic sources and multi-stream system analysis.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.MatrixCell
import jline.util.matrix.Matrix

/**
 * Computes the sum of two MMAPs, creating a superposition process.
 * This creates an independent superposition where arrivals from both processes
 * can occur simultaneously.
 *
 * @param mmap1 First MMAP
 * @param mmap2 Second MMAP
 * @return Superposed MMAP
 */
fun mmap_sum(mmap1: MatrixCell, mmap2: MatrixCell): MatrixCell {
    val order1 = mmap1[0].getNumRows()
    val order2 = mmap2[0].getNumRows()
    val classes1 = mmap1.size() - 2
    val classes2 = mmap2.size() - 2
    
    // Create compound state space (Kronecker product)
    val compoundOrder = order1 * order2
    val totalClasses = classes1 + classes2
    
    val sum = MatrixCell(totalClasses + 2)
    
    // Initialize matrices
    val D0 = Matrix(compoundOrder, compoundOrder)
    val D1 = Matrix(compoundOrder, compoundOrder)
    sum[0] = D0
    sum[1] = D1
    
    for (c in 0 until totalClasses) {
        sum[c + 2] = Matrix(compoundOrder, compoundOrder)
    }
    
    // Build compound generator matrix (Kronecker sum structure)
    val D0_1 = mmap1[0]
    val D0_2 = mmap2[0]
    
    for (i in 0 until order1) {
        for (j in 0 until order1) {
            for (k in 0 until order2) {
                for (l in 0 until order2) {
                    val compoundI = i * order2 + k
                    val compoundJ = j * order2 + l
                    
                    if (compoundI < compoundOrder && compoundJ < compoundOrder) {
                        if (i == j) {
                            // Second process transitions (first state unchanged)
                            D0[compoundI, compoundJ] += D0_2[k, l]
                        } else if (k == l) {
                            // First process transitions (second state unchanged)  
                            D0[compoundI, compoundJ] += D0_1[i, j]
                        }
                    }
                }
            }
        }
    }
    
    // Add arrivals from first MMAP
    var classOffset = 0
    for (c in 1 until mmap1.size()) {
        val D1_c = mmap1[c]
        val sumD_c = sum[classOffset + 2]
        
        for (i in 0 until order1) {
            for (j in 0 until order1) {
                for (k in 0 until order2) {
                    val compoundI = i * order2 + k
                    val compoundJ = j * order2 + k // Same second state
                    
                    if (compoundI < compoundOrder && compoundJ < compoundOrder) {
                        val arrivalRate = D1_c[i, j]
                        sumD_c[compoundI, compoundJ] = arrivalRate
                        D1[compoundI, compoundJ] += arrivalRate
                    }
                }
            }
        }
        classOffset++
    }
    
    // Add arrivals from second MMAP
    for (c in 1 until mmap2.size()) {
        val D2_c = mmap2[c]
        val sumD_c = sum[classOffset + 2]
        
        for (i in 0 until order1) {
            for (k in 0 until order2) {
                for (l in 0 until order2) {
                    val compoundI = i * order2 + k
                    val compoundJ = i * order2 + l // Same first state
                    
                    if (compoundI < compoundOrder && compoundJ < compoundOrder) {
                        val arrivalRate = D2_c[k, l]
                        sumD_c[compoundI, compoundJ] = arrivalRate
                        D1[compoundI, compoundJ] += arrivalRate
                    }
                }
            }
        }
        classOffset++
    }
    
    return sum
}

/**
 * Computes the sum of multiple MMAPs.
 *
 * @param mmaps List of MMAPs to sum
 * @return Superposed MMAP
 */
fun mmap_sum_multiple(mmaps: List<MatrixCell>): MatrixCell {
    require(mmaps.isNotEmpty()) {
        "Must provide at least one MMAP"
    }
    
    if (mmaps.size == 1) {
        return mmaps[0]
    }
    
    var result = mmaps[0]
    for (i in 1 until mmaps.size) {
        result = mmap_sum(result, mmaps[i])
    }
    
    return result
}

/**
 * Weighted sum of MMAPs with scaling factors.
 *
 * @param mmaps List of MMAPs to sum
 * @param weights Scaling factors for each MMAP
 * @return Weighted superposed MMAP
 */
fun mmap_sum_weighted(mmaps: List<MatrixCell>, weights: DoubleArray): MatrixCell {
    require(mmaps.size == weights.size) {
        "Number of weights must match number of MMAPs"
    }
    
    // Scale each MMAP by its weight
    val scaledMmaps = mmaps.mapIndexed { index, mmap ->
        val scaleFactor = Matrix.singleton(1.0 / weights[index]) // Scale mean inter-arrival time
        mmap_scale(mmap, scaleFactor) ?: mmap
    }
    
    return mmap_sum_multiple(scaledMmaps)
}

/**
 * Selective sum - only sum specific classes from each MMAP.
 *
 * @param mmaps List of MMAPs
 * @param classSelections List of class indices to include from each MMAP
 * @return Selective sum MMAP
 */
fun mmap_sum_selective(
    mmaps: List<MatrixCell>,
    classSelections: List<IntArray>
): MatrixCell {
    require(mmaps.size == classSelections.size) {
        "Must provide class selections for each MMAP"
    }
    
    // Create filtered MMAPs with only selected classes
    val filteredMmaps = mmaps.mapIndexed { index, mmap ->
        val selectedClasses = classSelections[index]
        val filtered = MatrixCell(selectedClasses.size + 2)
        
        // Copy D0 and initialize D1
        filtered[0] = mmap[0]
        val D1 = Matrix(mmap[0].getNumRows(), mmap[0].getNumCols())
        filtered[1] = D1
        
        // Copy selected class matrices
        for (i in selectedClasses.indices) {
            val classIndex = selectedClasses[i]
            if (classIndex + 2 < mmap.size()) {
                val classMatrix = mmap[classIndex + 2]
                filtered[i + 2] = classMatrix
                
                // Add to total arrivals
                for (row in 0 until classMatrix.getNumRows()) {
                    for (col in 0 until classMatrix.getNumCols()) {
                        D1[row, col] += classMatrix[row, col]
                    }
                }
            }
        }
        
        filtered
    }
    
    return mmap_sum_multiple(filteredMmaps)
}
/**
 * MMAP sum algorithms
 */
@Suppress("unused")
class MmapSum {
    companion object {
        // Class documentation marker for Dokka
    }
}