/**
 * @file Marked Markovian Arrival Process class-wise steady-state analysis
 * 
 * Computes steady-state probability vectors for each marked class in MMAP processes.
 * Essential for equilibrium analysis and performance metric calculations.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the steady-state probability vector for each class in an MMAP.
 * 
 * @param mmap MMAP represented as MatrixCell {D0, D1, D2, ..., Dm}
 * @return Matrix containing class probabilities (1 x m matrix)
 */
class Mmap_pie {
    companion object {
        @JvmStatic
        fun mmap_pie(mmap: MatrixCell): Matrix {
            if (mmap.size() < 2) {
                throw IllegalArgumentException("MMAP must have at least D0 and D1")
            }
            
            val m = mmap.size() - 2  // Number of classes
            
            // Handle single class case
            if (m <= 0) {
                val result = Matrix(1, 1)
                result.set(0, 0, 1.0)
                return result
            }
            
            // Get stationary distribution of the underlying MAP
            val mapCell = MatrixCell(2)
            mapCell[0] = mmap[0]  // D0
            mapCell[1] = mmap[1]  // Aggregate D1
            val piMap = map_pie(mapCell)
            
            // Compute class probabilities
            val classProbabilities = Matrix(1, m)
            var totalRate = 0.0
            
            // For each class, compute its arrival rate
            val classRates = DoubleArray(m)
            for (i in 0 until m) {
                val Di = mmap[2 + i]
                var classRate = 0.0
                
                // Sum all entries of Di weighted by stationary distribution
                for (j in 0 until Di.getNumRows()) {
                    for (k in 0 until Di.getNumCols()) {
                        classRate += piMap.get(j) * Di.get(j, k)
                    }
                }
                
                classRates[i] = classRate
                totalRate += classRate
            }
            
            // Normalize to get probabilities
            if (totalRate > 0) {
                for (i in 0 until m) {
                    classProbabilities.set(0, i, classRates[i] / totalRate)
                }
            } else {
                // Equal probabilities if no arrivals
                for (i in 0 until m) {
                    classProbabilities.set(0, i, 1.0 / m)
                }
            }
            
            return classProbabilities
        }
        
        /**
         * Overloaded version for Array<Matrix> input
         */
        @JvmStatic
        fun mmap_pie(mmap: Array<Matrix>): Matrix {
            return mmap_pie(MatrixCell(mmap))
        }
    }
}
/**
 * MMAP pie algorithms
 */
@Suppress("unused")
class MmapPieAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}