package jline.api.mam

import jline.util.matrix.MatrixCell
import jline.util.matrix.Matrix

/**
 * Modulates an MMAP by another MMAP, creating a compound arrival process.
 * The modulating process controls the rate of the modulated process.
 *
 * @param baseMmap The base MMAP to be modulated
 * @param modulatingMmap The MMAP that provides modulation
 * @param modulationFactor Strength of modulation (default: 1.0)
 * @return Modulated MMAP
 */
fun mmap_modulate(
    baseMmap: MatrixCell, 
    modulatingMmap: MatrixCell,
    modulationFactor: Double = 1.0
): MatrixCell {
    
    val baseOrder = baseMmap[0].numRows
    val modulatingOrder = modulatingMmap[0].numRows
    val baseClasses = baseMmap.size() - 2
    val modulatingClasses = modulatingMmap.size() - 2
    
    // Create compound state space (Kronecker product structure)
    val compoundOrder = baseOrder * modulatingOrder
    val totalClasses = baseClasses // Keep base classes, modulation affects rates
    
    val modulated = MatrixCell(totalClasses + 2)
    
    // Initialize matrices
    val D0 = Matrix(compoundOrder, compoundOrder)
    val D1 = Matrix(compoundOrder, compoundOrder) 
    modulated[0] = D0
    modulated[1] = D1
    
    for (c in 0 until totalClasses) {
        modulated[c + 2] = Matrix(compoundOrder, compoundOrder)
    }
    
    // Extract matrices from input MMAPs
    val baseD0 = baseMmap[0]
    val modulatingD0 = modulatingMmap[0]
    
    // Compute modulating process steady-state probabilities
    val modulatingPi = computeSteadyState(modulatingMmap)
    
    // Compute modulation rates for each state of modulating process
    val modulationRates = DoubleArray(modulatingOrder)
    for (i in 0 until modulatingOrder) {
        var totalRate = 0.0
        for (c in 1 until modulatingMmap.size()) {
            val Dc = modulatingMmap[c]
            for (j in 0 until modulatingOrder) {
                totalRate += Dc[i, j]
            }
        }
        modulationRates[i] = Math.max(0.1, totalRate * modulationFactor)
    }
    
    // Build compound generator matrix
    for (i in 0 until baseOrder) {
        for (j in 0 until baseOrder) {
            for (k in 0 until modulatingOrder) {
                for (l in 0 until modulatingOrder) {
                    val compoundI = i * modulatingOrder + k
                    val compoundJ = j * modulatingOrder + l
                    
                    if (compoundI < compoundOrder && compoundJ < compoundOrder) {
                        if (i == j) {
                            // Modulating process transitions (base state unchanged)
                            D0[compoundI, compoundJ] += modulatingD0[k, l]
                        } else if (k == l) {
                            // Base process transitions (modulating state unchanged)
                            D0[compoundI, compoundJ] += baseD0[i, j] * modulationRates[k]
                        }
                    }
                }
            }
        }
    }
    
    // Build class arrival matrices with modulation
    for (c in 0 until totalClasses) {
        val baseDc = baseMmap[c + 2]
        val modulatedDc = modulated[c + 2]
        
        for (i in 0 until baseOrder) {
            for (j in 0 until baseOrder) {
                for (k in 0 until modulatingOrder) {
                    for (l in 0 until modulatingOrder) {
                        val compoundI = i * modulatingOrder + k
                        val compoundJ = j * modulatingOrder + l
                        
                        if (compoundI < compoundOrder && compoundJ < compoundOrder && k == l) {
                            // Base arrivals modulated by current modulating state
                            val modulatedRate = baseDc[i, j] * modulationRates[k]
                            modulatedDc[compoundI, compoundJ] = modulatedRate
                            D1[compoundI, compoundJ] += modulatedRate
                        }
                    }
                }
            }
        }
    }
    
    return modulated
}

/**
 * Time-varying modulation of an MMAP.
 * The modulation factor varies according to a specified pattern.
 *
 * @param baseMmap The base MMAP to modulate  
 * @param modulationPattern Function that maps time to modulation factor
 * @param timeHorizon Time horizon for discretization
 * @param numSteps Number of time steps
 * @return Time-modulated MMAP (approximated as mixture)
 */
fun mmap_modulate_time_varying(
    baseMmap: MatrixCell,
    modulationPattern: (Double) -> Double,
    timeHorizon: Double = 10.0,
    numSteps: Int = 10
): MatrixCell {
    
    // Create mixture of MMAPs for different time periods
    val timeStepMmaps = mutableListOf<MatrixCell>()
    val weights = DoubleArray(numSteps) { 1.0 / numSteps }
    
    for (step in 0 until numSteps) {
        val time = timeHorizon * step / numSteps
        val modulationFactor = modulationPattern(time)
        
        // Create modulated MMAP for this time period
        val scaleFactor = Matrix.singleton(1.0 / modulationFactor) // Scale mean inter-arrival time
        val modulated = mmap_scale(baseMmap, scaleFactor)
        if (modulated != null) {
            timeStepMmaps.add(modulated)
        }
    }
    
    // Return mixture representing time-varying modulation
    return mmap_mixture_order2(timeStepMmaps, weights)
}

/**
 * Cross-modulation between two MMAPs.  
 * Each MMAP modulates the other, creating mutual influence.
 *
 * @param mmap1 First MMAP
 * @param mmap2 Second MMAP
 * @param crossModulationStrength Strength of cross-modulation
 * @return Cross-modulated MMAP system
 */
fun mmap_cross_modulate(
    mmap1: MatrixCell,
    mmap2: MatrixCell, 
    crossModulationStrength: Double = 0.5
): MatrixCell {
    
    // First, mmap1 modulates mmap2
    val mmap2ModulatedBy1 = mmap_modulate(mmap2, mmap1, crossModulationStrength)
    
    // Then, mmap2 modulates mmap1  
    val mmap1ModulatedBy2 = mmap_modulate(mmap1, mmap2, crossModulationStrength)
    
    // Combine using superposition with equal weights
    return mmap_super(mmap1ModulatedBy2, mmap2ModulatedBy1) ?: throw RuntimeException("Failed to create superposition")
}

/**
 * Compute steady-state probabilities for an MMAP (simplified)
 */
private fun computeSteadyState(mmap: MatrixCell): DoubleArray {
    val order = mmap[0].getNumRows()
    val pi = DoubleArray(order) { 1.0 / order } // Uniform distribution as approximation
    
    // In practice, would solve pi * Q = 0 where Q is the generator matrix
    // This is a simplified version
    return pi
}
/**
 * MMAP modulate algorithms
 */
@Suppress("unused")
class MmapModulateAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}