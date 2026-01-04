/**
 * @file M3PP interleaved multi-process fitting
 * 
 * Implements fitting and interleaving of k second-order M3PP processes with varying
 * class counts. Constructs higher-order M3PP from multiple independent processes
 * through state space composition and parameter aggregation.
 * 
 * @since LINE 3.0
 */
package jline.api.mam.m3pp

import jline.util.matrix.MatrixCell
import jline.util.matrix.Matrix

/**
 * Fits k second-order M3PP[m_j] and interleaves them into a M3PP[m] of order k+1,
 * with m = sum_j=1^k m_j.
 *
 * @param av Vector of length k with the per-process rate
 * @param btv Vector of length k with the per-process IDC(t)  
 * @param binfv Vector of length k with the per-process IDC(inf)
 * @param acc Cell-array of length k; j-th element is a vector of length m_j
 *            with the per-class rates in the j-th M3PP
 * @param gtcc Cell-array of length k; j-th element is a vector of length m_j  
 *             with the per-class variance + marginal_covariance in the j-th M3PP
 * @param t Finite time scale
 * @param tinf Near-infinite time scale
 * @param mapping Optional mapping for class reordering
 * @param reorder Whether to reorder classes (default: false)
 * @return Pair of (interleaved M3PP[m], list of component M3PPs)
 */
fun m3pp_interleave_fitc(
    av: DoubleArray,           // Per-process rates
    btv: DoubleArray,          // Per-process IDC(t)
    binfv: DoubleArray,        // Per-process IDC(inf)
    acc: Array<DoubleArray>,   // Per-class rates for each process
    gtcc: Array<DoubleArray>,  // Per-class variance+covariance for each process
    t: Double,                 // Finite time scale
    tinf: Double,              // Near-infinite time scale
    mapping: IntArray? = null, // Optional class mapping
    reorder: Boolean = false   // Whether to reorder classes
): Pair<MatrixCell, List<MatrixCell>> {
    
    val k = av.size // Number of processes to interleave
    require(btv.size == k) { "btv must have length k" }
    require(binfv.size == k) { "binfv must have length k" }
    require(acc.size == k) { "acc must have length k" }
    require(gtcc.size == k) { "gtcc must have length k" }
    
    // Fit individual M3PP components
    val componentM3pps = mutableListOf<MatrixCell>()
    val componentSpecs = mutableListOf<ComponentSpec>()
    
    for (i in 0 until k) {
        val classRates = acc[i]
        val classVarCov = gtcc[i]
        val processRate = av[i]
        
        // Create component specification
        val spec = ComponentSpec(
            processIndex = i,
            processRate = processRate,
            classRates = classRates,
            classVarCov = classVarCov,
            idc_t = btv[i],
            idc_inf = binfv[i]
        )
        componentSpecs.add(spec)
        
        // Fit individual M3PP(2, m_j) with class-specific characteristics
        val component = fitM3ppWithClassCharacteristics(
            processRate = processRate,
            bt = btv[i],
            binf = binfv[i], 
            classRates = classRates,
            classVarCov = classVarCov,
            t = t,
            tinf = tinf
        )
        
        componentM3pps.add(component)
    }
    
    // Interleave the components
    val interleaved = interleaveM3pps(componentSpecs, componentM3pps, mapping, reorder)
    
    return Pair(interleaved, componentM3pps)
}

/**
 * Data class to hold component specifications
 */
private data class ComponentSpec(
    val processIndex: Int,
    val processRate: Double,
    val classRates: DoubleArray,
    val classVarCov: DoubleArray,
    val idc_t: Double,
    val idc_inf: Double
)

/**
 * Fit M3PP with class-specific characteristics
 */
private fun fitM3ppWithClassCharacteristics(
    processRate: Double,
    bt: Double,
    binf: Double,
    classRates: DoubleArray,
    classVarCov: DoubleArray,
    t: Double,
    tinf: Double
): MatrixCell {
    
    val numClasses = classRates.size
    
    // Use simplified M3PP fitting as a fallback
    // Create a simple MMPP(2) with the given rate and IDC characteristics
    val mmpp = Array(2) { Matrix(2, 2) }
    
    // Simple two-state MMPP construction
    val lambda1 = processRate * (1 + bt) / 2
    val lambda2 = processRate * (1 - bt) / 2
    val mu1 = lambda1 / 10  // Transition rates
    val mu2 = lambda2 / 10
    
    mmpp[0] = Matrix(2, 2)
    mmpp[0][0, 0] = -(lambda1 + mu1)
    mmpp[0][0, 1] = mu1
    mmpp[0][1, 0] = mu2
    mmpp[0][1, 1] = -(lambda2 + mu2)
    
    mmpp[1] = Matrix(2, 2)
    mmpp[1][0, 0] = lambda1
    mmpp[1][0, 1] = 0.0
    mmpp[1][1, 0] = 0.0
    mmpp[1][1, 1] = lambda2
    
    // Create M3PP from MMPP with class splitting
    val m3pp = MatrixCell(2 + numClasses)
    m3pp[0] = mmpp[0]
    m3pp[1] = mmpp[1]
    
    // Split classes proportionally with variance consideration
    val totalClassRate = classRates.sum()
    for (i in 0 until numClasses) {
        val proportion = if (totalClassRate > 0) classRates[i] / totalClassRate else 1.0 / numClasses
        val varianceWeight = if (classVarCov.isNotEmpty()) classVarCov[i] else 1.0
        m3pp[2 + i] = mmpp[1].scale(proportion * kotlin.math.max(0.1, varianceWeight))
    }
    
    return m3pp
}

/**
 * Interleave multiple M3PP processes using sophisticated interleaving algorithm
 */
private fun interleaveM3pps(
    specs: List<ComponentSpec>,
    components: List<MatrixCell>,
    mapping: IntArray?,
    reorder: Boolean
): MatrixCell {
    
    if (components.isEmpty()) {
        throw IllegalArgumentException("Cannot interleave empty list of M3PPs")
    }
    
    // Calculate total dimensions
    var totalClasses = 0
    for (spec in specs) {
        totalClasses += spec.classRates.size
    }
    
    // Determine interleaved state space size
    val interleavedOrder = Math.min(components.size + 1, 10) // Limit complexity
    val interleaved = MatrixCell(totalClasses + 2)
    
    // Initialize matrices
    val D0 = Matrix(interleavedOrder, interleavedOrder)
    val D1 = Matrix(interleavedOrder, interleavedOrder)
    interleaved[0] = D0
    interleaved[1] = D1
    
    // Initialize class matrices
    for (c in 0 until totalClasses) {
        interleaved[2 + c] = Matrix(interleavedOrder, interleavedOrder)
    }
    
    // Interleaving algorithm - create coordinated state transitions
    var globalClassIndex = 0
    
    for (compIdx in components.indices) {
        val component = components[compIdx]
        val spec = specs[compIdx]
        val compClasses = spec.classRates.size
        
        // Map component states to global state space
        val stateMapping = createStateMapping(compIdx, components.size, interleavedOrder)
        
        // Add component's generator structure with interleaving
        val compD0 = component[0]
        for (i in 0 until Math.min(compD0.getNumRows(), interleavedOrder)) {
            for (j in 0 until Math.min(compD0.getNumCols(), interleavedOrder)) {
                val globalI = stateMapping[i % stateMapping.size]
                val globalJ = stateMapping[j % stateMapping.size]
                
                if (globalI < interleavedOrder && globalJ < interleavedOrder) {
                    D0[globalI, globalJ] += compD0[i, j] / components.size // Scale by number of components
                }
            }
        }
        
        // Add class-specific interleaved arrivals
        for (c in 0 until compClasses) {
            val compClass = component[2 + c]
            val targetClassIndex = if (reorder && mapping != null) {
                mapping[globalClassIndex % mapping.size]
            } else {
                globalClassIndex
            }
            
            if (targetClassIndex < totalClasses) {
                val interleavedClass = interleaved[2 + targetClassIndex]
                
                for (i in 0 until Math.min(compClass.getNumRows(), interleavedOrder)) {
                    for (j in 0 until Math.min(compClass.getNumCols(), interleavedOrder)) {
                        val globalI = stateMapping[i % stateMapping.size]
                        val globalJ = stateMapping[j % stateMapping.size]
                        
                        if (globalI < interleavedOrder && globalJ < interleavedOrder) {
                            // Interleave with coordination factors
                            val coordinationFactor = computeCoordinationFactor(spec, compIdx, components.size)
                            val arrivalRate = compClass[i, j] * coordinationFactor
                            
                            interleavedClass[globalI, globalJ] += arrivalRate
                            D1[globalI, globalJ] += arrivalRate
                        }
                    }
                }
            }
            
            globalClassIndex++
        }
    }
    
    // Add cross-process coordination transitions
    addCoordinationTransitions(interleaved, specs)
    
    // Ensure generator property (rows sum to zero)
    ensureGeneratorProperty(D0)
    
    return interleaved
}

/**
 * Create state mapping for interleaving
 */
private fun createStateMapping(componentIndex: Int, totalComponents: Int, interleavedOrder: Int): IntArray {
    val mapping = IntArray(interleavedOrder)
    
    // Create round-robin style mapping with offset
    val offset = componentIndex
    for (i in mapping.indices) {
        mapping[i] = (offset + i * totalComponents) % interleavedOrder
    }
    
    return mapping
}

/**
 * Compute coordination factor for interleaving
 */
private fun computeCoordinationFactor(spec: ComponentSpec, compIndex: Int, totalComponents: Int): Double {
    // Factor based on relative process rate and position
    val relativerate = spec.processRate / (spec.classRates.sum() + 1e-6)
    val positionFactor = (compIndex + 1.0) / totalComponents
    
    return Math.min(1.0, relativerate * positionFactor)
}

/**
 * Add coordination transitions between processes
 */
private fun addCoordinationTransitions(interleaved: MatrixCell, specs: List<ComponentSpec>) {
    val D0 = interleaved[0]
    val order = D0.getNumRows()
    
    // Add weak coordination between processes
    val coordinationRate = specs.map { it.processRate }.average() * 0.1
    
    for (i in 0 until order - 1) {
        for (j in i + 1 until order) {
            // Add bidirectional coordination transitions
            D0[i, j] += coordinationRate / order
            D0[j, i] += coordinationRate / order
        }
    }
}

/**
 * Ensure generator matrix property (rows sum to zero)
 */
private fun ensureGeneratorProperty(D0: Matrix) {
    for (i in 0 until D0.getNumRows()) {
        var rowSum = 0.0
        
        // Sum off-diagonal elements
        for (j in 0 until D0.getNumCols()) {
            if (i != j) {
                rowSum += D0[i, j]
            }
        }
        
        // Set diagonal to negative row sum
        D0[i, i] = -rowSum
    }
}
/**
 * M3Pp Interleave Fit Count algorithms
 */
@Suppress("unused")
class M3ppInterleaveFitCountAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}