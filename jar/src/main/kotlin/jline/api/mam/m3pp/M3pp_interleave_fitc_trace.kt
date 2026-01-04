/**
 * @file M3PP interleaved fitting from empirical trace data
 * 
 * Implements M3PP interleaving and fitting directly from empirical trace data 
 * with multiple arrival classes. Extracts statistical characteristics from 
 * inter-arrival time data for data-driven parameter estimation.
 * 
 * @since LINE 3.0
 */
package jline.api.mam.m3pp

import jline.util.matrix.MatrixCell
import jline.util.matrix.Matrix
import kotlin.math.*

/**
 * Interleaves k M3PP to fit a multi-class trace with m classes.
 * 
 * @param interArrivalTimes Array of inter-arrival times
 * @param classLabels Array of class labels corresponding to each arrival
 * @param t Finite time scale (optional, computed from trace if not provided)
 * @param tinf Near-infinite time scale (optional, computed from trace if not provided)
 * @param mapping m x k binary matrix mapping the m classes to k m3pp (optional, computed if not provided)
 * @return Pair of (fitted MMAP, list of component M3PPs)
 */
fun m3pp_interleave_fitc_trace(
    interArrivalTimes: DoubleArray,
    classLabels: IntArray,
    t: Double? = null,
    tinf: Double? = null,
    mapping: Array<BooleanArray>? = null
): Pair<MatrixCell, List<MatrixCell>> {
    
    require(interArrivalTimes.size == classLabels.size) {
        "Inter-arrival times and class labels must have the same length"
    }
    
    // Number of classes
    val uniqueLabels = classLabels.distinct().sorted()
    val m = uniqueLabels.size
    
    // Compute time scales if not provided
    val meanInterArrival = interArrivalTimes.average()
    val totalTime = interArrivalTimes.sum()
    val computedT = t ?: (10.0 * meanInterArrival)
    val computedTinf = tinf ?: maxOf(10.0 * computedT, totalTime / 100.0)
    
    // Compute counting processes
    val countsT = computeMulticlassCountsAtScale(interArrivalTimes, classLabels, computedT)
    val countsTinf = computeMulticlassCountsAtScale(interArrivalTimes, classLabels, computedTinf)
    
    // Overall arrival rate
    val arrivalRate = 1.0 / meanInterArrival
    
    // Per-class rates
    val classRates = DoubleArray(m) { i ->
        val classCount = classLabels.count { it == uniqueLabels[i] }
        (classCount.toDouble() / classLabels.size) * arrivalRate
    }
    
    // Compute or create mapping
    val finalMapping = mapping ?: computeCorrelationBasedMapping(countsTinf, m)
    val k = finalMapping[0].size // Number of M3PP processes
    
    // Validate mapping
    validateMapping(finalMapping, m)
    
    println("Fitting $m classes with $k M3PP(2,m_j) processes")
    
    // Create filters for each M3PP
    val filters = Array(k) { j ->
        BooleanArray(m) { i -> finalMapping[i][j] }
    }
    
    // Total rates for each M3PP
    val processRates = DoubleArray(k) { j ->
        classRates.filterIndexed { i, _ -> filters[j][i] }.sum()
    }
    
    // Per-class rates within each M3PP
    val classRatesPerProcess = Array(k) { j ->
        classRates.filterIndexed { i, _ -> filters[j][i] }.toDoubleArray()
    }
    
    // Compute IDC for each M3PP
    val idcT = DoubleArray(k)
    val idcTinf = DoubleArray(k)
    
    for (j in 0 until k) {
        // Aggregate counts for j-th M3PP
        val aggregateCountsT = DoubleArray(countsT.size) { windowIdx ->
            uniqueLabels.indices
                .filter { filters[j][it] }
                .sumOf { classIdx -> countsT[windowIdx][classIdx] }
        }
        
        val aggregateCountsTinf = DoubleArray(countsTinf.size) { windowIdx ->
            uniqueLabels.indices
                .filter { filters[j][it] }
                .sumOf { classIdx -> countsTinf[windowIdx][classIdx] }
        }
        
        // Compute IDC
        val meanT = aggregateCountsT.average()
        val varT = aggregateCountsT.map { (it - meanT) * (it - meanT) }.average()
        idcT[j] = varT / (processRates[j] * computedT)
        
        val meanTinf = aggregateCountsTinf.average()
        val varTinf = aggregateCountsTinf.map { (it - meanTinf) * (it - meanTinf) }.average()
        idcTinf[j] = varTinf / (processRates[j] * computedTinf)
    }
    
    // Compute per-class variance plus marginal covariance for each M3PP
    val gtc = Array(k) { j ->
        val numClassesInProcess = filters[j].count { it }
        val result = DoubleArray(numClassesInProcess)
        
        var resultIdx = 0
        for (i in 0 until m) {
            if (filters[j][i]) {
                // Per-class variance
                val classCounts = DoubleArray(countsT.size) { windowIdx -> countsT[windowIdx][i] }
                val classMean = classCounts.average()
                val classVariance = classCounts.map { (it - classMean) * (it - classMean) }.average()
                
                // Marginal covariance with other classes in the same M3PP
                val otherClassesCounts = DoubleArray(countsT.size) { windowIdx ->
                    uniqueLabels.indices
                        .filter { classIdx -> filters[j][classIdx] && classIdx != i }
                        .sumOf { classIdx -> countsT[windowIdx][classIdx] }
                }
                
                val covariance = if (otherClassesCounts.any { it != 0.0 }) {
                    val otherMean = otherClassesCounts.average()
                    classCounts.indices.map { idx ->
                        (classCounts[idx] - classMean) * (otherClassesCounts[idx] - otherMean)
                    }.average()
                } else {
                    0.0
                }
                
                result[resultIdx] = classVariance + covariance
                resultIdx++
            }
        }
        
        result
    }
    
    // Call the existing interleave fitting function
    return m3pp_interleave_fitc(
        av = processRates,
        btv = idcT,
        binfv = idcTinf,
        acc = classRatesPerProcess,
        gtcc = gtc,
        t = computedT,
        tinf = computedTinf,
        mapping = null
    )
}

/**
 * Compute multiclass counting processes at a given time scale
 */
private fun computeMulticlassCountsAtScale(
    interArrivalTimes: DoubleArray, 
    classLabels: IntArray, 
    scale: Double
): Array<DoubleArray> {
    
    val uniqueLabels = classLabels.distinct().sorted()
    val m = uniqueLabels.size
    val totalTime = interArrivalTimes.sum()
    val numWindows = maxOf(1, (totalTime / scale).toInt())
    
    val counts = Array(numWindows) { DoubleArray(m) }
    
    var currentTime = 0.0
    var windowIndex = 0
    
    for (i in interArrivalTimes.indices) {
        currentTime += interArrivalTimes[i]
        
        val targetWindowEnd = (windowIndex + 1) * scale
        if (currentTime <= targetWindowEnd) {
            val classIdx = uniqueLabels.indexOf(classLabels[i])
            if (classIdx >= 0) {
                counts[windowIndex][classIdx] += 1.0
            }
        } else {
            // Move to next window
            windowIndex++
            if (windowIndex >= numWindows) break
            
            val classIdx = uniqueLabels.indexOf(classLabels[i])
            if (classIdx >= 0) {
                counts[windowIndex][classIdx] = 1.0
            }
        }
    }
    
    return counts
}

/**
 * Compute correlation-based mapping using covariance analysis
 */
private fun computeCorrelationBasedMapping(
    counts: Array<DoubleArray>,
    m: Int,
    threshold: Double = 0.75
): Array<BooleanArray> {
    
    // Compute covariance matrix
    val covariance = computeCovarianceMatrix(counts)
    val variance = DoubleArray(m) { i ->
        val classCounts = DoubleArray(counts.size) { j -> counts[j][i] }
        val mean = classCounts.average()
        classCounts.map { (it - mean) * (it - mean) }.average()
    }
    
    // Group classes based on correlation
    val pool = (0 until m).toMutableSet()
    val groups = mutableListOf<List<Int>>()
    
    while (pool.isNotEmpty()) {
        // Find class with highest variance among remaining
        val pivot = pool.maxByOrNull { variance[it] } ?: break
        val currentGroup = mutableListOf(pivot)
        pool.remove(pivot)
        
        // Find correlated classes
        val toRemove = mutableSetOf<Int>()
        for (h in pool) {
            val correlation = if (variance[pivot] > 0 && variance[h] > 0) {
                covariance[pivot][h] / sqrt(variance[pivot] * variance[h])
            } else {
                0.0
            }
            
            if (correlation >= threshold) {
                currentGroup.add(h)
                toRemove.add(h)
            }
        }
        
        pool.removeAll(toRemove)
        groups.add(currentGroup)
    }
    
    // Create mapping matrix
    val k = groups.size
    val mapping = Array(m) { BooleanArray(k) }
    
    for (j in groups.indices) {
        for (classIdx in groups[j]) {
            mapping[classIdx][j] = true
        }
    }
    
    return mapping
}

/**
 * Compute covariance matrix for count data
 */
private fun computeCovarianceMatrix(counts: Array<DoubleArray>): Array<DoubleArray> {
    val m = counts[0].size
    val means = DoubleArray(m) { i ->
        counts.map { it[i] }.average()
    }
    
    val covariance = Array(m) { DoubleArray(m) }
    
    for (i in 0 until m) {
        for (j in 0 until m) {
            covariance[i][j] = counts.indices.map { k ->
                (counts[k][i] - means[i]) * (counts[k][j] - means[j])
            }.average()
        }
    }
    
    return covariance
}

/**
 * Validate mapping matrix
 */
private fun validateMapping(mapping: Array<BooleanArray>, m: Int) {
    require(mapping.size == m) { "Number of classes does not match mapping" }
    
    // Check that each class is mapped to exactly one M3PP
    for (i in 0 until m) {
        val mappingCount = mapping[i].count { it }
        require(mappingCount == 1) { "Invalid mapping: class $i mapped to $mappingCount processes" }
    }
}
/**
 * M3Pp Interleave Fit Count Trace algorithms
 */
@Suppress("unused")
class M3ppInterleaveFitCountTraceAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}