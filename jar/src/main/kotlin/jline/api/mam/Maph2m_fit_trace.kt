/**
 * @file Multi-class Absorbing Phase-type distribution trace-based fitting
 * 
 * Fits MAPH(2,m) from empirical trace data for multiclass service time modeling.
 * Essential for data-driven phase-type distribution modeling from real measurements.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.MatrixCell

/**
 * Fits a multi-class MAPH(2,m) model to trace data with class labels.
 *
 * @param interArrivalTimes Array of inter-arrival times
 * @param classLabels Array of class labels corresponding to each arrival (0-indexed)
 * @return Fitted MAPH(2,m) model
 */
fun maph2m_fit_trace(interArrivalTimes: DoubleArray, classLabels: IntArray): MatrixCell {
    require(interArrivalTimes.size == classLabels.size) {
        "Inter-arrival times and class labels must have the same length"
    }
    
    // Determine number of classes
    val maxClass = classLabels.maxOrNull() ?: 0
    val m = maxClass + 1
    
    // Compute overall moments
    val M1 = interArrivalTimes.average()
    val M2 = interArrivalTimes.map { it * it }.average()
    val M3 = interArrivalTimes.map { it * it * it }.average()
    
    // Compute class probabilities
    val classCounts = IntArray(m)
    for (label in classLabels) {
        if (label >= 0 && label < m) {
            classCounts[label]++
        }
    }
    val classProbs = classCounts.map { it.toDouble() / classLabels.size }.toDoubleArray()
    
    // Compute class-specific backward moments (conditional means)
    val classBackwardMoments = DoubleArray(m)
    for (i in 0 until m) {
        val classTimesIndexed = interArrivalTimes.indices
            .filter { classLabels[it] == i }
            .map { interArrivalTimes[it] }
        
        classBackwardMoments[i] = if (classTimesIndexed.isNotEmpty()) {
            classTimesIndexed.average()
        } else {
            M1 // Fallback to global mean if no samples for this class
        }
    }
    
    // Convert arrays to matrices
    val probMatrix = jline.util.matrix.Matrix(1, classProbs.size)
    for (i in classProbs.indices) {
        probMatrix[0, i] = classProbs[i]
    }
    
    val backMatrix = jline.util.matrix.Matrix(classBackwardMoments.size, 1)
    for (i in classBackwardMoments.indices) {
        backMatrix[i, 0] = classBackwardMoments[i]
    }
    
    return maph2m_fit(M1, M2, M3, probMatrix, backMatrix)
}

/**
 * Fits a multi-class MAPH(2,m) model to trace data with arrival time stamps and class labels.
 *
 * @param arrivalTimes Array of arrival time stamps (sorted)
 * @param classLabels Array of class labels corresponding to each arrival (0-indexed)
 * @return Fitted MAPH(2,m) model
 */
fun maph2m_fit_trace_timestamps(arrivalTimes: DoubleArray, classLabels: IntArray): MatrixCell {
    require(arrivalTimes.size == classLabels.size) {
        "Arrival times and class labels must have the same length"
    }
    require(arrivalTimes.size > 1) {
        "Need at least 2 arrivals to compute inter-arrival times"
    }
    
    // Convert timestamps to inter-arrival times
    val interArrivalTimes = DoubleArray(arrivalTimes.size - 1)
    for (i in 1 until arrivalTimes.size) {
        interArrivalTimes[i - 1] = arrivalTimes[i] - arrivalTimes[i - 1]
    }
    
    // Use the first n-1 class labels (since we have n-1 inter-arrival times)
    val adjustedClassLabels = classLabels.sliceArray(1 until classLabels.size)
    
    return maph2m_fit_trace(interArrivalTimes, adjustedClassLabels)
}
/**
 * MAPH 2m fit trace algorithms
 */
@Suppress("unused")
class Maph2mFitTraceAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}