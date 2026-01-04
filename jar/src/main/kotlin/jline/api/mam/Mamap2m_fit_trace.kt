/**
 * @file Markovian Arrival MAP with Marked arrivals trace-based fitting
 * 
 * Fits MAMAP(2,m) processes from empirical trace data with inter-arrival times and class labels.
 * Essential for data-driven modeling of real multiclass arrival processes.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.MatrixCell
import jline.util.matrix.Matrix

/**
 * Fits a MAMAP(2,m) to trace data.
 *
 * @param interArrivalTimes Array of inter-arrival times
 * @param classLabels Array of class labels for each arrival
 * @return Fitted MAMAP(2,m) model
 */
fun mamap2m_fit_trace(interArrivalTimes: DoubleArray, classLabels: IntArray): MatrixCell {
    require(interArrivalTimes.size == classLabels.size) {
        "Inter-arrival times and class labels must have same length"
    }
    
    // Extract basic moments
    val M1 = interArrivalTimes.average()
    val M2 = interArrivalTimes.map { it * it }.average()
    val M3 = interArrivalTimes.map { it * it * it }.average()
    
    // Extract class information
    val numClasses = (classLabels.maxOrNull() ?: 0) + 1
    val classCounts = IntArray(numClasses)
    for (label in classLabels) {
        if (label >= 0 && label < numClasses) {
            classCounts[label]++
        }
    }
    
    val classProbs = classCounts.map { it.toDouble() / classLabels.size }.toDoubleArray()
    val backwardMoments = computeClassBackwardMoments(interArrivalTimes, classLabels, numClasses)
    
    // Create proper matrices for the function call
    val dummyMap = MatrixCell(2)
    dummyMap[0] = Matrix(2, 2)
    dummyMap[1] = Matrix(2, 2)
    
    val F = doubleArrayOf(M2) // Forward moments
    val B = backwardMoments   // Backward moments
    
    val result = mamap2m_fit_fb_multiclass(dummyMap, classProbs, F, B)
    return result.first // Return just the MatrixCell part of the Triple
}

/**
 * Compute backward moments for each class
 */
private fun computeClassBackwardMoments(
    interArrivalTimes: DoubleArray,
    classLabels: IntArray,
    numClasses: Int
): DoubleArray {
    
    val backwardMoments = DoubleArray(numClasses)
    val lastArrivalTime = mutableMapOf<Int, Double>()
    var currentTime = 0.0
    
    for (i in interArrivalTimes.indices) {
        currentTime += interArrivalTimes[i]
        val currentClass = classLabels[i]
        
        if (currentClass >= 0 && currentClass < numClasses) {
            if (lastArrivalTime.containsKey(currentClass)) {
                val backwardTime = currentTime - lastArrivalTime[currentClass]!!
                backwardMoments[currentClass] = backwardTime
            } else {
                backwardMoments[currentClass] = currentTime
            }
            lastArrivalTime[currentClass] = currentTime
        }
    }
    
    // Fill in missing values with global mean
    val globalMean = interArrivalTimes.average()
    for (i in backwardMoments.indices) {
        if (backwardMoments[i] == 0.0) {
            backwardMoments[i] = globalMean
        }
    }
    
    return backwardMoments
}
/**
 * MAMAP 2m fit trace algorithms
 */
@Suppress("unused")
class Mamap2mFitTraceAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}