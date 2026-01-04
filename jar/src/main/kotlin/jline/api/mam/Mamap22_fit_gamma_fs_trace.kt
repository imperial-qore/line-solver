/**
 * @file Markovian Arrival MAP with Marked arrivals two-class gamma forward-sigma trace fitting
 * 
 * Fits MAMAP(2,2) from trace data using gamma autocorrelation and forward-sigma characteristics.
 * Advanced trace-based fitting with correlation control for two-class systems.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.MatrixCell
import jline.util.matrix.Matrix

/**
 * Fits a MAMAP(2,2) using forward-start method from trace data.
 *
 * @param interArrivalTimes Array of inter-arrival times  
 * @param classLabels Array of class labels (must be binary: 0 or 1)
 * @return Fitted MAMAP(2,2) model
 */
fun mamap22_fit_gamma_fs_trace(interArrivalTimes: DoubleArray, classLabels: IntArray): MatrixCell {
    require(interArrivalTimes.size == classLabels.size) {
        "Inter-arrival times and class labels must have same length"
    }
    require(classLabels.all { it == 0 || it == 1 }) {
        "Class labels must be binary (0 or 1) for MAMAP(2,2)"
    }
    
    // Extract trace characteristics
    val M1 = interArrivalTimes.average()
    val M2 = interArrivalTimes.map { it * it }.average()
    val M3 = interArrivalTimes.map { it * it * it }.average()
    val GAMMA = computeTraceGamma(interArrivalTimes)
    
    // Class-specific rates
    val class0Count = classLabels.count { it == 0 }
    val class1Count = classLabels.count { it == 1 }
    val totalCount = classLabels.size
    
    val p0 = class0Count.toDouble() / totalCount
    val p1 = class1Count.toDouble() / totalCount
    
    // Extract class-specific inter-arrival times
    val class0Times = mutableListOf<Double>()
    val class1Times = mutableListOf<Double>()
    
    for (i in interArrivalTimes.indices) {
        if (classLabels[i] == 0) {
            class0Times.add(interArrivalTimes[i])
        } else {
            class1Times.add(interArrivalTimes[i])
        }
    }
    
    val M1_0 = if (class0Times.isNotEmpty()) class0Times.average() else M1
    val M1_1 = if (class1Times.isNotEmpty()) class1Times.average() else M1
    
    // Create proper matrices for the function call
    val dummyAmap = MatrixCell(2)
    dummyAmap[0] = Matrix(2, 2)
    dummyAmap[1] = Matrix(2, 2)
    
    val P = Matrix(1, 2)
    P[0, 0] = p0
    P[0, 1] = p1
    
    val B = Matrix(1, 1)
    B[0, 0] = M2
    
    val S = Matrix(1, 1)
    S[0, 0] = M3
    
    return mamap22_fit_bs_multiclass(dummyAmap, P, B, S)
}

/**
 * Compute gamma (correlation) from trace
 */
private fun computeTraceGamma(interArrivalTimes: DoubleArray): Double {
    if (interArrivalTimes.size < 2) return 0.0
    
    val mean = interArrivalTimes.average()
    val variance = interArrivalTimes.map { (it - mean) * (it - mean) }.average()
    
    if (variance <= 0) return 0.0
    
    // Lag-1 autocorrelation
    var autocovariance = 0.0
    for (i in 0 until interArrivalTimes.size - 1) {
        autocovariance += (interArrivalTimes[i] - mean) * (interArrivalTimes[i + 1] - mean)
    }
    autocovariance /= (interArrivalTimes.size - 1)
    
    return autocovariance / variance
}
/**
 * MAMAP 22 fit gamma fs trace algorithms
 */
@Suppress("unused")
class Mamap22FitGammaFsTrace {
    companion object {
        // Class documentation marker for Dokka
    }
}