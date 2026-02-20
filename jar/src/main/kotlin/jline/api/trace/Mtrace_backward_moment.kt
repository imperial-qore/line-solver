package jline.api.trace

/**
 * Computes backward moments of a multi-class trace.
 * Backward moments characterize the time until the previous arrival.
 *
 * @param interArrivalTimes Array of inter-arrival times
 * @param classLabels Array of class labels for each arrival
 * @param order Moment order (1 for mean, 2 for second moment, etc.)
 * @param classIndex Optional specific class index (-1 for all classes)
 * @return Backward moment(s) - single value if classIndex specified, array if all classes
 */
fun mtrace_backward_moment(
    interArrivalTimes: DoubleArray,
    classLabels: IntArray,
    order: Int,
    classIndex: Int = -1
): Any {
    require(interArrivalTimes.size == classLabels.size) {
        "Inter-arrival times and class labels must have same length"
    }
    require(order >= 1) {
        "Moment order must be >= 1"
    }
    
    val numClasses = (classLabels.maxOrNull() ?: 0) + 1
    
    if (classIndex >= 0) {
        // Compute for specific class
        return computeClassBackwardMoment(interArrivalTimes, classLabels, classIndex, order)
    } else {
        // Compute for all classes
        return DoubleArray(numClasses) { c ->
            computeClassBackwardMoment(interArrivalTimes, classLabels, c, order)
        }
    }
}

/**
 * Compute backward moment for a specific class
 */
private fun computeClassBackwardMoment(
    interArrivalTimes: DoubleArray,
    classLabels: IntArray,
    targetClass: Int,
    order: Int
): Double {
    
    val backwardTimes = mutableListOf<Double>()
    var cumulativeTime = 0.0
    var lastClassTime = mutableMapOf<Int, Double>()
    
    for (i in interArrivalTimes.indices) {
        cumulativeTime += interArrivalTimes[i]
        val currentClass = classLabels[i]
        
        if (currentClass == targetClass) {
            // Find backward time to previous arrival of same class
            val backwardTime = if (lastClassTime.containsKey(targetClass)) {
                cumulativeTime - lastClassTime[targetClass]!!
            } else {
                cumulativeTime // Time since start if no previous arrival
            }
            backwardTimes.add(backwardTime)
        }
        
        lastClassTime[currentClass] = cumulativeTime
    }
    
    if (backwardTimes.isEmpty()) {
        return 0.0
    }
    
    // Compute moment of specified order
    return when (order) {
        1 -> backwardTimes.average()
        else -> backwardTimes.map { Math.pow(it, order.toDouble()) }.average()
    }
}

/**
 * Computes conditional backward moments given the forward recurrence time.
 *
 * @param interArrivalTimes Array of inter-arrival times
 * @param classLabels Array of class labels
 * @param conditioningClass Class to condition on
 * @param order Moment order
 * @return Conditional backward moment
 */
fun mtrace_backward_moment_conditional(
    interArrivalTimes: DoubleArray,
    classLabels: IntArray,
    conditioningClass: Int,
    order: Int
): Double {
    
    val conditionalBackwardTimes = mutableListOf<Double>()
    var cumulativeTime = 0.0
    val classArrivalTimes = mutableMapOf<Int, MutableList<Double>>()
    
    // Collect all arrival times by class
    for (i in interArrivalTimes.indices) {
        cumulativeTime += interArrivalTimes[i]
        val currentClass = classLabels[i]
        
        classArrivalTimes.computeIfAbsent(currentClass) { mutableListOf() }
            .add(cumulativeTime)
    }
    
    val conditioningArrivals = classArrivalTimes[conditioningClass] ?: return 0.0
    
    // For each conditioning arrival, compute backward times to other classes
    for (conditioningTime in conditioningArrivals) {
        for ((classIdx, arrivalTimes) in classArrivalTimes) {
            if (classIdx != conditioningClass) {
                // Find most recent arrival of this class before conditioning time
                val recentArrival = arrivalTimes.filter { it < conditioningTime }.maxOrNull()
                if (recentArrival != null) {
                    conditionalBackwardTimes.add(conditioningTime - recentArrival)
                }
            }
        }
    }
    
    if (conditionalBackwardTimes.isEmpty()) {
        return 0.0
    }
    
    return when (order) {
        1 -> conditionalBackwardTimes.average()
        else -> conditionalBackwardTimes.map { Math.pow(it, order.toDouble()) }.average()
    }
}
/**
 * Mtrace Backward Moment algorithms
 */
@Suppress("unused")
class MtraceBackwardMomentAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}