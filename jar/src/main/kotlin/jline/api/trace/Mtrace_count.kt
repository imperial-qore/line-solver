package jline.api.trace

/**
 * Computes count statistics from a multi-class trace over specified time windows.
 *
 * @param interArrivalTimes Array of inter-arrival times
 * @param classLabels Array of class labels
 * @param windowSize Size of counting window
 * @param classIndex Optional specific class index (-1 for all classes)
 * @return Count statistics
 */
fun mtrace_count(
    interArrivalTimes: DoubleArray,
    classLabels: IntArray,
    windowSize: Double,
    classIndex: Int = -1
): CountStatistics {
    
    require(interArrivalTimes.size == classLabels.size) {
        "Inter-arrival times and class labels must have same length"
    }
    require(windowSize > 0) {
        "Window size must be positive"
    }
    
    val totalTime = interArrivalTimes.sum()
    val numWindows = Math.max(1, (totalTime / windowSize).toInt())
    val numClasses = (classLabels.maxOrNull() ?: 0) + 1
    
    // Generate count processes
    val counts = if (classIndex >= 0) {
        arrayOf(generateCountProcess(interArrivalTimes, classLabels, windowSize, classIndex))
    } else {
        Array(numClasses) { c ->
            generateCountProcess(interArrivalTimes, classLabels, windowSize, c)
        }
    }
    
    // Compute statistics for each count process
    val countStats = counts.map { countProcess ->
        computeCountProcessStatistics(countProcess, windowSize)
    }.toTypedArray()
    
    return CountStatistics(
        windowSize = windowSize,
        numWindows = numWindows,
        classStats = countStats,
        totalTime = totalTime
    )
}

/**
 * Data class for count statistics
 */
data class CountStatistics(
    val windowSize: Double,
    val numWindows: Int,
    val classStats: Array<SingleClassCountStats>,
    val totalTime: Double
)

data class SingleClassCountStats(
    val mean: Double,
    val variance: Double,
    val idc: Double, // Index of dispersion for counts
    val skewness: Double,
    val counts: IntArray
)

/**
 * Generate count process for a specific class
 */
private fun generateCountProcess(
    interArrivalTimes: DoubleArray,
    classLabels: IntArray,
    windowSize: Double,
    targetClass: Int
): IntArray {
    
    val totalTime = interArrivalTimes.sum()
    val numWindows = Math.max(1, (totalTime / windowSize).toInt())
    val counts = IntArray(numWindows)
    
    var currentTime = 0.0
    var windowIndex = 0
    
    for (i in interArrivalTimes.indices) {
        currentTime += interArrivalTimes[i]
        
        // Determine which window this arrival belongs to
        val targetWindow = Math.min((currentTime / windowSize).toInt(), numWindows - 1)
        
        // Count arrivals in appropriate windows
        if (classLabels[i] == targetClass) {
            // Handle case where arrival spans multiple windows
            for (w in windowIndex..targetWindow) {
                if (w < numWindows) {
                    counts[w]++
                }
            }
        }
        
        windowIndex = targetWindow
    }
    
    return counts
}

/**
 * Compute statistics for a count process
 */
private fun computeCountProcessStatistics(
    counts: IntArray,
    windowSize: Double
): SingleClassCountStats {
    
    if (counts.isEmpty()) {
        return SingleClassCountStats(0.0, 0.0, 0.0, 0.0, counts)
    }
    
    val mean = counts.average()
    val variance = if (counts.size > 1) {
        counts.map { (it - mean) * (it - mean) }.average()
    } else {
        0.0
    }
    
    val idc = if (mean > 0) variance / mean else 0.0
    
    val skewness = if (variance > 0 && counts.size > 2) {
        val thirdMoment = counts.map { Math.pow(it - mean, 3.0) }.average()
        thirdMoment / Math.pow(variance, 1.5)
    } else {
        0.0
    }
    
    return SingleClassCountStats(
        mean = mean,
        variance = variance, 
        idc = idc,
        skewness = skewness,
        counts = counts
    )
}

/**
 * Compute multi-scale count statistics
 */
fun mtrace_count_multiscale(
    interArrivalTimes: DoubleArray,
    classLabels: IntArray,
    windowSizes: DoubleArray
): Array<CountStatistics> {
    
    return windowSizes.map { windowSize ->
        mtrace_count(interArrivalTimes, classLabels, windowSize)
    }.toTypedArray()
}

/**
 * Convert inter-arrival times to count process using native implementation equivalent
 */
fun mtrace_iat2counts(
    interArrivalTimes: DoubleArray,
    timeScale: Double
): IntArray {
    
    val totalTime = interArrivalTimes.sum()
    val numWindows = Math.max(1, (totalTime / timeScale).toInt())
    val counts = IntArray(numWindows)
    
    var currentTime = 0.0
    var windowIndex = 0
    
    for (iat in interArrivalTimes) {
        currentTime += iat
        val targetWindow = Math.min((currentTime / timeScale).toInt(), numWindows - 1)
        
        // Increment count for the appropriate window
        for (w in windowIndex..targetWindow) {
            if (w < numWindows) {
                counts[w]++
                break // Only count once per arrival
            }
        }
        
        windowIndex = targetWindow
    }
    
    return counts
}
/**
 * Mtrace Count algorithms
 */
@Suppress("unused")
class MtraceCountAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}