/**
 * @file Multi-class trace bootstrap resampling
 * 
 * Implements bootstrap resampling methods for multi-class empirical trace data. 
 * Provides statistical confidence intervals and uncertainty quantification for 
 * trace-derived parameters in queueing model fitting and validation.
 * 
 * @since LINE 3.0
 */
package jline.api.trace

import kotlin.random.Random

/**
 * Performs bootstrap resampling on a multi-class trace to estimate confidence intervals
 * and assess variability of trace statistics.
 *
 * @param interArrivalTimes Array of inter-arrival times
 * @param classLabels Array of class labels
 * @param numBootstraps Number of bootstrap samples (default: 1000)
 * @param blockSize Size of blocks for block bootstrap (default: auto-determined)
 * @param seed Random seed for reproducibility
 * @return Bootstrap results containing confidence intervals for key statistics
 */
fun mtrace_bootstrap(
    interArrivalTimes: DoubleArray,
    classLabels: IntArray,
    numBootstraps: Int = 1000,
    blockSize: Int? = null,
    seed: Long? = null
): BootstrapResults {
    
    require(interArrivalTimes.size == classLabels.size) {
        "Inter-arrival times and class labels must have same length"
    }
    
    val random = if (seed != null) Random(seed) else Random.Default
    val actualBlockSize = blockSize ?: determineOptimalBlockSize(interArrivalTimes)
    
    // Original statistics
    val originalStats = computeTraceStatistics(interArrivalTimes, classLabels)
    
    // Bootstrap samples
    val bootstrapStats = mutableListOf<TraceStatistics>()
    
    repeat(numBootstraps) {
        val (bootIAT, bootLabels) = if (actualBlockSize > 1) {
            blockBootstrapSample(interArrivalTimes, classLabels, actualBlockSize, random)
        } else {
            iidBootstrapSample(interArrivalTimes, classLabels, random)
        }
        
        val bootStats = computeTraceStatistics(bootIAT, bootLabels)
        bootstrapStats.add(bootStats)
    }
    
    // Compute confidence intervals
    val confidenceIntervals = computeConfidenceIntervals(bootstrapStats, 0.95)
    
    return BootstrapResults(
        original = originalStats,
        bootstrapSamples = bootstrapStats,
        confidenceIntervals = confidenceIntervals,
        blockSize = actualBlockSize
    )
}

/**
 * Data classes for bootstrap results
 */
data class BootstrapResults(
    val original: TraceStatistics,
    val bootstrapSamples: List<TraceStatistics>,
    val confidenceIntervals: Map<String, Pair<Double, Double>>,
    val blockSize: Int
)

data class TraceStatistics(
    val arrivalRate: Double,
    val scv: Double,
    val skewness: Double,
    val classProportions: DoubleArray,
    val lagCorrelations: DoubleArray
)

/**
 * Perform block bootstrap sampling
 */
private fun blockBootstrapSample(
    interArrivalTimes: DoubleArray,
    classLabels: IntArray,
    blockSize: Int,
    random: Random
): Pair<DoubleArray, IntArray> {
    
    val n = interArrivalTimes.size
    val numBlocks = (n + blockSize - 1) / blockSize
    
    val bootIAT = mutableListOf<Double>()
    val bootLabels = mutableListOf<Int>()
    
    // Resample blocks with replacement
    repeat(numBlocks) {
        val startIdx = random.nextInt(n - blockSize + 1)
        val endIdx = minOf(startIdx + blockSize, n)
        
        for (i in startIdx until endIdx) {
            bootIAT.add(interArrivalTimes[i])
            bootLabels.add(classLabels[i])
        }
    }
    
    // Truncate to original length
    val targetSize = minOf(bootIAT.size, n)
    return Pair(
        bootIAT.take(targetSize).toDoubleArray(),
        bootLabels.take(targetSize).toIntArray()
    )
}

/**
 * Perform IID bootstrap sampling
 */
private fun iidBootstrapSample(
    interArrivalTimes: DoubleArray,
    classLabels: IntArray,
    random: Random
): Pair<DoubleArray, IntArray> {
    
    val n = interArrivalTimes.size
    val bootIAT = DoubleArray(n)
    val bootLabels = IntArray(n)
    
    repeat(n) { i ->
        val idx = random.nextInt(n)
        bootIAT[i] = interArrivalTimes[idx]
        bootLabels[i] = classLabels[idx]
    }
    
    return Pair(bootIAT, bootLabels)
}

/**
 * Determine optimal block size using automatic bandwidth selection
 */
private fun determineOptimalBlockSize(interArrivalTimes: DoubleArray): Int {
    val n = interArrivalTimes.size
    
    // Compute autocorrelations
    val maxLag = minOf(50, n / 4)
    val autocorrs = DoubleArray(maxLag)
    
    val mean = interArrivalTimes.average()
    val variance = interArrivalTimes.map { (it - mean) * (it - mean) }.average()
    
    for (lag in 1 until maxLag) {
        var sumProduct = 0.0
        var count = 0
        
        for (i in 0 until n - lag) {
            sumProduct += (interArrivalTimes[i] - mean) * (interArrivalTimes[i + lag] - mean)
            count++
        }
        
        autocorrs[lag] = if (count > 0 && variance > 0) {
            sumProduct / (count * variance)
        } else {
            0.0
        }
    }
    
    // Find where autocorrelation becomes negligible
    var blockSize = 1
    for (lag in 1 until maxLag) {
        if (Math.abs(autocorrs[lag]) < 0.1) {
            blockSize = lag
            break
        }
    }
    
    return maxOf(1, minOf(blockSize, n / 10))
}

/**
 * Compute trace statistics
 */
private fun computeTraceStatistics(
    interArrivalTimes: DoubleArray,
    classLabels: IntArray
): TraceStatistics {
    
    val n = interArrivalTimes.size
    val numClasses = (classLabels.maxOrNull() ?: 0) + 1
    
    // Basic statistics
    val mean = interArrivalTimes.average()
    val variance = interArrivalTimes.map { (it - mean) * (it - mean) }.average()
    val arrivalRate = 1.0 / mean
    val scv = variance / (mean * mean)
    
    // Skewness
    val thirdMoment = interArrivalTimes.map { Math.pow(it - mean, 3.0) }.average()
    val skewness = if (variance > 0) {
        thirdMoment / Math.pow(variance, 1.5)
    } else {
        0.0
    }
    
    // Class proportions
    val classCounts = IntArray(numClasses)
    for (label in classLabels) {
        if (label >= 0 && label < numClasses) {
            classCounts[label]++
        }
    }
    val classProportions = classCounts.map { it.toDouble() / n }.toDoubleArray()
    
    // Lag correlations
    val maxLag = minOf(10, n / 4)
    val lagCorrelations = DoubleArray(maxLag)
    
    for (lag in 1 until maxLag) {
        var sumProduct = 0.0
        var count = 0
        
        for (i in 0 until n - lag) {
            sumProduct += (interArrivalTimes[i] - mean) * (interArrivalTimes[i + lag] - mean)
            count++
        }
        
        lagCorrelations[lag] = if (count > 0 && variance > 0) {
            sumProduct / (count * variance)
        } else {
            0.0
        }
    }
    
    return TraceStatistics(
        arrivalRate = arrivalRate,
        scv = scv,
        skewness = skewness,
        classProportions = classProportions,
        lagCorrelations = lagCorrelations
    )
}

/**
 * Compute confidence intervals for bootstrap samples
 */
private fun computeConfidenceIntervals(
    bootstrapStats: List<TraceStatistics>,
    confidence: Double
): Map<String, Pair<Double, Double>> {
    
    val alpha = 1.0 - confidence
    val lowerQuantile = alpha / 2.0
    val upperQuantile = 1.0 - alpha / 2.0
    
    val intervals = mutableMapOf<String, Pair<Double, Double>>()
    
    // Arrival rate
    val rates = bootstrapStats.map { it.arrivalRate }.sorted()
    intervals["arrivalRate"] = Pair(
        percentile(rates, lowerQuantile),
        percentile(rates, upperQuantile)
    )
    
    // SCV
    val scvs = bootstrapStats.map { it.scv }.sorted()
    intervals["scv"] = Pair(
        percentile(scvs, lowerQuantile),
        percentile(scvs, upperQuantile)
    )
    
    // Skewness
    val skews = bootstrapStats.map { it.skewness }.sorted()
    intervals["skewness"] = Pair(
        percentile(skews, lowerQuantile),
        percentile(skews, upperQuantile)
    )
    
    return intervals
}

/**
 * Compute percentile of sorted data
 */
private fun percentile(sortedData: List<Double>, p: Double): Double {
    if (sortedData.isEmpty()) return 0.0
    
    val n = sortedData.size
    val index = p * (n - 1)
    val lower = index.toInt()
    val upper = minOf(lower + 1, n - 1)
    val fraction = index - lower
    
    return sortedData[lower] + fraction * (sortedData[upper] - sortedData[lower])
}
/**
 * Mtrace Bootstrap algorithms
 */
@Suppress("unused")
class MtraceBootstrapAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}