/**
 * @file Single trace mean computation
 * 
 * Computes the arithmetic mean of empirical trace data. Fundamental statistical 
 * measure used in trace analysis for characterizing central tendency and 
 * parameterizing queueing models from measurement data.
 * 
 * @since LINE 3.0
 */
package jline.api.trace

/**
 * Computes the overall mean of the trace data.
 *
 * @param trace the array containing the trace data
 * @return the mean value of the trace data
 */
fun trace_mean(trace: DoubleArray): Double {
    var mean = 0.0
    for (i in trace.indices) {
        mean += trace[i]
    }
    return mean / trace.size
}
/**
 * Trace Mean algorithms
 */
@Suppress("unused")
class TraceMeanAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}