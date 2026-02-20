package jline.api.trace

import org.apache.commons.math3.stat.descriptive.moment.Skewness

/**
 * Computes the skewness of the trace data using Apache Commons Math.
 * This uses the same bias as MATLAB when calling skewness(data,0).
 *
 * @param trace the array containing the trace data
 * @return the skewness value of the trace data
 */
fun trace_skew(trace: DoubleArray): Double {
    val skewness = Skewness()
    return skewness.evaluate(trace)
}

/**
 * Computes the skewness of the trace data using Apache Commons Math.
 * This uses the same bias as MATLAB when calling skewness(data,0).
 *
 * @param trace the collection containing the trace data
 * @return the skewness value of the trace data
 */
fun trace_skew(trace: Collection<Double>): Double {
    return trace_skew(trace.toDoubleArray())
}

/**
 * Computes the skewness of the trace data using Apache Commons Math.
 * This uses the same bias as MATLAB when calling skewness(data,0).
 *
 * @param trace the list containing the trace data
 * @return the skewness value of the trace data
 */
fun trace_skew(trace: List<Double>): Double {
    return trace_skew(trace.toDoubleArray())
}
