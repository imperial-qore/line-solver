package jline.api.trace;

/**
 * Single trace mean computation.
 *
 * <p>Computes the arithmetic mean of empirical trace data. Fundamental statistical
 * measure used in trace analysis for characterizing central tendency and
 * parameterizing queueing models from measurement data.
 *
 * @since LINE 3.0
 */
public final class Trace_mean {
    private Trace_mean() {}

    /**
     * Computes the overall mean of the trace data.
     *
     * @param trace the array containing the trace data
     * @return the mean value of the trace data
     */
    public static double trace_mean(double[] trace) {
        double mean = 0.0;
        for (int i = 0; i < trace.length; i++) {
            mean += trace[i];
        }
        return mean / trace.length;
    }
}
