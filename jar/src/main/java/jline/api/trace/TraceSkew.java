package jline.api.trace;

import java.util.Collection;
import java.util.List;

import org.apache.commons.math3.stat.descriptive.moment.Skewness;

public final class TraceSkew {
    private TraceSkew() {}

    /**
     * Computes the skewness of the trace data using Apache Commons Math.
     * This uses the same bias as MATLAB when calling skewness(data,0).
     *
     * @param trace the array containing the trace data
     * @return the skewness value of the trace data
     */
    public static double trace_skew(double[] trace) {
        Skewness skewness = new Skewness();
        return skewness.evaluate(trace);
    }

    /**
     * Computes the skewness of the trace data using Apache Commons Math.
     * This uses the same bias as MATLAB when calling skewness(data,0).
     *
     * @param trace the collection containing the trace data
     * @return the skewness value of the trace data
     */
    public static double trace_skew(Collection<Double> trace) {
        double[] arr = new double[trace.size()];
        int i = 0;
        for (Double v : trace) {
            arr[i++] = v;
        }
        return trace_skew(arr);
    }

    /**
     * Computes the skewness of the trace data using Apache Commons Math.
     * This uses the same bias as MATLAB when calling skewness(data,0).
     *
     * @param trace the list containing the trace data
     * @return the skewness value of the trace data
     */
    public static double trace_skew(List<Double> trace) {
        double[] arr = new double[trace.size()];
        for (int i = 0; i < trace.size(); i++) {
            arr[i] = trace.get(i);
        }
        return trace_skew(arr);
    }
}
