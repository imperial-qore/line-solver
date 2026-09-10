/**
 * @file Multi-class trace mean computation
 *
 * @since LINE 3.0
 */
package jline.api.trace;

import jline.util.matrix.Matrix;

public final class Mtrace_mean {
    private Mtrace_mean() {}

    /**
     * Computes the mean of a trace, divided by types.
     *
     * @param trace  the array containing the trace data
     * @param ntypes the number of different types
     * @param type   an array indicating the type of each element in the trace
     * @return a matrix containing the mean values for each type
     */
    public static Matrix mtrace_mean(double[] trace, int ntypes, int[] type) {
        double[] mean = new double[ntypes];
        for (int c = 0; c < ntypes; c++) {
            double sum = 0.0;
            int ctr = 0;
            for (int i = 0; i < trace.length; i++) {
                if (type[i] == c) {
                    sum += trace[i];
                    ctr++;
                }
            }
            mean[c] = ctr > 0 ? sum / (double) ctr : Double.NaN;
        }
        return new Matrix(mean).transpose();
    }
}
