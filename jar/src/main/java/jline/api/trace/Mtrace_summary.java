/**
 * @file Multi-trace summary statistics
 *
 * @since LINE 3.0
 */
package jline.api.trace;

import jline.util.matrix.Matrix;
import jline.lib.kpctoolbox.trace.TraceAnalysis;

public final class Mtrace_summary {
    private Mtrace_summary() {}

    /**
     * Data holder representing a summary of multi-trace statistics.
     */
    public static final class MtraceSummary {
        public final double[] M;       // Moments [M1, M2, M3, M4, M5]
        public final double[] ACF;     // Autocorrelation function (lags 1-100)
        public final Matrix F1;        // Forward moment of order 1
        public final Matrix F2;        // Forward moment of order 2
        public final Matrix B1;        // Backward moment of order 1
        public final Matrix B2;        // Backward moment of order 2
        public final Matrix C1;        // Cross moment of order 1
        public final Matrix C2;        // Cross moment of order 2
        public final Matrix Pc;        // Class probabilities
        public final Matrix Pab;       // Transition probabilities

        public MtraceSummary(double[] M, double[] ACF, Matrix F1, Matrix F2,
                             Matrix B1, Matrix B2, Matrix C1, Matrix C2,
                             Matrix Pc, Matrix Pab) {
            this.M = M;
            this.ACF = ACF;
            this.F1 = F1;
            this.F2 = F2;
            this.B1 = B1;
            this.B2 = B2;
            this.C1 = C1;
            this.C2 = C2;
            this.Pc = Pc;
            this.Pab = Pab;
        }
    }

    /**
     * Computes comprehensive summary statistics for a multi-class trace.
     *
     * @param T the inter-arrival times
     * @param C the class labels
     * @return summary object containing all statistics
     */
    public static MtraceSummary mtrace_summary(double[] T, int[] C) {
        // Compute basic moments
        double[] M = new double[5];
        for (int i = 1; i <= 5; i++) {
            double sum = 0.0;
            for (int k = 0; k < T.length; k++) {
                sum += Math.pow(T[k], i);
            }
            M[i - 1] = sum / T.length;
        }

        // Compute autocorrelation function
        int[] lags = new int[100];
        for (int k = 0; k < 100; k++) {
            lags[k] = k + 1;
        }
        double[] ACF = TraceAnalysis.trace_acf(T, lags);

        // Compute forward moments
        Matrix F1 = Mtrace_forward_moment.mtrace_forward_moment(T, C, new int[]{1});
        Matrix F2 = Mtrace_forward_moment.mtrace_forward_moment(T, C, new int[]{2});

        // Compute backward moments
        Matrix B1 = new Matrix((double[]) Mtrace_backward_moment.mtrace_backward_moment(T, C, 1));
        Matrix B2 = new Matrix((double[]) Mtrace_backward_moment.mtrace_backward_moment(T, C, 2));

        // Compute cross moments
        Matrix C1 = Mtrace_cross_moment.mtrace_cross_moment(T, C, 1);
        Matrix C2 = Mtrace_cross_moment.mtrace_cross_moment(T, C, 2);

        // Compute class probabilities
        Matrix Pc = Mtrace_pc.mtrace_pc(T, C);

        // Compute transition probabilities
        Matrix Pab = Mtrace_sigma.mtrace_sigma(T, C);

        return new MtraceSummary(M, ACF, F1, F2, B1, B2, C1, C2, Pc, Pab);
    }
}
