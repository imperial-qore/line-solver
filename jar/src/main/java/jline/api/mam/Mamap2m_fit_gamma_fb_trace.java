/**
 * @file Markovian Arrival MAP with Marked arrivals gamma forward-backward trace fitting
 *
 * Performs approximate fitting of a marked trace, yielding a second-order
 * acyclic MMAP that matches the class probabilities, the forward and
 * backward moments.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.api.trace.Mtrace_backward_moment;
import jline.api.trace.Mtrace_forward_moment;
import jline.api.trace.Mtrace_pc;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mamap2m_fit_gamma_fb_trace {
    private Mamap2m_fit_gamma_fb_trace() {}

    /**
     * Performs approximate fitting of a marked trace, yielding a second-order
     * acyclic MMAP that matches the class probabilities, the forward and
     * backward moments.
     *
     * @param T Inter-arrival times as array
     * @param A Class marks of each job as array
     * @return Fitted MAMAP[m]
     */
    public static MatrixCell mamap2m_fit_gamma_fb_trace(double[] T, int[] A) {
        // Compute moments
        double M1 = average(T);
        double M2 = averageSquared(T);
        double M3 = averageCubed(T);
        double GAMMA = trace_gamma(T);

        // Compute class characteristics
        double[] P = Mtrace_pc.mtrace_pc(T, A).toArray1D();
        double[] F = Mtrace_forward_moment.mtrace_forward_moment(T, A, new int[]{1}, 1).getColumn(0).toArray1D();
        double[] B = (double[]) Mtrace_backward_moment.mtrace_backward_moment(T, A, 1);

        return Mamap2m_fit_gamma_fb_mmap.mamap2m_fit_gamma_fb(M1, M2, M3, GAMMA, P, F, B);
    }

    /**
     * Performs approximate fitting of a marked trace from Matrix inputs.
     */
    public static MatrixCell mamap2m_fit_gamma_fb_trace(Matrix T, Matrix A) {
        double[] tArray = T.toArray1D();
        double[] aFlat = A.toArray1D();
        int[] aArray = new int[A.getNumRows() * A.getNumCols()];
        for (int i = 0; i < aArray.length; i++) {
            aArray[i] = (int) aFlat[i];
        }
        return mamap2m_fit_gamma_fb_trace(tArray, aArray);
    }

    /**
     * Computes autocorrelation decay rate (gamma) from trace.
     */
    public static double trace_gamma(double[] T) {
        if (T.length < 3) return 0.0;

        double mean = average(T);
        double variance = 0.0;
        for (double v : T) {
            variance += (v - mean) * (v - mean);
        }
        variance /= T.length;

        if (variance < 1e-10) return 0.0;

        // Compute lag-1 autocorrelation
        double cov1 = 0.0;
        for (int i = 0; i < T.length - 1; i++) {
            cov1 += (T[i] - mean) * (T[i + 1] - mean);
        }
        cov1 /= (T.length - 1);

        double rho1 = cov1 / variance;

        // Compute lag-2 autocorrelation
        double cov2 = 0.0;
        for (int i = 0; i < T.length - 2; i++) {
            cov2 += (T[i] - mean) * (T[i + 2] - mean);
        }
        cov2 /= (T.length - 2);

        double rho2 = cov2 / variance;

        if (Math.abs(rho1) > 1e-10) {
            return Math.max(-0.99, Math.min(0.99, rho2 / rho1));
        } else {
            return 0.0;
        }
    }

    private static double average(double[] arr) {
        double s = 0.0;
        for (double v : arr) s += v;
        return s / arr.length;
    }

    private static double averageSquared(double[] arr) {
        double s = 0.0;
        for (double v : arr) s += v * v;
        return s / arr.length;
    }

    private static double averageCubed(double[] arr) {
        double s = 0.0;
        for (double v : arr) s += v * v * v;
        return s / arr.length;
    }
}
