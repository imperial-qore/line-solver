/**
 * @file MMAP mixture fitting from trace
 *
 * Fits a MMAP with m classes using a mixture of m^2 PH-distributions.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Mmap_mixture_fit_trace {
    private Mmap_mixture_fit_trace() {}

    /**
     * Fits a MMAP with m classes using a mixture of m^2 PH-distributions from trace data.
     *
     * Each PH distribution represents the probability distribution conditioned
     * on the fact that the last arrival was of class i and the next arrival is
     * of class j.
     *
     * @param T Inter-arrival times as array
     * @param A Class labels as array
     * @return Fitted MMAP and PH distributions for each transition
     */
    public static Ret.mamMMAPMixtureFit mmap_mixture_fit_trace(double[] T, int[] A) {
        // Compute two-step transition probabilities
        double[][][] sigma = jline.api.trace.Mtrace_sigma2.mtrace_sigma2(T, A);
        // Flatten 3D sigma to 2D matrix using marginal i->j (sum over h)
        int C = sigma.length;
        Matrix P2 = new Matrix(C, C);
        for (int i = 0; i < C; i++) {
            for (int j = 0; j < C; j++) {
                double s = 0.0;
                for (int h = 0; h < C; h++) {
                    s += sigma[i][j][h];
                }
                P2.set(i, j, s);
            }
        }

        // Compute cross moments of order 1, 2 and 3
        Matrix M1 = jline.api.trace.Mtrace_cross_moment.mtrace_cross_moment(T, A, 1);
        Matrix M2 = jline.api.trace.Mtrace_cross_moment.mtrace_cross_moment(T, A, 2);
        Matrix M3 = jline.api.trace.Mtrace_cross_moment.mtrace_cross_moment(T, A, 3);

        // Apply fitting algorithm
        return Mmap_mixture_fit.mmap_mixture_fit(P2, M1, M2, M3);
    }

    /**
     * Fits a MMAP with m classes using a mixture of m^2 PH-distributions from Matrix inputs.
     *
     * @param T Inter-arrival times as Matrix
     * @param A Class labels as Matrix
     * @return Fitted MMAP and PH distributions for each transition
     */
    public static Ret.mamMMAPMixtureFit mmap_mixture_fit_trace(Matrix T, Matrix A) {
        double[] tArray = T.toArray1D();
        int n = A.getNumRows() * A.getNumCols();
        double[] aArray1D = A.toArray1D();
        int[] aArray = new int[n];
        for (int i = 0; i < n; i++) {
            aArray[i] = (int) aArray1D[i];
        }
        return mmap_mixture_fit_trace(tArray, aArray);
    }
}
