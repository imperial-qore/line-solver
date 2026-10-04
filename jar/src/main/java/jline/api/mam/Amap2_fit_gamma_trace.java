/**
 * @file Acyclic Markovian Arrival Process trace-based fitting with autocorrelation
 *
 * Fits AMAP(2) from empirical traces while preserving autocorrelation characteristics.
 * Essential for data-driven modeling of correlated arrival processes from measurements.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.List;

import jline.api.trace.Trace_var;
import jline.util.matrix.MatrixCell;
import jline.util.Pair;

public final class Amap2_fit_gamma_trace {
    private Amap2_fit_gamma_trace() {}

    /**
     * Performs approximate fitting of a given trace, yielding a second-order
     * AMAP in canonical form.
     *
     * @param T The inter-arrival times
     * @return Pair of (best AMAP, all feasible AMAPs)
     */
    public static Pair<MatrixCell, List<MatrixCell>> amap2_fit_gamma_trace(double[] T) {
        double sum = 0.0;
        double sum2 = 0.0;
        double sum3 = 0.0;
        for (int i = 0; i < T.length; i++) {
            double t = T[i];
            sum += t;
            sum2 += t * t;
            sum3 += t * t * t;
        }
        double n = T.length;
        double M1 = sum / n;
        double M2 = sum2 / n;
        double M3 = sum3 / n;
        double GAMMA = Trace_var.trace_gamma(T)[0];

        return Amap2_fit_gamma.amap2_fit_gamma(M1, M2, M3, GAMMA);
    }
}
