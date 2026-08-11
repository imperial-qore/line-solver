package jline.api.mam;

import jline.io.Ret;

/**
 * Absorbing Phase-type distribution trace-based fitting.
 *
 * <p>Fits APH(2) distributions from empirical inter-arrival time traces.
 * Essential for data-driven phase-type distribution modeling from real measurements.
 *
 * @since LINE 3.0
 */
public final class Aph2_fit_trace {
    private Aph2_fit_trace() {}

    /**
     * Performs approximate fitting of a given trace, yielding a second-order
     * APH in canonical form.
     *
     * @param T The inter-arrival times
     * @return Fitted second-order phase-type distribution
     */
    public static Ret.mamAPH2Fit aph2_fit_trace(double[] T) {
        double sum1 = 0.0;
        double sum2 = 0.0;
        double sum3 = 0.0;
        int n = T.length;
        for (int i = 0; i < n; i++) {
            double t = T[i];
            double t2 = t * t;
            sum1 += t;
            sum2 += t2;
            sum3 += t2 * t;
        }
        double M1 = sum1 / n;
        double M2 = sum2 / n;
        double M3 = sum3 / n;

        return Aph2_fit.aph2_fit(M1, M2, M3);
    }
}
