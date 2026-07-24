/**
 * @file Markovian Arrival Process counting-process index of dispersion (IDC)
 *
 * Computes the time-dependent index of dispersion for counts (IDC) of a MAP,
 * I_a(t) = Var(A(t)) / E[A(t)], used by the Robust Queueing Network Analyzer.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.MatrixCell;

public final class Map_count_idc {
    private Map_count_idc() {}

    /**
     * Index of dispersion for counts (IDC) of a MAP at a set of time points.
     *
     * The IDC of the counting process A(t) associated to the MAP is
     * I_a(t) = Var(A(t)) / E[A(t)], t &gt; 0, i.e. the scaled variance-time
     * curve. It interpolates between I_a(0+) = SCV of the interarrival time
     * (renewal MAP) and the asymptotic value I_a(Inf) = map_idc(MAP).
     *
     * Reference: W. Whitt and W. You, "A Robust Queueing Network Analyzer Based
     * on Indices of Dispersion", eq. (1).
     *
     * @param MAP the Markovian Arrival Process stored in a MatrixCell
     * @param t   array of interval lengths (t &gt; 0)
     * @return array of IDC values, one per element of t
     */
    public static double[] map_count_idc(MatrixCell MAP, double[] t) {
        double[] m = Map_count_mean.map_count_mean(MAP, t);
        double[] v = Map_count_var.map_count_var(MAP, t);
        double[] I = new double[t.length];
        for (int i = 0; i < t.length; i++) {
            // For an orderly point process the counting IDC tends to 1 as t->0
            // (locally Poisson). Only hit at t==0 exactly.
            I[i] = (m[i] > 0) ? v[i] / m[i] : 1.0;
        }
        return I;
    }

    /**
     * Index of dispersion for counts (IDC) of a MAP at a single time point.
     *
     * @param MAP the Markovian Arrival Process stored in a MatrixCell
     * @param t   the interval length (t &gt; 0)
     * @return the IDC value at time t
     */
    public static double map_count_idc(MatrixCell MAP, double t) {
        return map_count_idc(MAP, new double[]{t})[0];
    }
}
