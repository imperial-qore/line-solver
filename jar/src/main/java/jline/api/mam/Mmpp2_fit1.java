/**
 * @file Markov Modulated Poisson Process single-parameter fitting
 *
 * Fits MMPP(2) models using simplified single-parameter approach for specific scenarios.
 * Provides efficient fitting for cases with limited statistical information.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import org.apache.commons.math3.util.FastMath;

import jline.util.matrix.MatrixCell;

public final class Mmpp2_fit1 {
    private Mmpp2_fit1() {}

    /**
     * Fits a 2-phase Markov Modulated Poisson Process (MMPP2) based on the specified parameters.
     *
     * @param mean the mean inter-arrival time
     * @param scv  the squared coefficient of variation (SCV)
     * @param skew the skewness of the inter-arrival times
     * @param idc  the index of dispersion for counts (IDC)
     * @return a MatrixCell representing the fitted MMPP2 transition matrices
     */
    public static MatrixCell mmpp2_fit1(double mean, double scv, double skew, double idc) {
        double E1 = mean;
        double E2 = (1 + scv) * FastMath.pow(E1, 2);
        double g2 = -(scv - idc) / (-1 + idc);
        double E3;
        if (skew == -1.0) {
            E3 = -1.0;
        } else {
            E3 = -(2 * FastMath.pow(E1, 3) - 3 * E1 * E2 - skew * FastMath.pow(E2 - FastMath.pow(E1, 2), 1.5));
        }

        return Map2_fit.map2_fit(E1, E2, E3, g2).MAP;
    }
}
