/**
 * @file Multi-regime feedback fluid queue (Kankaya-Akar)
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.List;

import jline.lib.butools.mam.Multiregime;
import jline.util.matrix.Matrix;

public final class Mfq_multiregime {
    private Mfq_multiregime() {}

    /** Multi-regime feedback fluid queue; returns {pdf, pdfd, cdf, cdfm}. */
    public static Matrix[] mfq_multiregime(List<Matrix> Q, List<double[]> R, List<Matrix> Qt, List<double[]> Rt,
                                           double[] T, double[] pdfpoints, double[] cdfpoints) {
        return Multiregime.multiregime(Q, R, Qt, Rt, T, pdfpoints, cdfpoints);
    }
}
