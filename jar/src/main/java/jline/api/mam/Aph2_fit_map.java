package jline.api.mam;

import jline.io.Ret;
import jline.util.matrix.MatrixCell;

/**
 * Absorbing Phase-type distribution fitting from MAP.
 *
 * <p>Fits APH(2) distributions by approximating arbitrary-order MAP processes.
 * Used for reducing MAP complexity while preserving key statistical properties.
 *
 * @since LINE 3.0
 */
public final class Aph2_fit_map {
    private Aph2_fit_map() {}

    /**
     * Performs approximate fitting of a MAP, yielding a second-order
     * APH in canonical form.
     *
     * @param map The MAP of arbitrary order to fit
     * @return Fitted second-order phase-type distribution
     */
    public static Ret.mamAPH2Fit aph2_fit_map(MatrixCell map) {
        double M1 = Map_mean.map_mean(map);
        double M2 = Map_moment.map_moment(map, 2);
        double M3 = Map_moment.map_moment(map, 3);

        return Aph2_fit.aph2_fit(M1, M2, M3);
    }
}
