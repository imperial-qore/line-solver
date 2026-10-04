/**
 * @file Acyclic Markovian Arrival Process fitting from MAP with autocorrelation
 *
 * Fits AMAP(2) by approximating arbitrary-order MAP with preserved correlation structure.
 * Used for reducing MAP complexity while maintaining temporal correlation patterns.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.List;

import jline.util.Pair;
import jline.util.matrix.MatrixCell;

public final class Amap2_fit_gamma_map {
    private Amap2_fit_gamma_map() {}

    /**
     * Performs approximate fitting of a given MAP, yielding a second-order
     * AMAP in canonical form.
     *
     * @param map The MAP (of arbitrary order) to fit
     * @return Pair of (best AMAP, all feasible AMAPs)
     */
    public static Pair<MatrixCell, List<MatrixCell>> amap2_fit_gamma_map(MatrixCell map) {
        double M1 = Map_mean.map_mean(map);
        double M2 = Map_moment.map_moment(map, 2);
        double M3 = Map_moment.map_moment(map, 3);
        double GAMMA = Map_gamma.map_gamma(map);

        return Amap2_fit_gamma.amap2_fit_gamma(M1, M2, M3, GAMMA);
    }
}
