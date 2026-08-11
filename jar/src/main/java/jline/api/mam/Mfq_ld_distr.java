/**
 * @file Stationary distribution of a level-dependent fluid queue
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.lib.butools.mam.LevelDependentFluidSolution;
import jline.lib.butools.mam.LevelDependentFluidStationary;
import jline.util.matrix.Matrix;

public final class Mfq_ld_distr {
    private Mfq_ld_distr() {}

    /** Stationary pdf/pdfd/cdf/cdfm at the requested points. */
    public static Matrix mfq_ld_distr(LevelDependentFluidSolution sol, double[] T, String what, double[] points) {
        return LevelDependentFluidStationary.stationaryDistr(sol, T, what, points);
    }
}
