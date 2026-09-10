/**
 * @file Stationary mean fluid level of a level-dependent fluid queue
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.lib.butools.mam.LevelDependentFluidSolution;
import jline.lib.butools.mam.LevelDependentFluidStationary;

public final class Mfq_ld_mean {
    private Mfq_ld_mean() {}

    /** Stationary mean fluid level E[X]. */
    public static double mfq_ld_mean(LevelDependentFluidSolution sol, double[] T) {
        return LevelDependentFluidStationary.stationaryMean(sol, T);
    }
}
