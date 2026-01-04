/**
 * @file Acyclic Markovian Arrival Process fitting from MAP with autocorrelation
 * 
 * Fits AMAP(2) by approximating arbitrary-order MAP with preserved correlation structure.
 * Used for reducing MAP complexity while maintaining temporal correlation patterns.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.MatrixCell

/**
 * Performs approximate fitting of a given MAP, yielding a second-order
 * AMAP in canonical form.
 *
 * @param map The MAP (of arbitrary order) to fit
 * @return Pair of (best AMAP, all feasible AMAPs) 
 */
fun amap2_fit_gamma_map(map: MatrixCell): Pair<MatrixCell?, List<MatrixCell>> {
    val M1 = map_mean(map)
    val M2 = map_moment(map, 2)
    val M3 = map_moment(map, 3)
    val GAMMA = map_gamma(map)
    
    return amap2_fit_gamma(M1, M2, M3, GAMMA)
}
/**
 * Amap2 Fit Gamma Map algorithms
 */
@Suppress("unused")
class Amap2FitGammaMap {
    companion object {
        // Class documentation marker for Dokka
    }
}