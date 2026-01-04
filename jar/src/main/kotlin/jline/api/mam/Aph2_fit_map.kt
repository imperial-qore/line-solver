/**
 * @file Absorbing Phase-type distribution fitting from MAP
 * 
 * Fits APH(2) distributions by approximating arbitrary-order MAP processes.
 * Used for reducing MAP complexity while preserving key statistical properties.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.io.Ret
import jline.util.matrix.MatrixCell

/**
 * Performs approximate fitting of a MAP, yielding a second-order
 * APH in canonical form.
 *
 * @param map The MAP of arbitrary order to fit
 * @return Fitted second-order phase-type distribution
 */
fun aph2_fit_map(map: MatrixCell): Ret.mamAPH2Fit {
    val M1 = map_mean(map)
    val M2 = map_moment(map, 2)
    val M3 = map_moment(map, 3)
    
    return aph2_fit(M1, M2, M3)
}
/**
 * APH 2 fit map algorithms
 */
@Suppress("unused")
class Aph2FitMapAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}