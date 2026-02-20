/**
 * @file Markovian Arrival Process asymptotic index of dispersion analysis
 * 
 * Computes asymptotic index of dispersion for MAP counting processes, measuring long-term
 * variability and correlation effects. Essential for traffic characterization and buffer analysis.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the asymptotic index of dispersion (IDC) for a Markovian Arrival Process (MAP).
 *
 *
 * The asymptotic index of dispersion is defined as:
 * <pre>
 * I = SCV(1 + 2 * sum(\rho_k))
</pre> *
 * where SCV is the squared coefficient of variation, and \rho_k is the lag-k autocorrelation coefficient
 * of inter-arrival times. It represents the limiting value of the index of dispersion for counts (IDC) and
 * the index of dispersion for intervals (IDI).
 *
 * @param D0 the hidden transition matrix of the MAP
 * @param D1 the visible transition matrix of the MAP
 * @return the asymptotic index of dispersion
 */

fun map_idc(D0: Matrix, D1: Matrix): Double {
    val e = Matrix(D0.length(), 1, D0.length())
    for (i in 0..<D0.length()) {
        e[i, 0] = 1
    }

    return 1 + 2 * (map_lambda(D0, D1) - map_pie(D0, D1).mult(map_infgen(D0, D1).add(1.0, e.mult(map_piq(D0, D1)))
            .inv()).mult(D1).mult(e)[0])
}

/**
 * Computes the asymptotic index of dispersion (IDC) for a MAP
 * stored in a MatrixCell that contains the MAP's transition matrices.
 *
 * @param MAP a MatrixCell containing the transition matrices D0 and D1 of the MAP
 * @return the asymptotic index of dispersion
 */
fun map_idc(MAP: MatrixCell): Double {
    return map_idc(MAP[0], MAP[1])
}
/**
 * MAP idc algorithms
 */
@Suppress("unused")
class MapIdcAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}