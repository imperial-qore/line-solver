/**
 * @file Markovian Arrival Process Erlang-k distribution fitting
 * 
 * Constructs MAP representations of Erlang-k processes with specified means and phases.
 * Used for modeling low-variability arrival processes with coefficient of variation less than 1.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Fits an Erlang-k process as a Markovian Arrival Process (MAP).
 *
 *
 * The Erlang-k process is characterized by k phases, each with an exponential distribution, resulting in a distribution
 * with a mean and a coefficient of variation smaller than 1. This method constructs a MAP that approximates the behavior
 * of an Erlang-k process with a specified mean and shape parameter k.
 *
 * @param mean the desired mean of the Erlang-k process
 * @param k    the shape parameter, representing the number of phases in the Erlang-k process
 * @return a MatrixCell containing the transition matrices `D0` and `D1` of the fitted MAP
 */
fun map_erlang(mean: Double, k: Int): MatrixCell {
    val mu = k.toDouble() / mean
    var MAP = MatrixCell()
    val D0 = Matrix(k, k, 2 * k - 1)
    D0[0, 0] = -mu
    for (i in 0..<k - 1) {
        D0[i, i + 1] = mu
    }
    val D1 = Matrix(k, k, 1)
    D1[k - 1, 0] = mu
    MAP[0] = D0
    MAP[1] = D1
    MAP = map_normalize(MAP)
    return MAP
}
/**
 * MAP erlang algorithms
 */
@Suppress("unused")
class MapErlangAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}