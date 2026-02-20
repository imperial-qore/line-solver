/**
 * @file Markovian Arrival Process autocorrelation function analysis
 * 
 * Computes autocorrelation function (ACF) values for MAP inter-arrival times at specified lags.
 * Essential for analyzing temporal correlation patterns in arrival processes for queueing systems.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the autocorrelation function (ACF) for a given MAP at multiple lags.
 *
 *
 * This method calculates the ACF values for the specified lags by computing the embedded DTMC matrix,
 * the stationary distribution vector, and the steady-state vector. The ACF values are normalized by
 * the squared coefficient of variation (SCV) of the inter-arrival times.
 *
 * @param D0   the hidden transition matrix of the MAP
 * @param D1   the visible transition matrix of the MAP
 * @param lags a matrix containing the lags at which to compute the ACF
 * @return a matrix containing the ACF values at the specified lags
 */

fun map_acf(D0: Matrix, D1: Matrix, lags: Matrix): Matrix {
    val P = map_embedded(D0, D1)
    val x = map_piq(D0, D1)
    x.scaleEq(map_lambda(D0, D1))

    val neg_D0 = D0.copy()
    neg_D0.scaleEq(-1.0)
    val y = neg_D0.inv().sumRows()
    val acfCoeffs = Matrix(1, lags.length(), lags.length())
    for (i in 0..<lags.length()) {
        acfCoeffs[i] = x.mult(Matrix.pow(P, lags[i].toInt())).mult(y)[0]
    }
    for (i in 0..<acfCoeffs.length()) {
        acfCoeffs[i] = acfCoeffs[i] - 1
    }

    acfCoeffs.scaleEq(1 / map_scv(D0, D1))

    return acfCoeffs
}

/**
 * Computes the autocorrelation function (ACF) for a given MAP at multiple lags using a MatrixCell.
 *
 *
 * This method is a convenience overload that extracts the D0 and D1 matrices from the provided MatrixCell
 * and calculates the ACF for the specified lags.
 *
 * @param MAP  the MatrixCell representing the MAP, containing the D0 and D1 matrices
 * @param lags a matrix containing the lags at which to compute the ACF
 * @return a matrix containing the ACF values at the specified lags
 */
fun map_acf(MAP: MatrixCell, lags: Matrix): Matrix {
    return map_acf(MAP[0], MAP[1], lags)
}

/**
 * Computes the autocorrelation function (ACF) for a given MAP using a default lag of 1.
 *
 *
 * This method calculates the ACF at lag 1 for the MAP described by the matrices D0 and D1.
 * It is a convenience method for cases where only the first lag's ACF is of interest.
 *
 * @param D0 the hidden transition matrix of the MAP
 * @param D1 the visible transition matrix of the MAP
 * @return a matrix containing the ACF value at lag 1
 */
fun map_acf(D0: Matrix, D1: Matrix): Matrix {
    val LAGS = Matrix(1, 1, 1)
    LAGS[0, 0] = 1
    return map_acf(D0, D1, LAGS)
}

/**
 * Computes the autocorrelation function (ACF) for a given MAP using a default lag of 1.
 *
 *
 * This method is a convenience overload that extracts the D0 and D1 matrices from the provided MatrixCell
 * and calculates the ACF at lag 1.
 *
 * @param MAP the MatrixCell representing the MAP, containing the D0 and D1 matrices
 * @return a matrix containing the ACF value at lag 1
 */
fun map_acf(MAP: MatrixCell): Matrix {
    return map_acf(MAP[0], MAP[1])
}

/**
 * Computes the autocorrelation function (ACF) for a given MAP at a specific lag.
 *
 *
 * This method calculates the ACF for the specified lag by first creating a singleton matrix with the lag value.
 * It then calls the overloaded `map_acf` method that accepts a matrix of lags.
 *
 * @param D0  the hidden transition matrix of the MAP
 * @param D1  the visible transition matrix of the MAP
 * @param lag the specific lag at which to compute the ACF
 * @return a matrix containing the ACF value at the specified lag
 */
fun map_acf(D0: Matrix, D1: Matrix, lag: Int): Matrix {
    return map_acf(D0, D1, Matrix.singleton(lag.toDouble()))
}

/**
 * MAP autocorrelation function algorithms
 */
@Suppress("unused")
class MapAcfAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
