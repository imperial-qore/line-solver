/**
 * @file Markovian Arrival Process autocorrelation function coefficients for counting processes
 * 
 * Computes ACFC values for MAP counting processes over time intervals, measuring temporal correlation
 * in event counts. Critical for analyzing bursty arrival patterns in network traffic modeling.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.Maths
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the autocorrelation function coefficients (ACFC) for a MAP counting process.
 *
 *
 * This method calculates the ACFC for a given MAP represented by matrices D0 and D1, at specified lags and over a specified timescale (u).
 * The ACFC provides insight into the temporal correlation of event counts over time.
 *
 * @param D0   the hidden transition matrix of the MAP
 * @param D1   the visible transition matrix of the MAP
 * @param lags an array of integers representing the lags at which to compute the ACFC
 * @param u    the length of the timeslot (timescale)
 * @return an array of doubles containing the ACFC values for the specified lags
 */

fun map_acfc(D0: Matrix, D1: Matrix, lags: IntArray, u: Double): DoubleArray {
    //        Example:
    //        MatrixCell MAP = mmpp2_fit(1,20,800,0.99);
    //        MAP.get(0).print();
    //        MAP.get(1).print();
    //
    //        int[] lags = {3, 4};
    //        double[] acfc = map_acfc(MAP.get(0), MAP.get(1), lags, 1.57);
    //        System.out.println(acfc[0]);
    //        System.out.println(acfc[1]);

    val n = D0.numCols
    val Q = map_infgen(D0, D1)
    val I = Matrix.eye(n)
    val piq = map_piq(D0, D1)
    val PRE = piq.mult(D1).mult(I.sub(Maths.matrixExp(Q.scale(u))))
    val e = Matrix.ones(n, 1)
    val inv2_epiqQ = (e.mult(piq).sub(Q)).square().inv()
    val POST = (I.sub(Maths.matrixExp(Q.scale(u)))).mult(inv2_epiqQ).mult(D1).mult(e)
    val vart = map_varcount(D0, D1, u)
    val acfCoeffs = DoubleArray(lags.size)
    for (i in lags.indices) {
        acfCoeffs[i] = PRE.mult(Maths.matrixExp(Q.scale((lags[i] - 1) * u))).mult(POST).scale(1 / vart).value()
    }
    return acfCoeffs
}

/**
 * Computes the autocorrelation function coefficients (ACFC) for a MAP counting process using a MatrixCell.
 *
 *
 * This method is a convenience overload that extracts the D0 and D1 matrices from the provided MatrixCell
 * and calculates the ACFC for the specified lags and timescale.
 *
 * @param MAP  the MatrixCell representing the MAP, containing the D0 and D1 matrices
 * @param lags an array of integers representing the lags at which to compute the ACFC
 * @param u    the length of the timeslot (timescale)
 * @return an array of doubles containing the ACFC values for the specified lags
 */

fun map_acfc(MAP: MatrixCell, lags: IntArray, u: Double): DoubleArray {
    return map_acfc(MAP[0], MAP[1], lags, u)
}

/**
 * MAP autocorrelation function coefficients algorithms
 */
@Suppress("unused")
class MapAcfcAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}